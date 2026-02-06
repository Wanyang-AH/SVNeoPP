import pandas as pd
import pysam
import gffutils
from Bio.Seq import Seq
import os
import sys
import argparse
import warnings

# 忽略非关键警告
warnings.filterwarnings("ignore")

"""
SVNeo: Structural Variation Neoantigen Predictor (Final Version)
Version: 1.2
Logic:
  1. Transcript: AnnotSV ID > Mapping (RefSeq->Ensembl) > Longest CDS.
  2. Fusion: Direct concatenation (Preserve Frameshift).
  3. Filtering: 
     - Frameshift: Keep all peptides downstream of mutation.
     - In-frame: Keep ONLY peptides overlapping the mutation junction/region (Self-antigen filtering).
  4. Sequence Extraction: Uses flanking windows to capture junction-spanning peptides.
"""

# =============================================================================
# 1. Helper Functions (辅助工具)
# =============================================================================

def strip_version(ensembl_id):
    """移除版本号 (e.g., ENST000001.2 -> ENST000001)"""
    if isinstance(ensembl_id, str):
        return ensembl_id.split('.')[0]
    return ensembl_id

def normalize_chrom_name(chrom, fasta_refs):
    """
    标准化染色体名称以匹配 FASTA 文件 (处理 chr1 vs 1, MT vs chrM)
    """
    str_chrom = str(chrom)
    if str_chrom in fasta_refs: return str_chrom
    
    # 尝试加 chr
    if not str_chrom.startswith('chr'):
        alt = f"chr{str_chrom}"
        if alt in fasta_refs: return alt
        
    # 尝试去 chr
    if str_chrom.startswith('chr'):
        alt = str_chrom[3:]
        if alt in fasta_refs: return alt
        
    # 线粒体特殊处理
    if str_chrom in ['MT', 'M', 'chrM', 'chrMT']:
        for mt in ['MT', 'M', 'chrM', 'chrMT']:
            if mt in fasta_refs: return mt
            
    return None

def map_genomic_to_tx(pos, cds_exons):
    """
    将基因组坐标映射到转录本 CDS 的相对坐标
    """
    if not cds_exons: return None
    strand = cds_exons[0].strand
    
    # 按转录方向排序外显子
    sorted_exons = sorted(cds_exons, key=lambda x: x.start)
    if strand == '-':
        sorted_exons = sorted(cds_exons, key=lambda x: x.start, reverse=True)

    accumulated_len = 0
    
    # 越界处理
    if strand == '+' and pos < sorted_exons[0].start: return 0
    if strand == '-' and pos > sorted_exons[0].end: return 0

    for i, exon in enumerate(sorted_exons):
        # 坐标在外显子内
        if exon.start <= pos <= exon.end:
            offset = pos - exon.start if strand == '+' else exon.end - pos
            return accumulated_len + offset
        
        # 坐标在内含子内 (映射到上一个外显子末尾)
        if i + 1 < len(sorted_exons):
            next_exon = sorted_exons[i+1]
            in_intron = False
            if strand == '+':
                if exon.end < pos < next_exon.start: in_intron = True
            else:
                if next_exon.end < pos < exon.start: in_intron = True
            
            if in_intron:
                return accumulated_len + len(exon)
        
        accumulated_len += len(exon)

    return accumulated_len


PEPTIDE_COLUMNS = [
    "sample_name",
    "peptide",
    "length",
    "sv_type",
    "peptide_source",
    "gene_name",
    "transcript_id",
    "sv_id",
    "genomic_location",
]


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def normalize_annotsv_dataframe(df):
    """Normalize AnnotSV column names to a stable internal schema."""
    normalized = {}
    for col in df.columns:
        key = col.strip().replace(" ", "_").replace("-", "_")
        normalized[col] = key
    df = df.rename(columns=normalized)

    required = [
        "AnnotSV_ID",
        "Annotation_mode",
        "SV_type",
        "Overlapped_CDS_length",
        "Gene_name",
        "SV_chrom",
        "SV_start",
        "SV_end",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required AnnotSV columns: {', '.join(missing)}")

    for col in ["SV_start", "SV_end", "Overlapped_CDS_length"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    return df

# =============================================================================
# 2. Resource Loading (资源加载)
# =============================================================================

def setup_database(db_path, gtf_path, log_file):
    """加载或创建 gffutils 数据库"""
    if not os.path.exists(db_path):
        log_file.write(f"Database '{db_path}' not found. Creating from GTF...\n")
        print(f"Building GTF database from {gtf_path} (this may take a while)...")
        try:
            gffutils.create_db(
                gtf_path, 
                dbfn=db_path, 
                force=True, 
                keep_order=True, 
                merge_strategy='merge', 
                id_spec={'gene': 'gene_id', 'transcript': 'transcript_id'}, 
                disable_infer_genes=True, 
                disable_infer_transcripts=True
            )
            log_file.write("Database created successfully.\n")
        except Exception as e:
            log_file.write(f"FATAL ERROR: Failed to create GTF database: {e}\n")
            sys.exit(1)
    return gffutils.FeatureDB(db_path)

def build_id_mappers(filepath, log_file):
    """解析 ID 映射文件 (RefSeq -> Ensembl)"""
    if not filepath:
        log_file.write("No ID mapping file provided. Skip RefSeq->Ensembl mapping.\n")
        return {}, {}
    if not os.path.exists(filepath):
        log_file.write(f"ID mapping file not found: {filepath}. Skip mapping.\n")
        return {}, {}

    log_file.write(f"Parsing ID mapping file: {filepath}...\n")
    try:
        df = pd.read_csv(
            filepath,
            sep=r"\s+",
            engine="python",
            on_bad_lines="skip",
            names=["Gene_stable_ID", "Transcript_stable_ID", "Gene_name", "RefSeq_mRNA_ID"],
            usecols=[0, 2, 4, 6],
            comment="#",
            quoting=3,
        )
    except Exception as e:
        log_file.write(f"Failed to read mapping file: {e}\n")
        return {}, {}
    
    mapping_refseq_to_ensembl = {}
    count = 0
    for _, row in df.iterrows():
        tx_id = strip_version(row['Transcript_stable_ID'])
        refseq_id = strip_version(row['RefSeq_mRNA_ID'])
        if pd.notna(refseq_id) and pd.notna(tx_id):
            mapping_refseq_to_ensembl[refseq_id] = tx_id
            count += 1
    log_file.write(f"Mapping loaded. {count} RefSeq->Ensembl pairs.\n")
    return {}, mapping_refseq_to_ensembl

# =============================================================================
# 3. Transcript Selection (转录本选择)
# =============================================================================

def find_best_transcript_fallback(db, sv_row, log_file, ref_genome_refs):
    """回退策略：在区域内寻找 CDS 最长的转录本"""
    gene_name = sv_row['Gene_name']
    chrom = normalize_chrom_name(sv_row['SV_chrom'], ref_genome_refs)
    if not chrom: return None
    
    start, end = sv_row['SV_start'], sv_row['SV_end']
    log_file.write(f"  [Fallback] Searching for longest CDS transcript for {gene_name}...\n")
    
    best_tx = None
    max_len = -1
    
    # 尝试匹配 GTF 中的染色体名 (通常是不带 chr 或者带 chr)
    search_chroms = [chrom, chrom.replace('chr', '')]
    
    for c in search_chroms:
        try:
            found = list(db.region(region=(c, start, end), featuretype='transcript'))
            if not found: continue
            
            for tx in found:
                if gene_name in tx.attributes.get('gene_name', []):
                    cds_list = list(db.children(tx, featuretype='CDS'))
                    if cds_list:
                        cds_len = sum(len(c) for c in cds_list)
                        if cds_len > max_len:
                            max_len = cds_len
                            best_tx = tx
        except: continue
            
    if best_tx:
        log_file.write(f"  [Fallback Success] Selected {best_tx.id} (CDS: {max_len}bp)\n")
    else:
        log_file.write(f"  [Fallback Fail] No coding transcript found.\n")
        
    return best_tx

def get_transcript_info(sv_row, db, ref_genome, mappers, log_file):
    _, mapping_refseq_to_ensembl = mappers
    original_tx = sv_row.get('Tx', sv_row.get('transcript_id'))
    transcript = None
    
    # --- 1. AnnotSV ID ---
    if pd.notna(original_tx):
        query_id = original_tx
        # 映射 RefSeq (NM_/XM_) -> Ensembl
        if query_id.startswith(('NM_', 'XM_')):
            stripped = strip_version(query_id)
            if stripped in mapping_refseq_to_ensembl:
                query_id = mapping_refseq_to_ensembl[stripped]
                log_file.write(f"  [Mapping] {original_tx} -> {query_id}\n")
        
        try:
            transcript = db[query_id]
        except gffutils.exceptions.FeatureNotFoundError:
            try:
                transcript = db[strip_version(query_id)]
            except:
                pass 

    # --- 2. Fallback (Longest CDS) ---
    if not transcript:
        transcript = find_best_transcript_fallback(db, sv_row, log_file, ref_genome.references)
    
    if not transcript: return None, None, None
    
    # --- 提取序列 ---
    cds_exons = sorted([cds for cds in db.children(transcript, featuretype='CDS')], key=lambda x: x.start)
    if not cds_exons: return None, None, None
    
    wt_cds_seq = ""
    chrom = normalize_chrom_name(transcript.chrom, ref_genome.references)
    if not chrom: return None, None, None

    try:
        for cds in cds_exons: 
            wt_cds_seq += ref_genome.fetch(chrom, cds.start - 1, cds.end)
    except Exception as e:
        log_file.write(f"  [Error] Fetch seq failed: {e}\n")
        return None, None, None

    if transcript.strand == '-': 
        wt_cds_seq = str(Seq(wt_cds_seq).reverse_complement())
    
    return transcript, cds_exons, wt_cds_seq

# =============================================================================
# 4. Mutation Reconstruction (变异重构)
# =============================================================================

def reconstruct_mutant_mrna(sv_row, wt_cds_seq, cds_exons, log_file):
    sv_type = sv_row['SV_type']
    sv_start, sv_end = sv_row['SV_start'], sv_row['SV_end']

    tx_start = map_genomic_to_tx(sv_start, cds_exons)
    tx_end = map_genomic_to_tx(sv_end, cds_exons)

    # 内含子跳跃逻辑 (Exon Skipping/Dup)
    strand = cds_exons[0].strand
    affected_exons = [exon for exon in cds_exons if sv_start < exon.start and exon.end < sv_end]
    
    if affected_exons and sv_type in ['DEL', 'DUP', 'INV']:
        accumulated_len = 0
        sorted_exons = sorted(cds_exons, key=lambda x: x.start) if strand == '+' else sorted(cds_exons, key=lambda x: x.start, reverse=True)
        start_idx, end_idx, current_pos = -1, -1, 0
        for exon in sorted_exons:
            if exon in affected_exons:
                if start_idx == -1: start_idx = current_pos
                end_idx = current_pos + len(exon)
            current_pos += len(exon)
        
        if start_idx != -1:
            before, within, after = wt_cds_seq[:start_idx], wt_cds_seq[start_idx:end_idx], wt_cds_seq[end_idx:]
            if sv_type == 'DEL': return before + after
            elif sv_type == 'DUP': return before + within + within + after
            elif sv_type == 'INV': return before + str(Seq(within).reverse_complement()) + after

    # 内部断点逻辑
    if tx_start is None or tx_end is None: return None
    tx_start, tx_end = min(tx_start, tx_end), max(tx_start, tx_end)
    before, within, after = wt_cds_seq[:tx_start], wt_cds_seq[tx_start:tx_end], wt_cds_seq[tx_end:]
    
    if sv_type == 'DEL': return before + after
    elif sv_type == 'DUP': return before + within + within + after
    elif sv_type == 'INV': return before + str(Seq(within).reverse_complement()) + after
    elif sv_type == 'INS':
        ins = sv_row.get('SV_seq')
        return before + ins + after if (pd.notna(ins) and ins) else None
    return None

# =============================================================================
# 5. Peptide Generation (侧翼提取与过滤)
# =============================================================================

def generate_peptides_with_filtering(
    mt_protein,        # 完整的突变蛋白
    junction_idx,      # 变异起始位置 (索引)
    mutation_end_idx,  # 变异结束位置 (如果是 Frameshift 则为 None)
    sample_name, sv_row, pep_source, gene_name, tx_id, is_frameshift
):
    peptides = []
    # 仅生成 8-11 mer
    min_len, max_len = 8, 11
    
    # 1. 定义 "感兴趣区域" (ROI) 以减少计算量
    # 向前多取 (max_len - 1) 个氨基酸，确保能生成跨越接头的肽段
    flank_size = max_len - 1
    start_extract = max(0, junction_idx - flank_size)
    
    if is_frameshift or mutation_end_idx is None:
        # 移码：取到末尾 (限制一下最大长度以防万一，比如取接头后 100 个 AA，视需求而定)
        # 这里暂时取全部，如果蛋白太长(>2000)可能需要截断
        extract_seq = mt_protein[start_extract:]
    else:
        # In-frame：取到变异结束点后 flank_size
        end_extract = min(len(mt_protein), mutation_end_idx + flank_size)
        extract_seq = mt_protein[start_extract:end_extract]

    # 计算 ROI 内的相对坐标
    rel_junction_start = junction_idx - start_extract
    if mutation_end_idx is not None:
        rel_mutation_end = mutation_end_idx - start_extract
    else:
        rel_mutation_end = len(extract_seq) # Frameshift 视为一直变异到底

    # 2. 滑动窗口生成
    for length in range(min_len, max_len + 1):
        for i in range(len(extract_seq) - length + 1):
            pep_seq = extract_seq[i : i+length]
            
            # 肽段在 extract_seq 中的区间 [pep_start, pep_end)
            pep_start = i
            pep_end = i + length
            
            # --- [核心过滤: 自身抗原剔除] ---
            # 肽段必须与变异区域有 "物理重叠"
            # Overlap 条件: (肽段结束 > 变异开始) AND (肽段开始 < 变异结束)
            has_overlap = (pep_end > rel_junction_start) and (pep_start < rel_mutation_end)
            
            if not has_overlap:
                continue # 纯自身抗原，丢弃
                
            peptides.append({
                "sample_name": sample_name,
                "peptide": pep_seq,
                "length": length,
                "sv_type": sv_row['SV_type'],
                "peptide_source": pep_source, # Unified tag: 'frameshift' or 'X_in-frame'
                "gene_name": gene_name,
                "transcript_id": tx_id,
                "sv_id": sv_row['AnnotSV_ID'],
                "genomic_location": f"{sv_row['SV_chrom']}:{sv_row['SV_start']}-{sv_row.get('SV_end', '')}"
            })
    return peptides

# --- Process Single SV ---
def process_single_event(sv_row, db, ref_genome, mappers, log_file, tumor_col, normal_col, sample_name):
    log_file.write(f"\n--- SV: {sv_row['AnnotSV_ID']} ({sv_row['SV_type']}) Gene: {sv_row['Gene_name']} ")
    transcript, cds_exons, wt_cds_seq = get_transcript_info(sv_row, db, ref_genome, mappers, log_file)
    if not all([transcript, cds_exons, wt_cds_seq]): return None
    
    mt_cds_seq = reconstruct_mutant_mrna(sv_row, wt_cds_seq, cds_exons, log_file)
    if not mt_cds_seq or wt_cds_seq == mt_cds_seq: return None
    
    wt_pro = str(Seq(wt_cds_seq).translate(to_stop=True))
    mt_pro = str(Seq(mt_cds_seq).translate(to_stop=True))
    if wt_pro == mt_pro: return None
    
    # 找变异点
    j_idx = next((i for i, (a, b) in enumerate(zip(wt_pro, mt_pro)) if a != b), min(len(wt_pro), len(mt_pro)))
    
    # 判断 Frameshift
    is_fs = sv_row.get('Frameshift') == 'yes'
    if not is_fs and (len(mt_cds_seq) - len(wt_cds_seq)) % 3 != 0: is_fs = True
    
    # 确定变异结束点
    mutation_end_idx = None
    if not is_fs:
        # In-frame: 估算变异长度，至少为1 (接头)
        # 如果是 DEL，mt比wt短，我们认为接头点就是变异区 (长度1)
        # 如果是 INS/DUP，长度差就是变异区
        len_diff = max(1, len(mt_pro) - len(wt_pro))
        if sv_row['SV_type'] == 'DEL': len_diff = 1 
        mutation_end_idx = j_idx + len_diff

    log_file.write(f"  [Success] Type: {'Frameshift' if is_fs else 'In-frame'}\n")
    
    pep_source = "frameshift" if is_fs else f"{sv_row['SV_type']}_in-frame"
    
    return generate_peptides_with_filtering(
        mt_pro, j_idx, mutation_end_idx,
        sample_name, sv_row, pep_source, sv_row['Gene_name'], transcript.id, is_fs
    )

# --- Process Fusion ---
def process_fusion_event(sv_row, mate_row, db, ref_genome, mappers, log_file, tumor_col, normal_col, sample_name):
    gene1, gene2 = sv_row['Gene_name'], mate_row['Gene_name']
    log_file.write(f"\n--- Fusion: {gene1}::{gene2} ")
    
    tx1, cds1, seq1 = get_transcript_info(sv_row, db, ref_genome, mappers, log_file)
    tx2, cds2, seq2 = get_transcript_info(mate_row, db, ref_genome, mappers, log_file)
    if not all([tx1, cds1, seq1, tx2, cds2, seq2]): return None
    
    bp1 = map_genomic_to_tx(sv_row['SV_start'], cds1)
    bp2 = map_genomic_to_tx(mate_row['SV_start'], cds2)
    if bp1 is None or bp2 is None: return None
    
    head = seq1[:bp1]
    tail = seq2[bp2:]
    
    # 融合拼接 (保留移码)
    fusion_cds = head + tail
    mt_pro = str(Seq(fusion_cds).translate(to_stop=True))
    
    # 变异点: Head 的末尾
    j_idx = len(head) // 3
    
    # 判断 Frameshift
    is_fs = (len(head) % 3) != 0
    
    # 变异结束点
    mutation_end_idx = None
    if not is_fs:
        # In-frame Fusion: 变异区仅仅是接头点 (我们定义接头点宽度为1，确保只保留跨越它的肽段)
        mutation_end_idx = j_idx + 1
    
    log_file.write(f"  [Success] Type: {'Frameshift' if is_fs else 'In-frame'}\n")
    
    pep_source = "frameshift" if is_fs else "fusion_in-frame"

    return generate_peptides_with_filtering(
        mt_pro, j_idx, mutation_end_idx,
        sample_name, sv_row, pep_source, f"{gene1}-{gene2}", f"{tx1.id}-{tx2.id}", is_fs
    )

# =============================================================================
# 6. Main Entry
# =============================================================================

def process_single_pair(annotsv_path, output_csv_path, db, ref_genome, mappers, log_file, sample_name):
    log_file.write(f"\n{'='*20}\nProcessing Sample: {sample_name}\n{'='*20}\n")
    ensure_parent_dir(output_csv_path)

    try:
        df = pd.read_csv(annotsv_path, sep='\t', low_memory=False)
        df = normalize_annotsv_dataframe(df)
    except Exception as e:
        log_file.write(f"Failed to read/normalize TSV: {e}\n")
        pd.DataFrame(columns=PEPTIDE_COLUMNS).to_csv(output_csv_path, index=False)
        return

    sv_types = ['DEL', 'DUP', 'INV', 'INS', 'BND']
    candidates = df[
        (df['Annotation_mode'] == 'split') &
        (df['SV_type'].isin(sv_types)) &
        (df['Overlapped_CDS_length'] > 0) &
        (df['Gene_name'].notna())
    ].copy()

    log_file.write(f"Found {len(candidates)} CDS-overlapping events.\n")

    all_peptides = []
    processed_fusion_ids = set()
    t_col, n_col = None, None

    for _, row in candidates.iterrows():
        res = None
        sv_id, sv_type = row['AnnotSV_ID'], row['SV_type']

        if sv_type == 'BND' and sv_id not in processed_fusion_ids:
            base_id = sv_id.rsplit('_', 1)[0]
            mate_id = f"{base_id}_2" if sv_id.endswith('_1') else f"{base_id}_1"
            mate_rows = candidates[candidates['AnnotSV_ID'] == mate_id]
            if not mate_rows.empty:
                mate_row = mate_rows.iloc[0]
                if row['Gene_name'] != mate_row['Gene_name']:
                    res = process_fusion_event(row, mate_row, db, ref_genome, mappers, log_file, t_col, n_col, sample_name)
                processed_fusion_ids.add(sv_id); processed_fusion_ids.add(mate_id)
        elif sv_type in ['DEL', 'DUP', 'INV', 'INS']:
            res = process_single_event(row, db, ref_genome, mappers, log_file, t_col, n_col, sample_name)

        if res:
            all_peptides.extend(res)

    if all_peptides:
        out_df = pd.DataFrame(all_peptides).drop_duplicates(subset=['peptide'])
        out_df = out_df.reindex(columns=PEPTIDE_COLUMNS)
        log_file.write(f"Generated {len(out_df)} peptides.\n")
        print(f"Success: {len(out_df)} peptides -> {output_csv_path}")
    else:
        out_df = pd.DataFrame(columns=PEPTIDE_COLUMNS)
        log_file.write("No valid peptides generated.\n")
        print(f"Warning: No peptides found for {sample_name}")

    out_df.to_csv(output_csv_path, index=False)


def run_pipeline(annotsv_file, output_csv, fasta, gtf, id_map, db_path, log_path, sample_name):
    ensure_parent_dir(output_csv)
    ensure_parent_dir(log_path)

    with open(log_path, "w") as log_file:
        log_file.write(
            "Inputs:\n"
            f"  annotsv={annotsv_file}\n"
            f"  fasta={fasta}\n"
            f"  gtf={gtf}\n"
            f"  id_map={id_map}\n"
            f"  db={db_path}\n"
            f"  output={output_csv}\n"
        )
        try:
            print("Initializing resources...")
            db = setup_database(db_path, gtf, log_file)
            ref_genome = pysam.FastaFile(fasta)
            mappers = build_id_mappers(id_map, log_file)
            process_single_pair(annotsv_file, output_csv, db, ref_genome, mappers, log_file, sample_name)
        except Exception:
            import traceback
            log_file.write(f"\nFATAL ERROR:\n{traceback.format_exc()}\n")
            raise

def main():
    parser = argparse.ArgumentParser(
        description="SVNeo: Structural Variation Neoantigen Predictor",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--file', required=True, help="Input AnnotSV .tsv file")
    parser.add_argument('-o', '--output_dir', default="./output", help="Output directory")
    parser.add_argument('--output_csv', default="", help="Output CSV file path")
    parser.add_argument('-g', '--fasta', required=True, help="Reference Genome FASTA")
    parser.add_argument('-t', '--gtf', required=True, help="Gene Annotation GTF")
    parser.add_argument('-m', '--id_map', default="", help="ID Mapping file (optional)")
    parser.add_argument('-d', '--db', required=True, help="gffutils Database File (.db)")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    base_name = os.path.basename(args.file).split('.')[0]
    output_csv_path = args.output_csv or os.path.join(args.output_dir, f"{base_name}.peptides.csv")
    log_path = os.path.join(os.path.dirname(output_csv_path), "svneo.log")
    run_pipeline(
        annotsv_file=args.file,
        output_csv=output_csv_path,
        fasta=args.fasta,
        gtf=args.gtf,
        id_map=args.id_map,
        db_path=args.db,
        log_path=log_path,
        sample_name=base_name,
    )


def run_from_snakemake():
    annotsv_file = str(snakemake.input.tsv)
    output_csv = str(snakemake.output.csv)
    fasta = str(snakemake.input.fasta)
    gtf = str(snakemake.input.gtf)
    id_map = str(snakemake.input.id_map) if "id_map" in snakemake.input else ""
    db_path = str(snakemake.input.db)
    log_path = str(snakemake.log[0]) if snakemake.log else os.path.join(os.path.dirname(output_csv), "svneo.log")
    pair_id = getattr(snakemake.wildcards, "pair_id", None)
    sample_name = str(pair_id) if pair_id else os.path.basename(annotsv_file).split(".")[0]
    run_pipeline(
        annotsv_file=annotsv_file,
        output_csv=output_csv,
        fasta=fasta,
        gtf=gtf,
        id_map=id_map,
        db_path=db_path,
        log_path=log_path,
        sample_name=sample_name,
    )

if __name__ == "__main__":
    if "snakemake" in globals():
        run_from_snakemake()
    else:
        main()
