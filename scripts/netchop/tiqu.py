import pandas as pd
import argparse
import os
import sys
import re

"""
Script: merge_netchop_results_upgrade.py
Description: 
    Merges NetChop output with FASTA metadata based on sequence order.
    Calculates advanced metrics (N_score, C_score, max_internal, cleavage_count)
    to match the logic of the user's previous analysis pipeline.
"""

def parse_fasta_ordered(fasta_path):
    """
    按顺序解析 FASTA 文件，提取完整元数据。
    """
    data_list = []
    try:
        with open(fasta_path, 'r') as f:
            header = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    header = line[1:] 
                elif line and header:
                    # 解析 Header: >SV_ID|Gene|Source|Type
                    parts = header.split('|')
                    
                    sv_id = parts[0] if len(parts) > 0 else "NA"
                    gene = parts[1] if len(parts) > 1 else "NA"
                    source = parts[2] if len(parts) > 2 else "NA"
                    sv_type = parts[3] if len(parts) > 3 else "NA"
                    
                    data_list.append({
                        'Ident': sv_id,  # 使用 Ident 作为主键，方便和旧流程兼容
                        'Gene': gene,
                        'Mutation_Type': source, # 对应旧脚本的 Source_Tool
                        'sample': "NA", # 如果 Header 里没写 sample，后续可以从文件名提取，暂设NA
                        'sv_type': sv_type,
                        'peptide_sequence': line,
                        'peptide_length': len(line)
                    })
                    header = None 
    except Exception as e:
        print(f"[错误] 解析 FASTA 失败: {e}")
        sys.exit(1)
    return data_list

def parse_netchop_ordered(netchop_path):
    """
    按顺序解析 NetChop 输出块，计算每个肽段的 4 个关键指标。
    """
    results = []
    current_scores = [] # 存储 (pos, score, is_cleavage)
    
    try:
        with open(netchop_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # 忽略注释
                if not line or line.startswith('#') or "NetChop" in line:
                    continue
                
                # 检测分隔符 ---- (标志着一个肽段块的结束或开始)
                if line.startswith('----'):
                    if current_scores:
                        results.append(calculate_advanced_metrics(current_scores))
                        current_scores = [] 
                    continue
                
                # 检测数据行 (以数字开头)
                if re.match(r'^\d+\s+', line):
                    parts = line.split()
                    # NetChop 3.1 format: pos(0), AA(1), C(2), score(3), Ident(4)
                    if len(parts) >= 4:
                        try:
                            pos = int(parts[0])
                            is_cleavage = (parts[2] == 'S') # 'S' 表示 Cleavage Site
                            score = float(parts[3])
                            current_scores.append((pos, score, is_cleavage))
                        except ValueError:
                            pass

            # 处理最后一个块
            if current_scores:
                results.append(calculate_advanced_metrics(current_scores))
                
    except Exception as e:
        print(f"[错误] 解析 NetChop 失败: {e}")
        sys.exit(1)
    return results

def calculate_advanced_metrics(scores_list):
    """
    模拟旧脚本的逻辑计算指标：
    - N_score (Pos 1)
    - C_score (Pos End)
    - max_internal_score (Internal Max)
    - internal_cleavage_count (Internal 'S' count)
    """
    if not scores_list:
        return {'N_score': 0.0, 'C_score': 0.0, 'max_internal_score': 0.0, 'internal_cleavage_count': 0}
    
    # 1. N_score (第一个位置)
    n_score = scores_list[0][1] # score value
    
    # 2. C_score (最后一个位置)
    c_score = scores_list[-1][1]
    
    # 3. 内部指标 (剔除头尾)
    internal_items = scores_list[1:-1]
    
    if internal_items:
        internal_scores = [x[1] for x in internal_items]
        # 内部最大分数
        max_internal = max(internal_scores)
        # 内部切割计数 (旧脚本逻辑是看有多少个符合条件的行，这里我们用阈值0.5或标记S)
        # NetChop 默认阈值 0.5 对应 'S'
        internal_count = sum(1 for x in internal_items if x[2] is True) 
    else:
        max_internal = 0.0
        internal_count = 0
        
    return {
        'N_score': n_score,
        'C_score': c_score,
        'max_internal_score': max_internal,
        'internal_cleavage_count': internal_count
    }

def main():
    parser = argparse.ArgumentParser(description="Merge NetChop results (Sequence Order) and Calculate Risk Metrics.")
    parser.add_argument('--fasta', required=True, help="FASTA file with full metadata")
    parser.add_argument('--netchop', required=True, help="NetChop output file")
    parser.add_argument('--output', required=True, help="Output CSV file")
    
    args = parser.parse_args()
    
    # 1. 读取
    fasta_data = parse_fasta_ordered(args.fasta)
    netchop_data = parse_netchop_ordered(args.netchop)
    
    print(f"FASTA 序列数: {len(fasta_data)}")
    print(f"NetChop 结果块数: {len(netchop_data)}")
    
    # 2. 校验
    if len(fasta_data) != len(netchop_data):
        print(f"[警告] 数量不一致! 取最小值进行截断合并...")
        min_len = min(len(fasta_data), len(netchop_data))
        fasta_data = fasta_data[:min_len]
        netchop_data = netchop_data[:min_len]
        
    # 3. 合并
    merged_list = []
    for meta, metrics in zip(fasta_data, netchop_data):
        # 如果需要从文件名自动提取 sample (e.g., l056_AnnotSV...), 可以在这里做
        # 但为了简单，我们假设 FASTA Header 里不一定有 sample，先保持 NA 或从文件名传参
        item = {**meta, **metrics}
        merged_list.append(item)
        
    df = pd.DataFrame(merged_list)
    
    # 4. 格式化列名 (为了兼容旧的 summrry.py)
    # 旧脚本需要的列: N_score, C_score, max_internal_score, internal_cleavage_count
    # 以及 Metadata: Ident, peptide_sequence, peptide_length, sample, Gene, Source_Tool
    
    # 确保列存在且顺序合理
    desired_cols = [
        'Ident', 'peptide_sequence', 'peptide_length', 'Gene', 'Mutation_Type', 'sv_type',
        'N_score', 'C_score', 'max_internal_score', 'internal_cleavage_count'
    ]
    
    # 只保留存在的列
    final_cols = [c for c in desired_cols if c in df.columns]
    df = df[final_cols]
    
    # 5. 输出
    df.to_csv(args.output, index=False)
    print(f"结果已保存: {args.output}")

if __name__ == "__main__":
    main()