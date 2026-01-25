import csv

def load_samples(datatype, status=None):
    """
    返回包含给定 datatype（以及可选 status）的字典列表。
    """
    samples = []
    datatype = (datatype or "").strip().lower()
    status = (status or "").strip().lower() or None
    with open(config["samples"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            row_datatype = (row.get("datatype") or "").strip().lower()
            row_status = (row.get("status") or "").strip().lower()
            if row_datatype != datatype:
                continue
            if status and row_status != status:
                continue
            samples.append(row)
    return samples

def build_index(rows):
    """
    返回每一行sample_id 为 key 的字典索引
    """
    idx = {}
    for row in rows:
        sid = row["sample_id"]
        if sid in idx:
            raise ValueError(f"Duplicate sample_id: {sid}")
        idx[sid] = row
    return idx


wgs_samples  = load_samples("wgs")
rna_samples  = load_samples("rna")
prot_samples = load_samples("prot")

wgs_index  = build_index(wgs_samples)
rna_index  = build_index(rna_samples)
prot_index = build_index(prot_samples)

WGS_IDS  = sorted(wgs_index.keys())
RNA_IDS  = sorted(rna_index.keys())
PROT_IDS = sorted(prot_index.keys())

# status-aware subsets（供需要 tumor/normal 过滤的规则使用）
tumor_wgs_samples = load_samples("wgs", status="tumor")
normal_wgs_samples = load_samples("wgs", status="normal")
tumor_rna_samples = load_samples("rna", status="tumor")
normal_rna_samples = load_samples("rna", status="normal")

tumor_wgs_index = build_index(tumor_wgs_samples)
normal_wgs_index = build_index(normal_wgs_samples)
tumor_rna_index = build_index(tumor_rna_samples)
normal_rna_index = build_index(normal_rna_samples)

TUMOR_WGS_IDS = sorted(tumor_wgs_index.keys())
NORMAL_WGS_IDS = sorted(normal_wgs_index.keys())
TUMOR_RNA_IDS = sorted(tumor_rna_index.keys())
NORMAL_RNA_IDS = sorted(normal_rna_index.keys())
