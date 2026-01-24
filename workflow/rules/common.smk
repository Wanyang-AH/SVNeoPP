import csv

def load_samples(datatype):
    """
    返回包含{datatype}的字典列表，key为样本表中的列名。
    """
    samples = []
    with open(config["samples"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if (row.get("datatype") or "").strip().lower() == datatype:
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

def zip_fields(index, ids):
    """
    返回一个字典，包含给定样本ID列表对应的pair_id、datatype和sample_id字段。
    以wgs为例：
    {
    "pair_id":   ["l041", "l041", ...],
    "datatype":  ["wgs",  "wgs",  ...],
    "sample_id": ["l041_tumor_wgs_1", "l041_normal_wgs_1", ...],
    }
    """
    return dict(
        pair_id=[index[sid]["pair_id"] for sid in ids],
        datatype=[index[sid]["datatype"] for sid in ids],
        sample_id=ids,
    )


wgs_samples  = load_samples("wgs")
rna_samples  = load_samples("rna")
prot_samples = load_samples("prot")

wgs_index  = build_index(wgs_samples)
rna_index  = build_index(rna_samples)
prot_index = build_index(prot_samples)

WGS_IDS  = sorted(wgs_index.keys())
RNA_IDS  = sorted(rna_index.keys())
PROT_IDS = sorted(prot_index.keys())
