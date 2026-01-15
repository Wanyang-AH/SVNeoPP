import csv

def load_samples(datatype):
    samples = []
    with open(config["samples"], newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if (row.get("datatype") or "").strip().lower() == datatype:
                samples.append(row)
    return samples

def build_index(rows):
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
