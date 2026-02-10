# SVNeoPP

SVNeoPP (Structural Variant Neoantigen Prediction and Prioritization) is a Snakemake workflow for structural-variant neoantigen analysis.
Current implemented chain includes:

- WGS/RNA preprocessing and alignment (`fastp`, `bwa`, `star`, `kallisto`)
- SV calling/annotation (`SvABA`, `AnnotSV`)
- Peptide generation from AnnotSV (`annotsv2pep`)
- netMHCpan input preparation, binding prediction, parse, and filter

## 1. Install Snakemake

Recommended (conda/mamba):

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake=8
conda activate snakemake
```

Optional check:

```bash
snakemake --version
```

## 2. Install Project from GitHub and Run Dry-Run

```bash
git clone git@github.com:Wanyang-AH/SVNeoPP.git
cd SVNeoPP
```

Full-project dry-run (requires your local data paths to be correctly configured):

```bash
snakemake -n --configfile config/config.yaml
```

Result-backtracking dry-run for `l041` (no full raw-data setup required):

This test starts from provided intermediate results and backtracks only the required upstream steps for netMHCpan.
Please ensure these files exist first:

- `results/annotsv2pep/l041/l041_AnnotSV.peptides.csv`
- `results/hla/optitype/l041/l041_tumor_rna_1_result.tsv`

Then run:

```bash
snakemake -n netmhcpan_filter_l041 --configfile config/config.yaml \
  --allowed-rules netchop_prepare_input netchop_run netchop_parse netchop_filter netmhcpan_prepare_input netmhcpan_run netmhcpan_parse netmhcpan_filter netmhcpan_filter_l041
```

## 3. Important Files and Locations

### `Snakefile`
- Path: `Snakefile`
- Role: entry point of the workflow, includes rules in `workflow/rules/*.smk`.

### `config/config.yaml`
- Path: `config/config.yaml`
- Role: central configuration for references, parameters, and output settings.
- Important:
- `references` paths must match real local software/resource locations.
- netMHCpan settings are under `params.netmhc`.

### `config/samples.tsv`
- Path: `config/samples.tsv`
- Role: sample sheet used by `workflow/rules/common.smk`.
- Important:
- Raw files must be placed at the exact paths listed in `r1`/`r2` columns.
- If `samples.tsv` paths are wrong, full workflow targets will fail.
- Required columns:
  - `pair_id`
  - `sample_id`
  - `replicate`
  - `status`
  - `datatype`
  - `r1`
  - `r2`
  - `platform`

### Rules and scripts
- Rules: `workflow/rules/*.smk`
- Python scripts: `scripts/*.py`
- netMHCpan-related scripts:
  - `scripts/netmhcpan/netmhcpan_prepare_input.py`
  - `scripts/netmhcpan/netmhcpan_parse.py`
  - `scripts/netmhcpan/netmhcpan_filter.py`

## 4. Resources Not Uploaded to GitHub and How to Configure

Large data/resources are expected under `resources/` and are configured in `config/config.yaml`.

### netMHCpan resource paths

- `references.netmhcpan_home`: `resources/netMHCpan-4.1/Linux_x86_64`
- `references.netmhcpan_bin`: `resources/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan`

### netMHCpan layout expectation

`references.netmhcpan_home` should point to netMHCpan runtime home (with `data/`), and `references.netmhcpan_bin` should be executable:

```text
resources/netMHCpan-4.1/
  Linux_x86_64/
    bin/netMHCpan
    data/
```

Workflow command internally uses:

```bash
NETMHCpan={netmhcpan_home} {netmhcpan_bin} -p {input_pep} -a {alleles} -BA -rdir {netmhcpan_home}
```

For full-pipeline runs, all other reference paths in `config/config.yaml` must also be configured.

## 5. Test Case Using Uploaded `results/` Data (L041)

This repository includes an L041 intermediate peptide file for fast module testing:

- `results/annotsv2pep/l041/l041_AnnotSV.peptides.csv`

You also need the L041 OptiType result:

- `results/hla/optitype/l041/l041_tumor_rna_1_result.tsv`

### Dry-run test (netMHCpan, l041)

```bash
snakemake -n results/netmhcpan/l041/l041.filtered.csv --configfile config/config.yaml
```

Or strict module-only dry-run:

```bash
snakemake -n netmhcpan_filter_l041 --configfile config/config.yaml \
  --allowed-rules netchop_prepare_input netchop_run netchop_parse netchop_filter netmhcpan_prepare_input netmhcpan_run netmhcpan_parse netmhcpan_filter netmhcpan_filter_l041
```

### Real run (single sample test)

```bash
snakemake -j 1 netmhcpan_filter_l041 --configfile config/config.yaml \
  --allowed-rules netchop_prepare_input netchop_run netchop_parse netchop_filter netmhcpan_prepare_input netmhcpan_run netmhcpan_parse netmhcpan_filter netmhcpan_filter_l041
```

### Expected outputs

- `results/netmhcpan/l041/l041.input.pep`
- `results/netmhcpan/l041/l041.input.map.csv`
- `results/netmhcpan/l041/l041.alleles.txt`
- `results/netmhcpan/l041/l041.raw.txt`
- `results/netmhcpan/l041/l041.xls`
- `results/netmhcpan/l041/l041.all.csv`
- `results/netmhcpan/l041/l041.filtered.csv`
- `logs/netmhcpan/l041.prepare.log`
- `logs/netmhcpan/l041.run.log`
- `logs/netmhcpan/l041.parse.log`
- `logs/netmhcpan/l041.filter.log`

## Notes

- Before committing, review `.gitignore` and update ignore rules according to your project data/output policy.
