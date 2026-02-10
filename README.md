# SVNeoPP

SVNeoPP (**S**tructural **V**ariant **Neo**antigen **P**rediction and **P**rioritization) is a **Snakemake** workflow for **structural-variant–derived neoantigen** analysis, designed to be modular, reproducible, and easy to run from either raw data or intermediate results.

## What’s implemented (current chain)

The current implemented chain includes:

* WGS/RNA preprocessing and alignment (`fastp`, `bwa`, `STAR`, `kallisto`)
* SV calling/annotation (`SvABA`, `AnnotSV`)
* Peptide generation from AnnotSV annotations (`annotsv2pep`)
* HLA typing integration (e.g. OptiType outputs as allele input)
* Antigen processing prescreen (NetChop rules are included in the workflow targets)
* netMHCpan input preparation → binding prediction → parse → filter

Planned / TODO (not yet fully integrated):

* `MHCflurry` (binding/presentation prediction)
* `DeepImmuno` (immunogenicity analysis)
* `FragPipe` (mass-spectrometry database search)

> Note: This repository expects some **large resources / licensed tools** (e.g. netMHCpan) to be provided locally and configured via `config/config.yaml`.

---

## Quickstart (recommended: run the included L041 module test)

This is the fastest way to verify your environment and the netMHCpan-related module, without setting up full raw data paths.

### 0) Install Snakemake

Recommended via conda/mamba:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake=8
conda activate snakemake

snakemake --version
```

Snakemake official resources (for reference, not required to run):

* Documentation: `https://snakemake.readthedocs.io/en/stable/`
* GitHub: `https://github.com/snakemake/snakemake`

### 1) Clone the repository

```bash
git clone git@github.com:Wanyang-AH/SVNeoPP.git
cd SVNeoPP
```

### 2) Confirm the L041 test inputs exist

This repository provides an intermediate peptide file for fast testing, and expects an OptiType result for HLA alleles. Please ensure these files exist first:

* `results/annotsv2pep/l041/l041_AnnotSV.peptides.csv`
* `results/hla/optitype/l041/l041_tumor_rna_1_result.tsv`

### 3) Configure netMHCpan paths in `config/config.yaml`

Make sure you have netMHCpan installed locally (license-restricted; do not commit it). The workflow expects something like:

```text
resources/netMHCpan-4.1/
  Linux_x86_64/
    bin/netMHCpan
    data/
```

And in `config/config.yaml`:

* `references.netmhcpan_home`: `resources/netMHCpan-4.1/Linux_x86_64`
* `references.netmhcpan_bin`: `resources/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan`

Internally, the workflow uses a command pattern like:

```bash
NETMHCpan={netmhcpan_home} {netmhcpan_bin} -p {input_pep} -a {alleles} -BA -rdir {netmhcpan_home}
```

### 4) Run a dry-run of the module (L041)

Option A (simple target dry-run):

```bash
snakemake -n results/netmhcpan/l041/l041.filtered.csv --configfile config/config.yaml
```

Option B (strict module-only dry-run; backtrack only what’s needed for netMHCpan):

```bash
snakemake -n netmhcpan_filter_l041 --configfile config/config.yaml \
  --allowed-rules \
    netchop_prepare_input netchop_run netchop_parse netchop_filter \
    netmhcpan_prepare_input netmhcpan_run netmhcpan_parse netmhcpan_filter \
    netmhcpan_filter_l041
```

### 5) Real run (single-sample module test)

```bash
snakemake -j 1 netmhcpan_filter_l041 --configfile config/config.yaml \
  --allowed-rules \
    netchop_prepare_input netchop_run netchop_parse netchop_filter \
    netmhcpan_prepare_input netmhcpan_run netmhcpan_parse netmhcpan_filter \
    netmhcpan_filter_l041
```

---

## Running the full workflow (your own samples)

### 1) Fill in `config/samples.tsv`

* Path: `config/samples.tsv`
* Role: sample sheet used by `workflow/rules/common.smk`

Required columns (must exist):

* `pair_id`
* `sample_id`
* `replicate`
* `status`
* `datatype`
* `r1`
* `r2`
* `platform`

Important notes:

* Raw files must be placed at the exact paths listed in `r1` / `r2`.
* If the paths are wrong, full workflow targets will fail.

### 2) Configure `config/config.yaml`

* Path: `config/config.yaml`
* Role: central configuration for references, parameters, and output settings.

Key points:

* `references` paths must match real local software/resource locations.
* netMHCpan settings are under `params.netmhc`.

### 3) Dry-run the full workflow

```bash
snakemake -n --configfile config/config.yaml
```

### 4) Execute (example)

Pick an appropriate number of cores for your machine:

```bash
snakemake --cores 16 --configfile config/config.yaml
```

Recommended Snakemake flags you may consider (optional, depending on your environment):

* `--rerun-incomplete` (avoid partial outputs)
* `--keep-going` (continue other jobs if one fails)
* `--printshellcmds` (print commands for debugging)
* `--reason` (show why jobs are triggered)
* `--latency-wait 60` (help on network filesystems)

---

## Expected outputs (netMHCpan module, L041)

For the L041 module test, expected outputs include:

* `results/netmhcpan/l041/l041.input.pep`
  Peptide sequences used as input to netMHCpan.
* `results/netmhcpan/l041/l041.input.map.csv`
  Mapping table from peptide → source records (e.g., originating SV annotation / identifiers).
* `results/netmhcpan/l041/l041.alleles.txt`
  HLA alleles used for prediction (e.g., parsed from OptiType result).
* `results/netmhcpan/l041/l041.raw.txt`
  Raw text output from netMHCpan (stdout/combined output).
* `results/netmhcpan/l041/l041.xls`
  netMHCpan `.xls` output file.
* `results/netmhcpan/l041/l041.all.csv`
  Parsed/merged results in CSV form.
* `results/netmhcpan/l041/l041.filtered.csv`
  Filtered/prioritized results (workflow-defined thresholds).

Logs (one per step):

* `logs/netmhcpan/l041.prepare.log`
* `logs/netmhcpan/l041.run.log`
* `logs/netmhcpan/l041.parse.log`
* `logs/netmhcpan/l041.filter.log`

> Tip: When debugging, always read the corresponding `logs/netmhcpan/...*.log` first.

---

## Repository layout

Key locations:

* `Snakefile`
  Workflow entry point; includes rules in `workflow/rules/*.smk`.
* `workflow/rules/*.smk`
  Snakemake rule modules.
* `scripts/*.py`
  Helper scripts used by rules.
* netMHCpan-related scripts:

  * `scripts/netmhcpan/netmhcpan_prepare_input.py`
  * `scripts/netmhcpan/netmhcpan_parse.py`
  * `scripts/netmhcpan/netmhcpan_filter.py`
* `config/config.yaml`
  Main configuration.
* `config/samples.tsv`
  Sample sheet.
* `resources/`
  Large resources (not uploaded) and licensed tools.

---

## Resources not uploaded to GitHub

Large data/resources are expected under `resources/` and configured via `config/config.yaml`.

### netMHCpan

Required configuration (example):

* `references.netmhcpan_home`: `resources/netMHCpan-4.1/Linux_x86_64`
* `references.netmhcpan_bin`: `resources/netMHCpan-4.1/Linux_x86_64/bin/netMHCpan`

Expected layout:

```text
resources/netMHCpan-4.1/
  Linux_x86_64/
    bin/netMHCpan
    data/
```

---

## Reproducibility & visualization (optional but recommended)

If you have Graphviz installed, you can visualize the workflow:

```bash
snakemake --dag --configfile config/config.yaml | dot -Tpng > dag.png
snakemake --rulegraph --configfile config/config.yaml | dot -Tpng > rulegraph.png
```

Generate an HTML report (Snakemake built-in):

```bash
snakemake --report snakemake_report.html --configfile config/config.yaml
```

---

## Troubleshooting (common issues)

* **`snakemake -n` fails immediately**

  * Check `config/config.yaml` paths, especially `references.*` (must exist locally).
* **netMHCpan cannot run / cannot find data**

  * Verify `references.netmhcpan_home` points to the directory containing `data/`, and `references.netmhcpan_bin` is executable.
* **Full workflow fails but module test works**

  * Most commonly due to incorrect `config/samples.tsv` file paths (`r1`/`r2`).
* **Permission issues on `bin/netMHCpan`**

  * Ensure it is executable (`chmod +x .../bin/netMHCpan`).

---

## Roadmap

Planned modules (not yet integrated):

* `MHCflurry` for additional binding/presentation scoring
* `DeepImmuno` for immunogenicity scoring
* `FragPipe` for mass-spectrometry database search integration

---


## License

This project is licensed under the terms in the repository root `LICENSE` file.
Also note that some third-party tools/resources (e.g., netMHCpan) are license-restricted and must not be redistributed.

---

## Notes

* Before committing, review `.gitignore` and update ignore rules according to your project data/output policy.

---


