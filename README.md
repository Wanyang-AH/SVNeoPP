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
* MHCflurry input preparation → prediction → parse → filter
* DeepImmuno input preparation → prediction → parse → filter

Planned / TODO (not yet fully integrated):

* `FragPipe` (mass-spectrometry database search)

> Note: This repository expects some **large resources / licensed tools** (e.g. netMHCpan) to be provided locally and configured via `config/config.yaml`.

---

## Quickstart (recommended: run the included L041 module test)

This is the fastest way to verify your environment and the netMHCpan / MHCflurry / DeepImmuno modules, without setting up full raw data paths.

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

### 3) Configure netMHCpan / MHCflurry / DeepImmuno resources in `config/config.yaml`

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
* `references.mhcflurry_model_tarball`: `resources/mhcflurry/models_class1_presentation.20200611.tar.bz2`
* `params.mhcflurry.docker_image`: Docker image ID containing MHCflurry CLI
* `references.deepimmuno_home`: `resources/DeepImmuno`
* `params.deepimmuno.python`: Python executable used to run DeepImmuno scripts (default: `python`)

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

Option C (strict module-only dry-run; MHCflurry on L041):

```bash
snakemake -n mhcflurry_filter_l041 --configfile config/config.yaml --rerun-triggers mtime
```

Option D (strict module-only dry-run; DeepImmuno on L041):

```bash
snakemake -n deepimmuno_filter_l041 --configfile config/config.yaml --rerun-triggers mtime
```

### 5) Real run (single-sample module test)

```bash
snakemake -j 1 netmhcpan_filter_l041 --configfile config/config.yaml \
  --allowed-rules \
    netchop_prepare_input netchop_run netchop_parse netchop_filter \
    netmhcpan_prepare_input netmhcpan_run netmhcpan_parse netmhcpan_filter \
    netmhcpan_filter_l041
```

MHCflurry single-sample module run:

```bash
snakemake -j 1 mhcflurry_filter_l041 --configfile config/config.yaml --rerun-triggers mtime
```

DeepImmuno single-sample module run:

```bash
snakemake -j 1 deepimmuno_filter_l041 --configfile config/config.yaml --rerun-triggers mtime
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
* MHCflurry settings are under `params.mhcflurry`.
* DeepImmuno settings are under `params.deepimmuno`.

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

## Expected outputs (module tests, L041)

The module tests produce filtered result tables under `results/<module>/l041/`, for example:

* `results/netmhcpan/l041/l041.filtered.csv`
* `results/mhcflurry/l041/l041.filtered.csv`
* `results/deepimmuno/l041/l041.filtered.csv`

Step-level logs are written to `logs/<module>/l041.*.log`.

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
* MHCflurry-related scripts:

  * `scripts/mhcflurry/prepare_input.py`
  * `scripts/mhcflurry/parse.py`
  * `scripts/mhcflurry/filter.py`
* DeepImmuno-related scripts:

  * `scripts/deepimmuno/prepare_input.py`
  * `scripts/deepimmuno/run.py`
  * `scripts/deepimmuno/parse.py`
  * `scripts/deepimmuno/filter.py`
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

### MHCflurry

Required configuration (example):

* `references.mhcflurry_model_tarball`: `resources/mhcflurry/models_class1_presentation.20200611.tar.bz2`
* `params.mhcflurry.docker_image`: Docker image ID containing MHCflurry CLI
* GitHub: `https://github.com/openvax/mhcflurry`

### DeepImmuno

Required configuration (example):

* `references.deepimmuno_home`: `resources/DeepImmuno`
* `params.deepimmuno.python`: `python`
* GitHub: `https://github.com/frankligy/DeepImmuno`

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

## Roadmap

Planned modules (not yet integrated):

* `FragPipe` for mass-spectrometry database search integration

---


## License

This project is licensed under the terms in the repository root `LICENSE` file.
Also note that some third-party tools/resources (e.g., netMHCpan) are license-restricted and must not be redistributed.

---

## Notes

* Before committing, review `.gitignore` and update ignore rules according to your project data/output policy.

---
