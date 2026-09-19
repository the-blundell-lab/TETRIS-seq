# Installation

The pipeline orchestrates several external bioinformatics tools (via Bash) plus
a set of custom Python scripts. Below are the exact tool versions used in the
manuscript and where the pipeline expects to find them.

## 1. Python environments

Two Python environments were used:

| Purpose | Python | Definition |
|---|---|---|
| Sequencing pipeline (`Watson_code_*.py` invoked by the `.sh` scripts) | 3.7.4 | [`environment_sequencing.yml`](../environment_sequencing.yml) - also provides BWA, samtools, tabix and Java 8 |
| Data analysis + figures (`noise_correction_model/`) | 3.11.13 | [`environment_analysis.yml`](../environment_analysis.yml) / [`requirements.txt`](../requirements.txt) |

```bash
# sequencing pipeline
conda env create -f environment_sequencing.yml
conda activate tetris-seq-pipeline

# analysis & figures (separate environment)
conda env create -f environment_analysis.yml
conda activate tetris-seq-analysis
```

Python package versions used for analysis/figures: Biopython 1.85, NumPy 2.0.2,
Matplotlib 3.10.3, pandas 2.3.0, SciPy 1.16.0, scikit-learn 1.7.1,
pyfaidx 0.8.1.4, pysam 0.23.3.

## 2. External command-line tools

Install the following at the exact versions used in the manuscript. Several are
JARs or large binaries and are **not** included in this repository.

| Tool | Version | Notes |
|---|---|---|
| BWA | 0.7.17 | read alignment |
| Samtools | 1.10 and 1.11 | both used |
| Picard | 2.18.15 | `picard.jar` (conda environment, or your own copy) |
| GATK | 3.8 | `GenomeAnalysisTK.jar` (3.8-1-0-gf15c1c3ef) - download separately |
| fgbio | 1.3.0 | `fgbio.jar` (conda environment, or your own copy) |
| VarDictJava | 1.8.2 | SNV/indel calling (installed by the conda environment, called as `vardict-java`) |
| ANNOVAR | June 2020 release | with `humandb/` databases |
| Pindel | v0.3 | FLT3-ITD detection |

## 3. Reference genome

`Homo_sapiens_assembly19.fasta` (Broad b37 / GRCh37), together with its `.fai`
and `.dict` indexes and the BWA index files. This is ~3 GB and is not included
here; download it from the Broad resource bundle (or the Zenodo archive linked
in the main README).

## 3b. dbSNP interval files

Both panel scripts annotate called positions against dbSNP, and expect a directory
of per-chromosome interval files (UCSC dbSNP153 "common" track, hg19), named:

```
UCSC_dbSNP153common_hg19_chr1.interval
UCSC_dbSNP153common_hg19_chr2.interval
...
```

Each file is the UCSC table dump for that chromosome, in 0-based coordinates
(the pipeline adds 1 when building its lookup). Put them in `$EXTERNAL_TOOLS/dbSNP`, or set the location in `config/config.sh`.

They can be regenerated from the UCSC Table Browser (assembly hg19, group
"Variation", track "Common SNPs(153)", output format "all fields from selected
table"), one query per chromosome.

## 4. Tell the pipeline where everything is

Copy the example config and edit the paths:

```bash
cp config/config.sh.example config/config.sh
$EDITOR config/config.sh
```

Two separate locations are used:

| variable | holds | default |
|---|---|---|
| `PIPELINE_TOOLS` | the repo's own `Watson_code_*.py` scripts and panel BED/CSV resources | [`../pipeline_tools/`](../pipeline_tools), found relative to `scripts/` |
| `EXTERNAL_TOOLS` | what you install: the JARs, reference genome, `annovar/`, `pindel/`, `dbSNP/` | `$HOME/Pipeline_tools` |

So in normal use you only set `EXTERNAL_TOOLS`; `PIPELINE_TOOLS` looks after itself, and the repo's
scripts stay current when you `git pull`.

Check what the pipeline can find before you run anything:

```bash
scripts/check_setup.sh
```

## 5. Script completeness

All custom scripts referenced by `Watson_code_environment_setup_v2.6.sh` are now
present in `pipeline_tools/` (27 Python scripts). One version note:

- The environment setup originally referenced
  `Watson_code_CNV_phased_SNP_plotting_v1.2.py`, but only **v1.1** exists; the
  reference now points to `Watson_code_CNV_phased_SNP_plotting_v1.1.py`.

The only things not stored in the repository are the external tools and large
reference files described in sections 2–3 (Picard, GATK, fgbio, Pindel,
ANNOVAR `humandb/`, the reference genome, and the optional `dbSNP/` directory).
Install these separately and point to them via `config/config.sh`.
