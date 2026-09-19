# Installation

The pipeline orchestrates several external bioinformatics tools (via Bash) plus
a set of custom Python scripts. Below are the exact tool versions used in the
manuscript and where the pipeline expects to find them.

## 1. Python environments

**conda is required.** We recommend [Miniforge](https://github.com/conda-forge/miniforge),
which uses conda-forge by default: Anaconda's own `defaults` channels now sit behind a
terms-of-service prompt and a licence that many institutions cannot accept.

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
source ~/miniforge3/etc/profile.d/conda.sh
echo 'source ~/miniforge3/etc/profile.d/conda.sh' >> ~/.bashrc
```

If you already have Miniconda or Anaconda and hit
`CondaToSNonInteractiveError`, point conda at the community channels instead:

```bash
printf 'channels:\n  - conda-forge\n  - bioconda\nchannel_priority: strict\n' > ~/.condarc
```

| Purpose | Python | Definition |
|---|---|---|
| Sequencing pipeline (`Watson_code_*.py` invoked by the `.sh` scripts) | 3.7 | [`environment_sequencing.yml`](../environment_sequencing.yml) - also provides BWA, samtools, tabix, Java 8, Picard, fgbio and VarDictJava |
| Data analysis + figures (`noise_correction_model/`) | 3.11.13 | [`environment_analysis.yml`](../environment_analysis.yml) / [`requirements.txt`](../requirements.txt) |

```bash
# sequencing pipeline
mamba env create -f environment_sequencing.yml     # or: conda env create -f ...
conda activate tetris-seq-pipeline

# analysis & figures (separate environment)
mamba env create -f environment_analysis.yml
```

## 2. External command-line tools

Install the following at the exact versions used in the manuscript. Several are
JARs or large binaries and are **not** included in this repository.

| Tool | Version | Notes |
|---|---|---|
| BWA | 0.7.17 | read alignment |
| Samtools | 1.10 and 1.11 | both used |
| Picard | 2.18.15 | `picard.jar` (conda environment, or your own copy) |
| GATK | 3.8 | `GenomeAnalysisTK.jar`, build 3.8-1-0-gf15c1c3ef - [download](https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2) |
| fgbio | 1.3.0 | `fgbio.jar` (conda environment, or your own copy) |
| VarDictJava | 1.8.2 | SNV/indel calling (installed by the conda environment, called as `vardict-java`) |
| ANNOVAR | June 2020 release | with `humandb/` databases - [registration form](https://www.openbioinformatics.org/annovar/annovar_download_form.php) |
| Pindel | v0.3 | FLT3-ITD detection - [source](https://github.com/genome/pindel) |

## 3. Reference genome

`Homo_sapiens_assembly19.fasta` (Broad b37 / GRCh37) with its `.fai`, `.dict`
and BWA index files (~3 GB download, ~8 GB once indexed). One command:

```bash
conda activate tetris-seq-pipeline        # needs bwa for the indexing step
scripts/fetch_reference.sh "$EXTERNAL_TOOLS"
```

The script downloads from the Broad's public bucket, checks the FASTA is
complete, and builds the BWA index only if one is not already present (that step
takes about an hour). If you already hold an indexed copy, put it in
`$EXTERNAL_TOOLS` and the script will confirm and exit.

## 3b. dbSNP interval files

Both panel scripts annotate called positions against dbSNP, and expect a
directory of per-chromosome files (dbSNP build 153 "common", hg19), named:

```
UCSC_dbSNP153common_hg19_chr1.interval ... _chr22.interval, _chrX.interval
```

Columns are UCSC's `dbSnp153Common` fields (chrom, chromStart, chromEnd, rsID,
ref, altCount, alts, shiftBases, freqSourceCount, minorAlleleFreq, ...), with a
header line and 0-based coordinates (the pipeline adds 1 when building its
lookup).

Two ways to obtain them:

**Build them yourself** (needs `bigBedToBed`, ~1.4 GB download, ~10 GB scratch):

```bash
conda install -c bioconda ucsc-bigbedtobed
pipeline_tools/make_dbSNP_intervals.sh "$EXTERNAL_TOOLS/dbSNP"
```

**Or download the prepared copy** from the Zenodo archive linked in the main
README, and unpack it into `$EXTERNAL_TOOLS/dbSNP`.

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
