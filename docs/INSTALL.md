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
| Command-line tools for the pipeline | 3.7 | [`environment_sequencing.yml`](../environment_sequencing.yml) - provides BWA, samtools, tabix, Java 8, Picard, fgbio, VarDictJava and Pindel. **Its Python cannot run the pipeline scripts — see section 1b.** |
| Data analysis + figures (`noise_correction_model/`) | 3.11.13 | [`environment_analysis.yml`](../environment_analysis.yml) / [`requirements.txt`](../requirements.txt) |

[`environment_sequencing.lock.yml`](../environment_sequencing.lock.yml) records the
exact versions resolved on a working Linux install, for when a result needs
reproducing precisely rather than merely re-running.

```bash
# sequencing pipeline
mamba env create -f environment_sequencing.yml     # or: conda env create -f ...
conda activate tetris-seq-pipeline

# analysis & figures (separate environment)
mamba env create -f environment_analysis.yml
```

## 1b. Python for the pipeline scripts (required)

The consensus callers store the UMI-grouped reads with Python's `shelve`, which
needs a real dbm backend — `gdbm` or `ndbm`. **conda-forge's Python does not ship
the gdbm extension** (it is GPL, so conda-forge omits it, as Debian does), and
when no backend is available Python silently falls back to `dbm.dumb`. That
fallback rewrites its index on every sync, so consensus calling appears to run —
at a few percent CPU, with the temp files barely growing — and never finishes.

So the pipeline uses two environments together: **conda for the command-line
tools, and a virtualenv built from a system Python for the scripts themselves**.

```bash
scripts/setup_python_env.sh
```

That finds a Python with a working backend, creates `./tetris-py` (bootstrapping
pip if your distribution lacks `python3-venv`), installs `pysam`, `pyfaidx`,
`numpy`, `pandas`, `scipy`, `biopython`, `matplotlib`, `networkx`, `fuzzywuzzy`,
`python-Levenshtein`, `xlsxwriter` and `openpyxl`, and verifies the result. Then activate both, conda first so the virtualenv's `python` wins:

```bash
conda activate tetris-seq-pipeline
source tetris-py/bin/activate
```

Check the split is right — `python` from the virtualenv, everything else from
conda:

```bash
which python bwa samtools vardict-java pindel
```

`scripts/check_setup.sh` reports the backend in use and fails if it is
`dbm.dumb`. You can also check any interpreter directly:

```bash
python -c "import dbm.gnu"      # or dbm.ndbm; either is fine, dbm.dumb is not
```

If no suitable Python exists on the machine, install your distribution's
bindings (`sudo apt install python3-gdbm`, `sudo dnf install python3-gdbm`) and
re-run the setup script.

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
| ANNOVAR | June 2020 release | registration required; databases via `scripts/fetch_annovar_databases.sh` (section 3c) |
| Pindel | 0.2.5b9 | FLT3-ITD detection (installed by the conda environment; [source](https://github.com/genome/pindel)) |

## 3. Reference genome

`Homo_sapiens_assembly19.fasta` (Broad b37 / GRCh37) with its `.fai`, `.dict`
and BWA index files. One command:

```bash
scripts/fetch_reference.sh "$EXTERNAL_TOOLS"
```

By default this downloads a prebuilt copy (FASTA, `.fai`, `.dict` and the BWA
index, 4 GB compressed, ~8 GB unpacked) from the Zenodo archive that accompanies
this repository:

> **TETRIS-seq reference data** — <https://doi.org/10.5281/zenodo.22846473>

The script verifies the md5 before unpacking, so no indexing step is needed.

If you would rather build it yourself, `--from-broad` downloads the FASTA from
the Broad's public bucket and runs `bwa index` locally (~3 GB download, about an
hour of indexing; needs the conda environment active for `bwa`):

```bash
scripts/fetch_reference.sh "$EXTERNAL_TOOLS" --from-broad
```

If you already hold an indexed copy, put it in `$EXTERNAL_TOOLS` and the script
will confirm and exit.

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

One command:

```bash
pipeline_tools/make_dbSNP_intervals.sh "$EXTERNAL_TOOLS/dbSNP"
```

By default this downloads the prepared files (693 MB, md5-verified) from the
same Zenodo archive as the reference genome:

> **TETRIS-seq reference data** — <https://doi.org/10.5281/zenodo.22846473>

To regenerate them from UCSC instead, use `--rebuild` (needs `bigBedToBed`,
~1.4 GB download, ~10 GB scratch space):

```bash
conda install -c bioconda ucsc-bigbedtobed
pipeline_tools/make_dbSNP_intervals.sh "$EXTERNAL_TOOLS/dbSNP" --rebuild
```

## 3c. ANNOVAR databases

The SNV panel annotates against six hg19 databases:

```
refGene  cosmic92_coding  cosmic92_noncoding  exac03  gnomad_genome  clinvar_20200316
```

ANNOVAR itself needs a free academic
[registration](https://www.openbioinformatics.org/annovar/annovar_download_form.php);
unpack it to `$EXTERNAL_TOOLS/annovar`, then:

```bash
scripts/fetch_annovar_databases.sh
```

That downloads refGene, exac03, gnomad_genome and clinvar_20200316 into
`humandb/`. The two COSMIC databases are licensed and cannot be redistributed,
so they are built from your own COSMIC v92 (GRCh37) download - register at
[cancer.sanger.ac.uk](https://cancer.sanger.ac.uk/cosmic/register), gunzip the
four files into one directory and pass it in:

```bash
scripts/fetch_annovar_databases.sh --cosmic-dir /path/to/cosmic_v92
```

| COSMIC file | builds |
|---|---|
| `CosmicMutantExport.tsv` + `CosmicCodingMuts.vcf` | `hg19_cosmic92_coding.txt` |
| `CosmicNCV.tsv` + `CosmicNonCodingVariants.vcf` | `hg19_cosmic92_noncoding.txt` |

A newer COSMIC release works just as well, but the COSMIC columns in the
published variant call tables came from v92; if you use another release, change
the `-protocol` names in `scripts/Watson_code_SNV_panel_v1.7.sh` to match.

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

All custom scripts referenced by `Watson_code_environment_setup_v2.6.sh` are
present in `pipeline_tools/` (27 Python scripts).

The only things not stored in the repository are the external tools and large
reference files described in sections 2–3 (GATK, ANNOVAR `humandb/`, the
reference genome and the `dbSNP/` directory).
Install these separately and point to them via `config/config.sh`.
