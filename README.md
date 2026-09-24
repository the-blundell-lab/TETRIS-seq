# TETRIS-seq

**Duplex error-corrected sequencing pipeline for SNV, indel, copy-number,
FLT3-ITD and translocation detection.**

This repository contains the bioinformatic pipeline for TETRIS-seq, the
sequencing strategy used in Watson *et al.*, *"Evolutionary dynamics in the
decades preceding acute myeloid leukaemia"*. TETRIS-seq applies duplex
error-corrected sequencing to a panel covering clonal-haematopoiesis and
AML-associated alterations (gene mutations, chromosomal rearrangements and
mosaic chromosomal alterations) together with a custom in silico de-noising
algorithm.

The pipeline builds single-strand (SSCS) and duplex (DCS) consensus reads, calls
and annotates variants, and applies a beta-binomial noise-correction model to
separate true low-frequency variants from sequencing error.

**Citation.** If you use this code, please cite the paper and the software
archive; the details are in [`CITATION.cff`](CITATION.cff).

![TETRIS-seq pipeline. Two capture panels are processed from separate FASTQ libraries. Per-sample stages are automated; cross-sample and longitudinal stages are triggered manually.](docs/figures/pipeline.png)

---

## Repository layout

```
TETRIS-seq/
├── scripts/                  Top-level pipeline entry points (Bash)
│   ├── run_sample.sh                           convenience driver: runs all automated per-sample stages for one sample
│   ├── Watson_code_SNV_panel_v1.7.sh
│   ├── Watson_code_CNV_panel_v2.3.sh
│   ├── Watson_code_SNV_panel_for_FLT3_calling_v1.1.sh   consensus reads for the FLT3-ITD arm
│   ├── Watson_code_Pindel_FLT3_caller.sh
│   └── Watson_code_environment_setup_v2.6.sh   (paths; sourced by the others)
├── pipeline_tools/           Custom Python scripts + panel BED/CSV resources
├── noise_correction_model/   Position-specific beta-binomial error model
│   └── data_files/           Small reference tables (COSMIC sites, metadata)
├── mCA_caller/               Mosaic chromosomal alteration (mCA) caller: PON, unphased and phased
├── chromosomal_rearrangement_caller/   Custom translocation / rearrangement caller
├── config/
│   └── config.sh.example     Copy to config.sh and set your tool paths
├── examples/                 Example sample-index sheet
├── docs/
│   ├── INSTALL.md              Tool versions + installation
│   ├── pipeline_overview.md    Stages, data flow, code by variant class
│   ├── running_the_pipeline.md Step-by-step run order (Stage A automated, Stage B manual)
│   ├── data_layout.md          Input/output directory structure
│   └── adapting_to_a_new_panel.md  Study-specific assumptions to change for reuse
├── environment_sequencing.yml   Conda env, Python 3.7.4 (pipeline)
├── environment_analysis.yml     Conda env, Python 3.11.13 (analysis/figures)
├── requirements.txt             pip equivalent of the analysis env
├── CITATION.cff
└── LICENSE                      BSD 3-Clause
```

## Requirements

Linux or macOS, and **conda** (Miniconda, Miniforge or Anaconda) to build the
environments. If you do not have it:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init bash && exec bash
```

The conda environment provides BWA, samtools, tabix, Java 8, Picard, fgbio,
VarDictJava and Pindel. The pipeline's own Python scripts run in a separate
virtualenv built by `scripts/setup_python_env.sh`, because conda-forge's Python
lacks the gdbm extension that consensus calling depends on — see
[section 1b of docs/INSTALL.md](docs/INSTALL.md). Four things are installed separately because of their size or
licensing: **GATK 3.8**, **ANNOVAR** (plus its `humandb/` databases), the
**b37 reference genome** and the **dbSNP interval files** - see
[docs/INSTALL.md](docs/INSTALL.md).

The reference genome (with a prebuilt BWA index) and the dbSNP interval files
are archived on Zenodo and downloaded by `scripts/fetch_reference.sh` and
`pipeline_tools/make_dbSNP_intervals.sh`:
**TETRIS-seq reference data**, <https://doi.org/10.5281/zenodo.22846473>.

## Quick start

```bash
# 1. Create the environments (see docs/INSTALL.md for the external tools)
conda env create -f environment_sequencing.yml      # Python 3.7.4
conda env create -f environment_analysis.yml        # Python 3.11.13

# 2. Point the pipeline at your tool installs and reference genome
cp config/config.sh.example config/config.sh
$EDITOR config/config.sh

# 3. Activate both environments, then check what the pipeline can find
source /path/to/TETRIS-seq/scripts/activate.sh   # both environments, in order
scripts/check_setup.sh

# 4. Run the automated per-sample stages.
#
#    The two panels are separate captures, sequenced as separate libraries, so
#    each is its own run over its own pair of FASTQ files. They share no input.

#    Output goes to <sample>_<UDI>/ relative to where you run these, so work
#    from a data directory outside the repository.
cd /path/to/your/working/directory

#    4a. SNV / indel panel (and FLT3-ITD, which re-uses the SNV panel's BAM)
/path/to/TETRIS-seq/scripts/run_sample.sh -o <snv_library_name> -t 8 --no-cnv \
    --snv-csv ../examples/Sample_indexes_example.csv \
    --snv-r1 <snv_read1.fastq> --snv-r2 <snv_read2.fastq>

#    4b. mCA / rearrangement panel - a different capture, a different library
/path/to/TETRIS-seq/scripts/run_sample.sh -o <cnv_library_name> -t 8 --no-snv --no-flt3 \
    --cnv-csv ../examples/Sample_indexes_example.csv \
    --cnv-r1 <cnv_read1.fastq> --cnv-r2 <cnv_read2.fastq>

#    (If you have both FASTQ pairs to hand, one call with all six arguments and
#     no --no-* flags does the same work in one go.)
```

The wrapper handles one sample per call, so loop it over samples to do a whole
lane. Consensus calling is memory-hungry, so how many you can run at once
depends on the machine. Once the automated stages are done,
the cross-sample and manual steps (noise model, PON and mCA calling,
rearrangement calling, curation) are triggered by hand. The run order is set out
below.

See **[docs/running_the_pipeline.md](docs/running_the_pipeline.md)** for the
step-by-step run order (automated Stage A vs manual Stage B),
**[docs/pipeline_overview.md](docs/pipeline_overview.md)** for the
full data flow and code by variant class,
**[docs/data_layout.md](docs/data_layout.md)** for the input/output directory
structure, **[docs/INSTALL.md](docs/INSTALL.md)** for tool versions and
installation, and **[docs/adapting_to_a_new_panel.md](docs/adapting_to_a_new_panel.md)**
for the study-specific assumptions (genome build, naming, hardcoded coordinates)
and what would need changing to reuse this with a different panel.

## Software versions used

**External tools:** Picard 2.18.15, BWA 0.7.17, Samtools 1.10 and 1.11,
GATK 3.8,
ANNOVAR (June 2020). The reference genome is `Homo_sapiens_assembly19.fasta`
(Broad b37/GRCh37).

**Python:** 3.7.4 for the pipeline and 3.11.13 for analysis. Key packages are
Biopython 1.85, NumPy 2.0.2, pandas 2.3.0, SciPy 1.16.0, Matplotlib 3.10.3,
scikit-learn 1.7.1, pyfaidx 0.8.1.4 and pysam 0.23.3.

## Large files

The following are not stored in this repository and are installed separately
(see [docs/INSTALL.md](docs/INSTALL.md)):

- Third-party JARs / binaries (Picard, GATK, fgbio, Pindel, ANNOVAR).
- The reference genome (`Homo_sapiens_assembly19.fasta`; 4 GB download,
  ~8 GB unpacked).
- The dbSNP files used for annotation (693 MB).

Running the pipeline also writes several GB per sample, into `TEMP/` and the
sample's output folder; see
[docs/running_the_pipeline.md](docs/running_the_pipeline.md).

## Manuscript analysis

The downstream analysis for Watson et al. lives in a separate repository,
[preAML_evolutionary_dynamics](https://github.com/the-blundell-lab/preAML_evolutionary_dynamics):
post-processing of these calls into variant trajectories, the fitness and
acquisition-age inference, phylogenies, Muller plots and the paper's figures.
The data-availability statement for the study's variant calls is there too.

The sequencing data this pipeline was written for, and the per-sample call files
it produces, are individual-level participant data: they are held under
controlled access in the European Genome-phenome Archive (EGA) and released to
approved researchers via a Data Access Committee. The somatic calls themselves
are open (supplementary tables and Zenodo); see the analysis repository for the
full statement.

## License

BSD 3-Clause. See [LICENSE](LICENSE).
