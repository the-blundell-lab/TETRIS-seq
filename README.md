# TETRIS-seq

**Duplex error-corrected sequencing pipeline for SNV, indel, copy-number,
FLT3-ITD and translocation detection.**

This repository contains the bioinformatic pipeline for TETRIS-seq, the
sequencing strategy used in Watson *et al.*, *"Evolutionary dynamics in the
decades preceding acute myeloid leukaemia"* (Nature, 2026). TETRIS-seq applies
duplex error-corrected sequencing to a panel covering clonal-haematopoiesis and
AML-associated alterations (gene mutations, chromosomal rearrangements and
mosaic chromosomal alterations, or mCAs) together with a custom in silico
de-noising algorithm, reaching over 1500× duplex consensus depth.

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
├── noise_correction_model/   Beta-binomial model + analysis/figure notebooks
│   └── data_files/           Small reference tables (COSMIC sites, metadata)
├── mCA_caller/               Mosaic chromosomal alteration (mCA) caller + simulations
├── chromosomal_rearrangement_caller/   Custom translocation / rearrangement caller
├── config/
│   └── config.sh.example     Copy to config.sh and set your tool paths
├── examples/                 Example sample-index sheet
├── docs/
│   ├── INSTALL.md              Tool versions + installation
│   ├── pipeline_overview.md    Stages, data flow, code by variant class
│   ├── pipeline_schematic.md   Colour-coded flow diagram (automated vs manual)
│   ├── running_the_pipeline.md Step-by-step run order (Stage A automated, Stage B manual)
│   ├── data_layout.md          Input/output directory structure
│   └── adapting_to_a_new_panel.md  Study-specific assumptions to change for reuse
├── environment_sequencing.yml   Conda env, Python 3.7.4 (pipeline)
├── environment_analysis.yml     Conda env, Python 3.11.13 (analysis/figures)
├── requirements.txt             pip equivalent of the analysis env
├── CITATION.cff
└── LICENSE                      BSD 3-Clause
```

## Quick start

```bash
# 1. Create the environments (see docs/INSTALL.md for the external tools)
conda env create -f environment_sequencing.yml      # Python 3.7.4
conda env create -f environment_analysis.yml        # Python 3.11.13

# 2. Point the pipeline at your tool installs and reference genome
cp config/config.sh.example config/config.sh
$EDITOR config/config.sh

# 3. Run the automated per-sample stages for one sample
#    (SNV and CNV are SEPARATE captures with SEPARATE FASTQ files)
conda activate tetris-seq-pipeline
cd scripts
./run_sample.sh -o <library_name> -t 8 \
    --snv-csv ../examples/Sample_indexes_example.csv \
    --snv-r1 <snv_read1.fastq> --snv-r2 <snv_read2.fastq> \
    --cnv-csv ../examples/Sample_indexes_example.csv \
    --cnv-r1 <cnv_read1.fastq> --cnv-r2 <cnv_read2.fastq>
```

Process one sample at a time, since consensus calling is memory-hungry, and loop
the wrapper over samples to do a whole lane. Once the automated stages are done,
the cross-sample and manual steps (noise model, PON and mCA calling,
rearrangement calling, curation) are triggered by hand. The run order is set out
below.

See **[docs/running_the_pipeline.md](docs/running_the_pipeline.md)** for the
step-by-step run order (automated Stage A vs manual Stage B),
**[docs/pipeline_schematic.md](docs/pipeline_schematic.md)** for the colour-coded
flow diagram, **[docs/pipeline_overview.md](docs/pipeline_overview.md)** for the
full data flow and code by variant class,
**[docs/data_layout.md](docs/data_layout.md)** for the input/output directory
structure, **[docs/INSTALL.md](docs/INSTALL.md)** for tool versions and
installation, and **[docs/adapting_to_a_new_panel.md](docs/adapting_to_a_new_panel.md)**
for the study-specific assumptions (genome build, naming, hardcoded coordinates)
and what would need changing to reuse this with a different panel.

## Software versions used

**External tools:** Picard 2.18.15, BWA 0.7.17, Samtools 1.10 and 1.11,
GATK 3.8, fgbio 1.3.0, Pindel v0.3, VarDictJava, snpEff/SnpSift,
ANNOVAR (June 2020). The reference genome is `Homo_sapiens_assembly19.fasta`
(Broad b37/GRCh37).

**Python:** 3.7.4 for the pipeline and 3.11.13 for analysis. Key packages are
Biopython 1.85, NumPy 2.0.2, pandas 2.3.0, SciPy 1.16.0, Matplotlib 3.10.3,
scikit-learn 1.7.1, pyfaidx 0.8.1.4 and pysam 0.23.3.

## Large files and data availability

The following are not stored in this repository and should be downloaded
separately (see [docs/INSTALL.md](docs/INSTALL.md)):

- Third-party JARs / binaries (Picard, GATK, fgbio, Pindel, ANNOVAR).
- The reference genome (`Homo_sapiens_assembly19.fasta`, ~3 GB).
- The full per-sample model output (~9 GB).

A small example dataset is bundled in this repository for testing. As stated in
the paper's data-availability statement, somatic SNV, indel and FLT3-ITD calls
are provided in Supplementary Table 6 and somatic mCAs in Supplementary Table 7,
with the corresponding somatic VCFs deposited at Zenodo
(DOI: 10.5281/zenodo.22262496). Germline variant calls, the raw sequencing reads
and the exact participant timing data are not open: they are available through
the study's controlled-access route (EGA), as set out in the paper.

## License

BSD 3-Clause. See [LICENSE](LICENSE).
