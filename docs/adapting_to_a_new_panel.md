# Assumptions & adapting to a new panel

This pipeline is the exact code used in the manuscript for the TETRIS-seq
TWIST panels on the UKCTOCS samples. It is not a turnkey, panel-agnostic
tool. Several things are specific to this study's panel, genome build, library
prep and sequencing-core naming. This page lists those assumptions so you can
judge what is reusable and what would need editing for a different panel.

## Hard assumptions (would require code edits to change)

### 1. Genome build: hg19 / b37 (GRCh37)
The build is assumed throughout. It is not a configurable parameter beyond
the reference FASTA path:
- Reference `Homo_sapiens_assembly19.fasta` (Broad b37) in `config/config.sh`.
- Pindel is called with `-R hg19` (`scripts/Watson_code_Pindel_FLT3_caller.sh`).
- dbSNP filtering reads `…/UCSC_dbSNP153common_hg19_<chr>.interval`
  (`pipeline_tools/Watson_code_*VCF_SNP_calling*.py`).
- All hardcoded coordinates (below) are hg19.

Running on a different build (e.g. hg38) requires editing these and replacing
all hardcoded coordinates.

### 2. FASTQ filename convention: `<library>.<UDI>.…fastq`
Sample and library names are parsed from the FASTQ filename by splitting on `.`:
- `pipeline_tools/Watson_code_sample_name_from_UDI_index_v1.py`:
  UDI = `fastqname.split('.')[1]`, looked up in the sample-index CSV.
- `pipeline_tools/Watson_code_library_name_v1.py`:
  library = `fastqname.split('.')[0]` (e.g. `SLX-20124`).

Differently named FASTQs (e.g. the Illumina default
`Sample_S1_L001_R1_001.fastq.gz`) will produce wrong sample/library names.

### 3. Duplex-UMI read-name structure
Consensus calling assumes the two UMIs are encoded in the read name as
`UMI1-UMI2:…`:
- `UMI = qname.split(':')[0]`; `UMI1 = UMI.split('-')[0]`, `UMI2 = UMI.split('-')[1]`
  (`pipeline_tools/Watson_code_SSCS_calling_2.1.py`, `Watson_code_DCS_calling_1.4.py`).

A different UMI scheme / read-name layout breaks SSCS/DCS calling.

### 4. Hardcoded genomic coordinate tables
Some scripts contain panel/gene-specific coordinates baked into the code (hg19).
These are not derived from the panel BED and must be edited for new genes:
- `pipeline_tools/Watson_code_translocation_breakpoint_coverage_v1.2.py`
- `pipeline_tools/Watson_code_SNV_panel_read_coverage.py` (per-gene exon tables, e.g. RAD21)
- `pipeline_tools/Watson_code_KMT2A-PTD_calling_v1.9.py`
- `chromosomal_rearrangement_caller/` (targeted rearrangement regions)

### 5. Gene-specific callers
`Watson_code_Pindel_FLT3_caller.sh` (FLT3-ITD) and the KMT2A-PTD caller are
specific to those genes by design; they are only meaningful if the new panel
covers FLT3 / KMT2A.

## Panel-specific files (swappable, but study-specific)

These are passed in via `config/config.sh` / arguments, so they can in principle
be replaced for a new panel without code changes. Note that the supplied versions
are specific to the TWIST v2 panels:

| File (in `pipeline_tools/` unless noted) | Purpose |
|---|---|
| `TWIST_SNV_panel-TE-92048328_hg19.bed` (+ `_GATK_compatible.bed`) | SNV panel target regions |
| `TWIST_CNV_panel_TE-95031423_h19.bed` (+ `_GATK_compatible.bed`) | mCA/translocation panel target regions |
| `TWIST_CNV_panel_targeted_SNPs_positions.csv` | SNP positions for BAF / mCA calling |
| `TWIST_v2_transcript_details.csv` | transcripts used for annotation |
| `KMT2A_coordinates_hg19.csv` | KMT2A exon coordinates for PTD calling |
| `chromosome_ideogram_hg19.txt` | ideogram for plots |
| `chromosomal_rearrangement_caller/Translocation_regions_of_interest.bed` | rearrangement target regions |
| dbSNP `UCSC_dbSNP153common_hg19_*.interval` (provide separately) | common-SNP filtering |

## What IS already configurable

Tool paths, reference path, the panel BED files above, the sample-index CSV, and
the SSCS/DCS family-size / quality thresholds are all set via `config/config.sh`
and the script options (`-h`). See [INSTALL.md](INSTALL.md) and
[pipeline_overview.md](pipeline_overview.md).

## Checklist for a new panel (same hg19 build, same UMI design, same FASTQ naming)

1. Replace the panel BED files and point `config.sh` at them.
2. Replace the targeted-SNP, transcript-details and ideogram files.
3. Update the hardcoded coordinate tables (§4) for any new genes.
4. Confirm FLT3/KMT2A coverage if you want those callers.
5. Provide an hg19 reference + dbSNP intervals.

A different genome build, UMI scheme or FASTQ-naming convention additionally
requires the code edits described in §1–3.
