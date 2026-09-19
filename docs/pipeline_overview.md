# Pipeline overview

TETRIS-seq uses **two capture panels**, processed separately, and analysis runs
at **two levels**: per individual sample, and across groups of samples.

| Panel | Consensus used | Detects |
|---|---|---|
| **SNV panel** | SSCS **and** DCS | SNVs, indels, FLT3-ITDs |
| **mCA / translocation panel** (the "CNV panel") | SSCS | mCAs, KMT2A-PTDs, chromosomal rearrangements |

From raw paired-end FASTQs both panels first build single-strand (SSCS) and, for
the SNV panel, duplex (DCS) consensus reads. `<m>`/`<d>` below are the SSCS/DCS
minimum family sizes (default 3).

## Top-level entry points (`scripts/`)

| Script | Panel / purpose |
|---|---|
| `Watson_code_SNV_panel_v1.7.sh` | SNV panel: SSCS + DCS, SNV/indel calling and annotation |
| `Watson_code_CNV_panel_v2.3.sh` | mCA/translocation panel: SSCS, read depths, BAFs, KMT2A-PTD |
| `Watson_code_SNV_panel_for_FLT3_calling_v1.1.sh` | SSCS + DCS consensus calling for the FLT3-ITD arm (no SAM flag filtering) |
| `Watson_code_Pindel_FLT3_caller.sh` | FLT3-ITD detection (Pindel) on the consensus BAMs from the script above |
| `Watson_code_environment_setup_v2.6.sh` | sourced by the above; defines all tool/resource paths |

Example (SNV panel; `-h` for all options):
```bash
cd scripts
./Watson_code_SNV_panel_v1.7.sh -s ../examples/Sample_indexes_example.csv \
    -o <run_label> -f <read1.fastq> -g <read2.fastq> -t 8
```

---

## SNV panel: SSCS and DCS

### Per individual sample
- **SNV / indel calling (DCS and SSCS):** two callers are run,
  **VarDictJava** and a **custom variant caller**
  (`pipeline_tools/Watson_code_*VCF_SNP_calling*`), followed by ANNOVAR
  annotation (`Watson_code_DCS_annotation`, `Watson_code_SSCS_annotation`).
- **FLT3-ITD calling** with **Pindel** (`Watson_code_Pindel_FLT3_caller.sh`),
  on **both SSCS and DCS**. This arm uses its own consensus BAMs, built by
  `Watson_code_SNV_panel_for_FLT3_calling_v1.1.sh`: the SNV panel's consensus
  callers filter reads on their SAM flags, which discards the discordant and
  heavily soft-clipped reads that carry an ITD.

### Across groups of samples
- **In silico error model** fitted across **~40 samples** from the custom
  variant-caller results
  (`noise_correction_model/Duplex_Error_Model_initial_variant_calling_v5.py`).
- **Post-processing** of the error-corrected SNVs and the VarDictJava indels
  (`noise_correction_model/*post-model*` and trajectory notebooks).

## mCA / translocation panel: SSCS

### Per individual sample
- **KMT2A-PTD calling** on SSCS (`pipeline_tools/Watson_code_KMT2A-PTD_calling_v1.9.py`).
- **mCA calling:** compute **BAF and LRR** values (initially **not** PON-normalised)
  from the CNV panel (`Watson_code_sample_normalised_read_depths`,
  `Watson_code_CNV_SNP_plotting`); call mCAs **once the PON is available**
  (`mCA_caller/`), with **longitudinal phasing / mCA calling** if an mCA is
  present at the final timepoint. This phasing and earlier-timepoint calling is
  done entirely in `mCA_caller/` (the phased notebook,
  `mCA_caller/watson_code_mCA_caller_phased_v5.ipynb`).
- **Chromosomal rearrangement calling** on the **mapped BAM (pre-SSCS)**
  (`chromosomal_rearrangement_caller/`). If a rearrangement is called as real, a
  **dedicated SSCS BAM is built over the breakpoint regions only** and the caller
  is re-run on it to get a more accurate VAF estimate. That consensus comes from
  `chromosomal_rearrangement_caller/Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1_specific_regions.py`
  — the same no-SAM-flag-filter caller the FLT3 arm uses, with a `--regions`
  restriction. It is **not** the CNV panel's own SSCS BAM: that one is called by
  `Watson_code_SSCS_calling_2.1.py`, which keeps only properly paired reads and so
  discards the discordant pairs that span a breakpoint.

### Across groups of samples
- **Panel of normals (PON)** created for mCA calling
  (`mCA_caller/watson_code_create_PON.ipynb`): the PON is built from QC-passing
  final-timepoint controls, and each sample's LRR is computed against it with
  **leave-one-out** normalisation for samples that are themselves in the PON.

---

## Manual curation steps (not automated in this repository)

A few steps described in the manuscript were performed by **manual curation** and
are intentionally **not** implemented as code here. The scripts produce the
intermediate values; the final calls were made by hand:

- **FLT3-ITD calls.** Pindel is run independently on the SSCS and DCS BAMs
  (`Watson_code_Pindel_FLT3_caller.sh`), producing per-arm candidate calls. The
  manuscript's requirement that an ITD be detected in **both** the SSCS and DCS
  files to be called real was applied by **manual curation** of these outputs.
- **KMT2A-PTD cell fraction.** The caller
  (`pipeline_tools/Watson_code_KMT2A-PTD_calling_v1.9.py`) computes and reports
  the exon 3:27 read-depth ratio *R*. The cell fraction `2 × (R − 1)` was
  **calculated manually** for the sample(s) in which a PTD was present.
- **mCA germline exclusion.** CN-LOH events at ~100% cell fraction at both the
  earliest and latest available timepoints were classified as likely germline and
  excluded by **manual curation**. (The automated caller additionally uses a
  fixed known-germline-region skip list.)

## Consensus / alignment stages (both panels)

```
FASTQ (R1, R2)
   │  Picard FastqToSam            reads packed into an unaligned BAM, so that
   ▼                               per-read tags can be carried through
Unmapped BAM
   │  fgbio ExtractUmisFromBam     read structure 3M2S146T: 3 bp inline UMI ->
   │                               ZA/ZB, concatenated into RX; 2 bp spacer
   │                               dropped; 146 bp template kept
   │  Picard MarkIlluminaAdapters  3' adapter read-through marked in XT
   ▼
Unmapped BAM (tagged: RX, XT)
   │
   │  ── alignment round trip, run as one pipe ──
   │  Picard SamToFastq        stream out interleaved FASTQ, adapter bases
   │                           (CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X) masked
   │                           to X so BWA cannot align them
   │  | bwa mem -p             align to hg19
   │  | Picard MergeBamAlignment
   │                           merge the alignments BACK ONTO the unmapped BAM
   │                           -> the RX/ZA/ZB tags survive alignment
   ▼
mapped + merged BAM        (coordinate-sorted, indexed; also the BAM the
   │                        rearrangement caller reads, since it is the only
   │                        one still holding raw discordant / split reads)
   │  Picard CollectInsertSizeMetrics + panel read-coverage plots  [QC only]
   ▼
SSCS calling (Watson_code_SSCS_calling_2.1.py)  ──►  SSCS BAM
   │                                                   │
   │  the SSCS BAM as emitted is the input to both branches below
   │                                                   │
   ├───────────────────────────────────────────────────┤
   │                                                   │
   ▼                                                   ▼
DCS calling (Watson_code_DCS_calling_1.4.py)      re-align SSCS BAM
   │   duplex consensus (SNV panel only)           (unmap, sort, mark adapters,
   ▼                                                BWA, MergeBamAlignment)
re-align DCS BAM (same cycle)                          │
   │                                                   ▼
   ▼                                            ClipBam + GATK indel realignment
ClipBam + GATK indel realignment                       │
   │                                                   ▼
   ▼                                            SSCS variant calling
DCS variant calling
```

The `SamToFastq | bwa mem | MergeBamAlignment` round trip is used because an aligner
cannot read a BAM, and a plain FASTQ round trip would discard the UMI tags that were
just written. Merging the alignments back onto the unaligned BAM keeps every tag. The
same three-tool pipe is used again on each consensus BAM: a consensus read is a new
sequence, so it has to be re-aligned rather than inherit the coordinates of the reads
it was built from.

The re-alignment / ClipBam / indel-realignment cycle is applied **once per consensus
level**, and each realigned BAM feeds only that level's variant calling. DCS reads are
built from the SSCS BAM *before* it is re-aligned. The FLT3-ITD arm takes its consensus
BAMs straight from MergeBamAlignment, with no ClipBam and no indel realignment.

### One consensus caller, two arms

The FLT3-ITD arm and the chromosomal rearrangement arm share a single consensus
caller: `pipeline_tools/Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1.py`.
It is the standard SSCS caller with the SAM flag filter removed (and the adapter
clip position, `XT`, carried onto the consensus read); the grouping key, the
consensus threshold and the family-size rules are unchanged. Both arms need this
because the reads that carry an ITD, or that span a rearrangement breakpoint, are
frequently not flagged as properly paired, and `Watson_code_SSCS_calling_2.1.py`
keeps only flags 99/147/83/163.

The two arms differ only in scope:

| Arm | Script | Scope |
|---|---|---|
| FLT3-ITD | `pipeline_tools/Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1.py` | the whole panel |
| Rearrangements | `chromosomal_rearrangement_caller/Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1_specific_regions.py` | only the breakpoint regions passed to `--regions`, and only once a hit is found in the raw BAM |

The `_specific_regions` variant differs from the other by one change: instead of
sweeping the whole BAM, it iterates the requested chrom/start/end triples. The
consensus logic is identical.

## Code locations by variant class

| Variant class | Code |
|---|---|
| SNVs / indels | `pipeline_tools/Watson_code_*VCF_SNP_calling*`, VarDictJava; error model in `noise_correction_model/` |
| FLT3-ITD | `scripts/Watson_code_SNV_panel_for_FLT3_calling_v1.1.sh` then `scripts/Watson_code_Pindel_FLT3_caller.sh` + `pipeline_tools/Watson_code_FLT3_*` |
| KMT2A-PTD | `pipeline_tools/Watson_code_KMT2A-PTD_calling_v1.9.py` |
| mCAs (BAF/LRR + PON) | per-sample depths/BAFs `pipeline_tools/Watson_code_sample_normalised_read_depths*`, `*CNV_SNP_plotting*`; PON + LRR `mCA_caller/watson_code_create_PON.ipynb`; calling in `mCA_caller/` |
| Chromosomal rearrangements | `chromosomal_rearrangement_caller/` |

See [data_layout.md](data_layout.md) for input/output directory structure and the
reorganisation step between pipeline output and the analysis code.
