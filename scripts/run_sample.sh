#!/bin/bash
###############################################################################
# run_sample.sh
#
# Reference driver for ONE sample of TETRIS-seq.
#
# This is a *documented convenience wrapper*. It runs the fully-automated,
# per-sample stages of the pipeline for a single sample, in the canonical
# order, then prints the ordered list of MANUAL / CROSS-SAMPLE steps that are
# triggered by hand once you have processed as many samples as you want.
#
# IMPORTANT — two separate captures:
#   The SNV panel and the CNV/mCA panel are SEPARATE libraries with SEPARATE
#   FASTQ files. They are supplied independently below (--snv-* and --cnv-*)
#   and are never combined. The FLT3-ITD stage starts from the SNV panel's
#   mapped merged BAM, so it uses the SNV FASTQ / CSV.
#
# IMPORTANT — one sample at a time:
#   The consensus-calling steps are memory-hungry; running a whole lane at once
#   may exceed available RAM. Process a lane by calling this once per sample
#   (a shell loop or a cluster job array) — not all in parallel unless you know
#   you have the memory.
#
# It does NOT replace the individual scripts — it just chains them in order.
# See docs/running_the_pipeline.md and docs/pipeline_overview.md.
###############################################################################

set -euo pipefail

P="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_TOOLS="${PIPELINE_TOOLS:-$P/../pipeline_tools}"

SNV_PANEL="$P/Watson_code_SNV_panel_v1.7.sh"
CNV_PANEL="$P/Watson_code_CNV_panel_v2.3.sh"
FLT3_PREP="$P/Watson_code_SNV_panel_for_FLT3_calling_v1.1.sh"
FLT3_CALLER="$P/Watson_code_Pindel_FLT3_caller.sh"
SAMPLE_SCRIPT="$PIPELINE_TOOLS/Watson_code_sample_name_from_UDI_index_v1.py"
SAMPLE_UDI_SCRIPT="$PIPELINE_TOOLS/Watson_code_sample_name_from_UDI_index_with_UDI_v1.py"

# ---- defaults --------------------------------------------------------------
SNV_CSV=""
SNV_R1=""
SNV_R2=""
CNV_CSV=""        # defaults to SNV_CSV if left blank
CNV_R1=""
CNV_R2=""
LIBRARY=""
THREADS=8
RUN_SNV=1
RUN_CNV=1
RUN_FLT3=1

usage() {
  cat << EOF
Usage: $0 -o <library> [SNV inputs] [CNV inputs] [options]

Runs the automated stages for ONE sample, in order:
  1. SNV panel   (SSCS + DCS; SNV/indel calling + annotation)   -- SNV FASTQ
  2. FLT3-ITD    (own SSCS + DCS consensus, then Pindel)         -- SNV FASTQ
  3. CNV panel   (SSCS; read depths, BAFs, KMT2A-PTD, coverage) -- CNV FASTQ

The SNV and CNV panels are SEPARATE captures with SEPARATE FASTQ files.

REQUIRED:
   -o, --library <name>   library/run label (e.g. SLX_20125)

SNV panel inputs (needed unless --no-snv AND --no-flt3):
       --snv-csv <csv>    sample_indexes.csv for the SNV library
       --snv-r1  <fq>     SNV panel read-1 FASTQ
       --snv-r2  <fq>     SNV panel read-2 FASTQ

CNV panel inputs (needed unless --no-cnv):
       --cnv-csv <csv>    sample_indexes.csv for the CNV library
                          (default: same as --snv-csv)
       --cnv-r1  <fq>     CNV panel read-1 FASTQ
       --cnv-r2  <fq>     CNV panel read-2 FASTQ

OPTIONS:
   -t, --threads <cpus>   threads (default: $THREADS)
       --no-snv           skip the SNV panel
       --no-cnv           skip the CNV panel
       --no-flt3          skip the FLT3-ITD caller
   -h, --help             show this help

To process a whole lane, call this once per sample, e.g.:

   for r1 in /snv_fastqs/*_R1_*.fastq.gz; do
       s=\$(basename "\$r1" | sed 's/_R1_.*//')          # sample stem
       $0 -o SLX_20125 -t 8 \\
          --snv-csv snv_samples.csv --snv-r1 "\$r1" --snv-r2 "\${r1/_R1_/_R2_}" \\
          --cnv-csv cnv_samples.csv \\
          --cnv-r1 /cnv_fastqs/\${s}_R1_.fastq.gz \\
          --cnv-r2 /cnv_fastqs/\${s}_R2_.fastq.gz
   done
EOF
}

# ---- arg parsing (long flags) ----------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--library) LIBRARY="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    --snv-csv)    SNV_CSV="$2"; shift 2;;
    --snv-r1)     SNV_R1="$2";  shift 2;;
    --snv-r2)     SNV_R2="$2";  shift 2;;
    --cnv-csv)    CNV_CSV="$2"; shift 2;;
    --cnv-r1)     CNV_R1="$2";  shift 2;;
    --cnv-r2)     CNV_R2="$2";  shift 2;;
    --no-snv)     RUN_SNV=0;  shift;;
    --no-cnv)     RUN_CNV=0;  shift;;
    --no-flt3)    RUN_FLT3=0; shift;;
    -h|--help)    usage; exit 0;;
    *) echo "Error: unknown argument '$1'"; echo; usage; exit 1;;
  esac
done

# CNV CSV defaults to the SNV CSV when the same sample sheet is used for both
[[ -z "$CNV_CSV" ]] && CNV_CSV="$SNV_CSV"

# ---- validation ------------------------------------------------------------
[[ -z "$LIBRARY" ]] && { usage; echo; echo "Error: -o/--library is required"; exit 1; }

need_snv_inputs=0
[[ $RUN_SNV -eq 1 || $RUN_FLT3 -eq 1 ]] && need_snv_inputs=1

if [[ $need_snv_inputs -eq 1 ]]; then
  [[ -z "$SNV_CSV" ]] && { echo "Error: --snv-csv is required for the SNV/FLT3 stages"; exit 1; }
  [[ -z "$SNV_R1"  ]] && { echo "Error: --snv-r1 is required for the SNV/FLT3 stages"; exit 1; }
  [[ -z "$SNV_R2"  ]] && { echo "Error: --snv-r2 is required for the SNV/FLT3 stages"; exit 1; }
  [[ -f "$SNV_CSV" ]] || { echo "Error: SNV CSV not found: $SNV_CSV"; exit 1; }
  [[ -f "$SNV_R1"  ]] || { echo "Error: SNV read1 not found: $SNV_R1"; exit 1; }
  [[ -f "$SNV_R2"  ]] || { echo "Error: SNV read2 not found: $SNV_R2"; exit 1; }
fi

if [[ $RUN_CNV -eq 1 ]]; then
  [[ -z "$CNV_CSV" ]] && { echo "Error: --cnv-csv (or --snv-csv) is required for the CNV stage"; exit 1; }
  [[ -z "$CNV_R1"  ]] && { echo "Error: --cnv-r1 is required for the CNV stage"; exit 1; }
  [[ -z "$CNV_R2"  ]] && { echo "Error: --cnv-r2 is required for the CNV stage"; exit 1; }
  [[ -f "$CNV_CSV" ]] || { echo "Error: CNV CSV not found: $CNV_CSV"; exit 1; }
  [[ -f "$CNV_R1"  ]] || { echo "Error: CNV read1 not found: $CNV_R1"; exit 1; }
  [[ -f "$CNV_R2"  ]] || { echo "Error: CNV read2 not found: $CNV_R2"; exit 1; }
fi

echo "=============================================================="
echo " TETRIS-seq — single sample"
echo "   library    : $LIBRARY"
echo "   threads    : $THREADS"
echo "   stages     : SNV=$RUN_SNV  FLT3=$RUN_FLT3  CNV=$RUN_CNV"
[[ $need_snv_inputs -eq 1 ]] && echo "   SNV FASTQ  : $(basename "$SNV_R1") / $(basename "$SNV_R2")"
[[ $RUN_CNV -eq 1 ]]        && echo "   CNV FASTQ  : $(basename "$CNV_R1") / $(basename "$CNV_R2")"
echo "=============================================================="

# ---- 1. SNV panel ----------------------------------------------------------
if [[ $RUN_SNV -eq 1 ]]; then
  echo; echo ">> 1. SNV panel"
  "$SNV_PANEL" -s "$SNV_CSV" -o "$LIBRARY" -f "$SNV_R1" -g "$SNV_R2" -t "$THREADS"
fi

# ---- 2. FLT3-ITD (starts from the SNV panel's mapped merged BAM) -----------
if [[ $RUN_FLT3 -eq 1 ]]; then
  echo; echo ">> 2. FLT3-ITD (SSCS + DCS consensus, then Pindel)"
  # derive sample names from the SNV FASTQ exactly as the panel scripts do
  SAMPLE=$(python "$SAMPLE_SCRIPT" --samples_csv_file "$SNV_CSV" --fastq "$SNV_R1")
  SAMPLE_UDI=$(python "$SAMPLE_UDI_SCRIPT" --samples_csv_file "$SNV_CSV" --fastq "$SNV_R1")
  # the FLT3 arm needs consensus reads called WITHOUT SAM flag filtering, so it
  # cannot use the SNV panel's SSCS/DCS BAMs -- these are called separately
  "$FLT3_PREP"   -s "$SAMPLE" -u "$SAMPLE_UDI" -o "$LIBRARY" -t "$THREADS"
  "$FLT3_CALLER" -s "$SAMPLE" -u "$SAMPLE_UDI" -o "$LIBRARY" -t "$THREADS"
fi

# ---- 3. CNV / mCA panel ----------------------------------------------------
if [[ $RUN_CNV -eq 1 ]]; then
  echo; echo ">> 3. CNV / mCA panel"
  "$CNV_PANEL" -s "$CNV_CSV" -o "$LIBRARY" -f "$CNV_R1" -g "$CNV_R2" -t "$THREADS"
fi

# ---- manual next steps -----------------------------------------------------
cat << 'EOF'

==============================================================
 Automated per-sample stages complete for this sample.

 The MANUAL / CROSS-SAMPLE steps below are run once you have
 processed as many samples as you want (a whole lane, or fewer).
 See docs/running_the_pipeline.md and docs/pipeline_overview.md.

 SNV panel — across your batch of samples:
   * Fit the position-specific noise (error) model:
       noise_correction_model/Duplex_Error_Model_initial_variant_calling_v5.py
   * Post-model processing + trajectories (manuscript analysis, separate repo):
       https://github.com/the-blundell-lab/preAML_evolutionary_dynamics
   * Curate FLT3-ITD calls (require detection in BOTH SSCS and DCS).

 CNV / mCA panel:
   * Build the Panel of Normals (once, from the 36 control samples):
       mCA_caller/watson_code_create_PON.ipynb
   * Call mCAs against the PON:
       mCA_caller/watson_code_mCA_caller_v13_use_this_one.ipynb
   * For any individual with an mCA: sequence all timepoints, phase the index
     sample, then call at earlier timepoints (all within mCA_caller/):
       mCA_caller/watson_code_mCA_caller_phased_v5.ipynb
   * Chromosomal rearrangements: call on the raw BAM, then re-call any hit on
     the SSCS BAM for accurate VAF:
       chromosomal_rearrangement_caller/
   * Curate: KMT2A-PTD cell fraction 2*(R-1); mCA germline exclusion.
==============================================================
EOF
