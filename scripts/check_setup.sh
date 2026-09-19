#!/bin/bash
###############################################################################
# check_setup.sh
#
# Preflight check: reports which tools, JARs, reference files and panel
# resources the pipeline can find, using the same resolution rules as
# Watson_code_environment_setup_v2.6.sh (config/config.sh, then defaults).
#
# Run this before processing any samples:
#
#   scripts/check_setup.sh
#
# Exit status is 0 if everything required was found, 1 otherwise.
###############################################################################

P="$(cd "$(dirname "$0")" && pwd)"
CONFIG_FILE="${TETRIS_CONFIG:-$P/../config/config.sh}"
[ -f "$CONFIG_FILE" ] && source "$CONFIG_FILE" && echo "config: $CONFIG_FILE" \
                                               || echo "config: none found (using defaults)"

: "${EXTERNAL_TOOLS:=$HOME/Pipeline_tools}"
: "${REF:=$EXTERNAL_TOOLS/Homo_sapiens_assembly19.fasta}"
: "${ANNOVAR_HOME:=$EXTERNAL_TOOLS/annovar}"
: "${PINDEL_DIR:=$EXTERNAL_TOOLS/pindel}"
: "${PIPELINE_TOOLS:=$P/../pipeline_tools}"

missing=0
ok()   { printf '  [ ok ]      %-34s %s\n' "$1" "$2"; }
bad()  { printf '  [MISSING]   %-34s %s\n' "$1" "$2"; missing=$((missing+1)); }
note() { printf '  [ note ]    %-34s %s\n' "$1" "$2"; }

check_cmd()  { if command -v "$1" >/dev/null 2>&1; then ok "$1" "$(command -v "$1")"; else bad "$1" "not on PATH"; fi; }
check_file() { if [ -f "$2" ]; then ok "$1" "$2"; else bad "$1" "$2"; fi; }
check_dir()  { if [ -d "$2" ]; then ok "$1" "$2"; else bad "$1" "$2"; fi; }

echo
echo "Command-line tools"
for c in bwa samtools tabix vardict-java; do check_cmd "$c"; done

if command -v python >/dev/null 2>&1; then ok "python" "$(command -v python) ($(python --version 2>&1))"
elif command -v python3 >/dev/null 2>&1; then ok "python3" "$(command -v python3) ($(python3 --version 2>&1))"
else bad "python" "not on PATH"; fi

if command -v java >/dev/null 2>&1; then
    jv=$(java -version 2>&1 | head -1)
    case "$jv" in
        *\"1.8*|*\"8.*) ok "java" "$jv" ;;
        *) note "java" "$jv - Picard 2.18, GATK 3.8 and fgbio 1.3 expect Java 8" ;;
    esac
else bad "java" "not on PATH"; fi

echo
echo "Java archives (\$EXTERNAL_TOOLS = $EXTERNAL_TOOLS)"
if [ -f "$EXTERNAL_TOOLS/picard.jar" ]; then ok "picard.jar" "$EXTERNAL_TOOLS/picard.jar"
elif ls ${CONDA_PREFIX:-/nonexistent}/share/picard*/picard.jar >/dev/null 2>&1; then ok "picard.jar" "$(ls ${CONDA_PREFIX}/share/picard*/picard.jar | head -1) (conda)"
else bad "picard.jar" "$EXTERNAL_TOOLS/picard.jar"; fi
if [ -f "$EXTERNAL_TOOLS/fgbio-1.3.0.jar" ]; then ok "fgbio-1.3.0.jar" "$EXTERNAL_TOOLS/fgbio-1.3.0.jar"
elif ls ${CONDA_PREFIX:-/nonexistent}/share/fgbio*/fgbio.jar >/dev/null 2>&1; then ok "fgbio" "$(ls ${CONDA_PREFIX}/share/fgbio*/fgbio.jar | head -1) (conda)"
else bad "fgbio-1.3.0.jar" "$EXTERNAL_TOOLS/fgbio-1.3.0.jar"; fi
check_file "GenomeAnalysisTK.jar"  "$EXTERNAL_TOOLS/GenomeAnalysisTK.jar"

echo
echo "Other third-party tools"
check_file "ANNOVAR table_annovar" "$ANNOVAR_HOME/table_annovar.pl"
check_dir  "ANNOVAR humandb/"      "$ANNOVAR_HOME/humandb"
for db in refGene cosmic92_coding cosmic92_noncoding exac03 gnomad_genome clinvar_20200316; do
    if [ -s "$ANNOVAR_HOME/humandb/hg19_${db}.txt" ]; then ok "  $db" "hg19_${db}.txt"
    else bad "  $db" "$ANNOVAR_HOME/humandb/hg19_${db}.txt"; fi
done
for prog in pindel pindel2vcf; do
    if [ -x "$PINDEL_DIR/$prog" ]; then ok "$prog" "$PINDEL_DIR/$prog"
    elif command -v "$prog" >/dev/null 2>&1; then ok "$prog" "$(command -v "$prog")"
    else bad "$prog" "$PINDEL_DIR/$prog"; fi
done

echo
echo "Reference genome"
check_file "reference FASTA"       "$REF"
check_file "  .fai index"          "${REF}.fai"
check_file "  .dict"               "${REF%.fasta}.dict"
check_file "  BWA index (.bwt)"    "${REF}.bwt"
check_dir  "dbSNP intervals"       "$EXTERNAL_TOOLS/dbSNP"

echo
echo "Panel resources and custom scripts (\$PIPELINE_TOOLS = $PIPELINE_TOOLS)"
for f in Watson_code_SSCS_calling_2.1.py \
         Watson_code_DCS_calling_1.4.py \
         Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1.py \
         TWIST_SNV_panel-TE-92048328_hg19.bed \
         TWIST_CNV_panel_TE-95031423_h19.bed \
         TWIST_v2_transcript_details.csv \
         chromosome_ideogram_hg19.txt; do
    check_file "$f" "$PIPELINE_TOOLS/$f"
done

echo
if [ "$missing" -eq 0 ]; then
    echo "All checks passed."
    exit 0
fi

echo "$missing item(s) missing. To obtain them:"
echo

if [ ! -s "$REF" ] || [ ! -s "${REF}.bwt" ]; then
    echo "  reference genome (4 GB, prebuilt BWA index, from Zenodo):"
    echo "      scripts/fetch_reference.sh \"$EXTERNAL_TOOLS\""
    echo
fi
if [ ! -d "$EXTERNAL_TOOLS/dbSNP" ]; then
    echo "  dbSNP interval files (693 MB, from Zenodo):"
    echo "      pipeline_tools/make_dbSNP_intervals.sh \"$EXTERNAL_TOOLS/dbSNP\""
    echo
fi
if [ ! -f "$EXTERNAL_TOOLS/GenomeAnalysisTK.jar" ]; then
    echo "  GATK 3.8 (15 MB; licence prevents redistribution, so download it yourself):"
    echo "      wget -O - https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 | tar -xjf - -C \"$EXTERNAL_TOOLS\" --strip-components=1 --wildcards '*/GenomeAnalysisTK.jar'"
    echo
fi
if [ ! -f "$ANNOVAR_HOME/table_annovar.pl" ]; then
    echo "  ANNOVAR (registration required):"
    echo "      https://www.openbioinformatics.org/annovar/annovar_download_form.php"
    echo "      then unpack it to $ANNOVAR_HOME and run scripts/fetch_annovar_databases.sh"
    echo
fi
if [ ! -x "$PINDEL_DIR/pindel" ] && ! command -v pindel >/dev/null 2>&1; then
    echo "  Pindel (FLT3-ITD calling) - provided by the conda environment:"
    if [ -n "${CONDA_PREFIX:-}" ]; then
        echo "      the environment is active but does not have it; update it with"
        echo "      conda env update -n \"$(basename "$CONDA_PREFIX")\" -f environment_sequencing.yml"
    else
        echo "      conda activate tetris-seq-pipeline"
    fi
    echo
fi
echo "Full details: docs/INSTALL.md"
exit 1
