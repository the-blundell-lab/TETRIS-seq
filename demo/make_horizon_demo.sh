#!/bin/bash
# Build a small, locus-restricted demo FASTQ pair from the full Horizon HD829
# capture (Myeloid_100A, SLX-20548 xGenUDI1).
#
# Every read at the 14 certified variant sites and across the FLT3-ITD region is kept and everything else is
# discarded, so UMI family structure is preserved exactly where the demo is
# checked. Reads are selected by locus, and all members of a UMI family share a
# start position, so families come out whole - unlike random downsampling,
# which collapses them.
#
#   ./make_horizon_demo.sh <R1.fq.gz> <R2.fq.gz> [outdir]
#
# Requires bwa, samtools and $EXTERNAL_TOOLS (or $REF) for the reference genome.

set -euo pipefail

R1="$1"; R2="$2"; OUT="${3:-Horizon_demo}"
REF="${REF:-${EXTERNAL_TOOLS:?set EXTERNAL_TOOLS or REF}/Homo_sapiens_assembly19.fasta}"
THREADS="${THREADS:-8}"
WINDOW="${WINDOW:-150}"

[ -f "$REF" ]     || { echo "reference not found: $REF"; exit 1; }
[ -f "$REF.bwt" ] || { echo "bwa index not found next to $REF"; exit 1; }

mkdir -p "$OUT"
BED="$OUT/horizon_demo_loci.bed"

# chromosome, 1-based start, 1-based end, label.  hg19 / b37 - bare chromosome
# names, to match Homo_sapiens_assembly19.fasta.  Point variants have start=end;
# the FLT3-ITD is a region, because the duplication is not at a fixed position.
# Each is padded by +/- $WINDOW bp.
while read -r chrom start end label; do
    printf '%s\t%s\t%s\t%s\n' "$chrom" "$((start - 1 - WINDOW))" "$((end + WINDOW))" "$label"
done > "$BED" <<'LOCI'
1 115256529 115256529 NRAS_p.Q61L
2 25457243 25457243 DNMT3A_p.R882C
2 198266713 198266713 SF3B1_p.G740E
2 209113113 209113113 IDH1_p.R132C
4 106164914 106164914 TET2_p.R1261H
7 148514471 148514471 EZH2_p.R413Q
9 5073770 5073770 JAK2_p.V617F
11 119148988 119148988 CBL_p.S403F
12 25398281 25398281 KRAS_p.G13D
13 28592642 28592642 FLT3_p.D835Y
13 28608024 28608351 FLT3_ITD_region
15 90631838 90631838 IDH2_p.R172K
17 7577559 7577559 TP53_p.S241F
20 31022903 31022903 ASXL1_p.W796C
21 36206711 36206711 RUNX1_p.M240I
LOCI

echo ">> 1/3  aligning, and keeping only the read names at the 15 target regions"
echo "        (the BAM is never written - it exists only in the pipe)"
echo "        bwa reports progress below as it goes; a full lane takes 1-2 hours"
# bwa's progress goes to stderr: keep a copy in bwa.log but also let it through,
# so it appears in whatever log this script is writing to.
bwa mem -t "$THREADS" "$REF" "$R1" "$R2" 2> >(tee "$OUT/bwa.log" >&2) \
  | samtools view -L "$BED" - \
  | cut -f1 \
  | sort -u > "$OUT/read_names.txt"

echo "        read names selected: $(wc -l < "$OUT/read_names.txt")"

echo ">> 2/3  pulling those read pairs out of the original FASTQs"
subset () {
    zcat "$1" | paste - - - - | awk -F'\t' -v N="$OUT/read_names.txt" '
        BEGIN { while ((getline n < N) > 0) keep[n] = 1 }
        { h = substr($1, 2); sub(/[ \t].*$/, "", h)
          if (h in keep) print $1 "\n" $2 "\n" $3 "\n" $4 }
    ' | gzip > "$2"
}
# The sample-name lookup takes the SECOND dot-field of the FASTQ name as the
# UDI index (pipeline_tools/Watson_code_sample_name_from_UDI_index_v1.py), so
# the demo files have to keep the <library>.<UDI>.<...>.r_N.fq.gz shape.
STEM="${STEM:-SLX-20548.xGenUDI1.DEMO.s_2}"

subset "$R1" "$OUT/$STEM.r_1.fq.gz"
subset "$R2" "$OUT/$STEM.r_2.fq.gz"

# matching sample sheet: sample name, then the UDI index as it appears above
printf '%s,%s\n' "${SAMPLE:-Horizon_HD829_demo}" "$(echo "$STEM" | cut -d. -f2)" \
    > "$OUT/horizon_samples.csv"

echo ">> 3/3  done"
ls -lh "$OUT/$STEM".r_*.fq.gz
echo "read pairs kept: $(( $(zcat "$OUT/$STEM.r_1.fq.gz" | wc -l) / 4 ))"
echo "sample sheet:    $OUT/horizon_samples.csv"
