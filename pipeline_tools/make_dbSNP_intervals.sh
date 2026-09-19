#!/bin/bash
###############################################################################
# make_dbSNP_intervals.sh
#
# Builds the per-chromosome dbSNP interval files the SNV and CNV panel scripts
# annotate against:
#
#   UCSC_dbSNP153common_hg19_chr1.interval ... _chr22.interval, _chrX.interval
#
# Source: UCSC's dbSnp153Common bigBed for hg19 (dbSNP build 153, "common"
# subset). NCBI dbSNP is public domain; UCSC permits redistribution of its
# table dumps, so the output may be shared.
#
# Usage:
#   pipeline_tools/make_dbSNP_intervals.sh [output_directory]
#
# Default output directory: $EXTERNAL_TOOLS/dbSNP (or ./dbSNP if unset).
#
# Requires bigBedToBed:  conda install -c bioconda ucsc-bigbedtobed
# Needs ~1.4 GB to download and ~10 GB of scratch space while converting.
###############################################################################

set -euo pipefail

OUTDIR="${1:-${EXTERNAL_TOOLS:-.}/dbSNP}"
URL="https://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp153Common.bb"
BB="$OUTDIR/dbSnp153Common.bb"
BED="$OUTDIR/dbSnp153Common.bed"

command -v bigBedToBed >/dev/null 2>&1 || {
    echo "bigBedToBed not found. Install it with:" >&2
    echo "    conda install -c bioconda ucsc-bigbedtobed" >&2
    exit 1
}

mkdir -p "$OUTDIR"

echo "1/3  downloading $URL (1.4 GB)"
[ -f "$BB" ] || wget -O "$BB" "$URL"

echo "2/3  converting bigBed to BED (this takes a while)"
[ -f "$BED" ] || bigBedToBed "$BB" "$BED"

echo "3/3  splitting by chromosome into $OUTDIR"
awk -v outdir="$OUTDIR" '
BEGIN {
    FS = OFS = "\t"
    split("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX", wanted, " ")
    for (i in wanted) keep[wanted[i]] = 1
    header = "chrom" OFS "chromStart" OFS "chromEnd" OFS "name" OFS "ref" OFS "altCount" OFS "alts" OFS "shiftBases" OFS "freqSourceCount" OFS "minorAlleleFreq" OFS "majorAllele" OFS "minorAllele"
}
($1 in keep) {
    f = outdir "/UCSC_dbSNP153common_hg19_" $1 ".interval"
    if (!(f in started)) { print header > f; started[f] = 1 }
    print >> f
}
' "$BED"

echo
echo "Done. Files written to $OUTDIR:"
ls -1 "$OUTDIR"/UCSC_dbSNP153common_hg19_*.interval | wc -l | sed 's/^/  /'
echo "  (remove $BB and $BED to reclaim space)"
