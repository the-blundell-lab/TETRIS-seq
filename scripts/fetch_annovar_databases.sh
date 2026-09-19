#!/bin/bash
###############################################################################
# fetch_annovar_databases.sh
#
# Obtains the six ANNOVAR (hg19) databases the SNV panel annotates against:
#
#     refGene  cosmic92_coding  cosmic92_noncoding  exac03  gnomad_genome  clinvar_20200316
#
#   scripts/fetch_annovar_databases.sh [annovar_directory] [--cosmic-dir DIR]
#
# Default directory: $ANNOVAR_HOME, else $EXTERNAL_TOOLS/annovar.
#
# Four of the six download automatically. The two COSMIC databases cannot be
# redistributed, so they are built from your own COSMIC v92 download - see the
# COSMIC section below, or run the script without --cosmic-dir and it will tell
# you exactly which files to fetch.
###############################################################################

set -euo pipefail

ANNOVAR_DIR=""
COSMIC_DIR=""
expect_cosmic_dir=0
for arg in "$@"; do
    if [ "$expect_cosmic_dir" = 1 ]; then COSMIC_DIR="$arg"; expect_cosmic_dir=0; continue; fi
    case "$arg" in
        --cosmic-dir) expect_cosmic_dir=1 ;;
        --cosmic-dir=*) COSMIC_DIR="${arg#*=}" ;;
        -h|--help) sed -n '3,17p' "$0"; exit 0 ;;
        *) ANNOVAR_DIR="$arg" ;;
    esac
done

ANNOVAR_DIR="${ANNOVAR_DIR:-${ANNOVAR_HOME:-${EXTERNAL_TOOLS:-$HOME/Pipeline_tools}/annovar}}"
HUMANDB="$ANNOVAR_DIR/humandb"
COSMIC_VERSION=92

if [ ! -f "$ANNOVAR_DIR/annotate_variation.pl" ]; then
    echo "ANNOVAR not found at $ANNOVAR_DIR"
    echo
    echo "ANNOVAR needs a (free, academic) registration before it can be downloaded:"
    echo "    https://www.openbioinformatics.org/annovar/annovar_download_form.php"
    echo "Unpack the tarball to $ANNOVAR_DIR, then run this script again."
    exit 1
fi

mkdir -p "$HUMANDB"
cd "$ANNOVAR_DIR"

###############################################################################
# 1. The four freely downloadable databases
###############################################################################
# refGene comes from ANNOVAR's own mirror; the rest from the annovar web source.
download_db () {
    local db="$1" source_flag="$2" marker="$3"
    if [ -s "$HUMANDB/$marker" ]; then
        echo "  [have] $db"
        return
    fi
    echo "  [get ] $db"
    # shellcheck disable=SC2086
    perl annotate_variation.pl -buildver hg19 -downdb $source_flag "$db" "$HUMANDB/"
}

echo "ANNOVAR databases -> $HUMANDB"
download_db refGene          "-webfrom annovar" hg19_refGene.txt
download_db exac03           "-webfrom annovar" hg19_exac03.txt
download_db gnomad_genome    "-webfrom annovar" hg19_gnomad_genome.txt
download_db clinvar_20200316 "-webfrom annovar" hg19_clinvar_20200316.txt

###############################################################################
# 2. COSMIC v92 - built locally from your own download
###############################################################################
coding_db="$HUMANDB/hg19_cosmic${COSMIC_VERSION}_coding.txt"
noncoding_db="$HUMANDB/hg19_cosmic${COSMIC_VERSION}_noncoding.txt"

if [ -s "$coding_db" ] && [ -s "$noncoding_db" ]; then
    echo "  [have] cosmic${COSMIC_VERSION}_coding, cosmic${COSMIC_VERSION}_noncoding"
    echo
    echo "All six databases are present."
    exit 0
fi

cosmic_explain () {
    cat <<EOF

COSMIC v${COSMIC_VERSION} (coding and non-coding) is licensed and cannot be redistributed,
so it is not downloaded for you. Register (free for academic use) at

    https://cancer.sanger.ac.uk/cosmic/register

then download the four GRCh37 v${COSMIC_VERSION} files into one directory:

    CosmicMutantExport.tsv.gz            CosmicCodingMuts.vcf.gz
    CosmicNCV.tsv.gz                     CosmicNonCodingVariants.vcf.gz

(Some releases name the VCFs CosmicCodingMuts.normal.vcf.gz and
CosmicNonCodingVariants.normal.vcf.gz; either name is accepted.) Gunzip them,
then re-run:

    scripts/fetch_annovar_databases.sh --cosmic-dir /path/to/those/files

Newer COSMIC releases work too, but the annotation columns in the published
variant call tables came from v${COSMIC_VERSION}; if you use another release, change the
-protocol names in scripts/Watson_code_SNV_panel_v1.7.sh to match.
EOF
}

if [ -z "$COSMIC_DIR" ]; then
    cosmic_explain
    exit 1
fi

if [ ! -d "$COSMIC_DIR" ]; then
    echo "--cosmic-dir $COSMIC_DIR does not exist"
    exit 1
fi

# Locate each input, allowing the .normal.vcf naming and gzipped copies.
find_cosmic () {
    local f
    for f in "$@"; do
        [ -s "$COSMIC_DIR/$f" ] && { echo "$COSMIC_DIR/$f"; return 0; }
    done
    return 1
}

mutant_export=$(find_cosmic CosmicMutantExport.tsv)                                   || true
coding_vcf=$(find_cosmic CosmicCodingMuts.vcf CosmicCodingMuts.normal.vcf)            || true
ncv_export=$(find_cosmic CosmicNCV.tsv)                                               || true
noncoding_vcf=$(find_cosmic CosmicNonCodingVariants.vcf CosmicNonCodingVariants.normal.vcf) || true

missing=0
for pair in "CosmicMutantExport.tsv:$mutant_export" "CosmicCodingMuts.vcf:$coding_vcf" \
            "CosmicNCV.tsv:$ncv_export" "CosmicNonCodingVariants.vcf:$noncoding_vcf"; do
    if [ -z "${pair#*:}" ]; then
        echo "missing in $COSMIC_DIR: ${pair%%:*}   (gunzip it if it is still .gz)"
        missing=1
    fi
done
if [ "$missing" = 1 ]; then
    cosmic_explain
    exit 1
fi

# prepare_annovar_user.pl ships with ANNOVAR, in the top directory.
prepare="$ANNOVAR_DIR/prepare_annovar_user.pl"
if [ ! -f "$prepare" ]; then
    echo "prepare_annovar_user.pl not found in $ANNOVAR_DIR"
    echo "It is part of the ANNOVAR download; re-unpack the tarball."
    exit 1
fi

if [ ! -s "$coding_db" ]; then
    echo "  [make] cosmic${COSMIC_VERSION}_coding   (this takes a while - CosmicMutantExport.tsv is ~15 GB)"
    perl "$prepare" -dbtype cosmic "$mutant_export" -vcf "$coding_vcf" > "$coding_db.tmp"
    mv "$coding_db.tmp" "$coding_db"
fi

if [ ! -s "$noncoding_db" ]; then
    echo "  [make] cosmic${COSMIC_VERSION}_noncoding"
    perl "$prepare" -dbtype cosmic "$ncv_export" -vcf "$noncoding_vcf" > "$noncoding_db.tmp"
    mv "$noncoding_db.tmp" "$noncoding_db"
fi

echo
echo "All six databases are present in $HUMANDB"
