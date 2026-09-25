#!/usr/bin/env python3
"""
Check a Horizon HD829 demo run against the 14 certified variants.

Usage:
    check_horizon_demo.py <sample>_SNV_watson_code_DCS_variants_MUFs_3_annotated.txt

Prints the measured duplex (DCS) depth and VAF at each certified position
alongside the values reported in the manuscript, and the difference between
them.

The certified variants are the myeloid reference standard's own specification;
the "manuscript" column is what TETRIS-seq measured for Myeloid_100A in the
full sequencing run (Supplementary Fig. 8). A demo built by locus-restricted
subsetting keeps every UMI family at these positions, so the VAFs should agree
closely; depths agree only if no reads were lost in subsetting.
"""
import csv
import sys

# gene, chrom, hg19 position, ref, alt, manuscript VAF %, depth, variant depth
CERTIFIED = [
    ("KRAS",   "chr12",  25398281, "C", "T", 37.266, 2399, 894),
    ("RUNX1",  "chr21",  36206711, "C", "T", 34.026, 2407, 819),
    ("NRAS",   "chr1",  115256529, "T", "A",  8.565, 2230, 191),
    ("DNMT3A", "chr2",   25457243, "G", "A",  5.309, 3014, 160),
    ("IDH1",   "chr2",  209113113, "G", "A",  5.063, 2153, 109),
    ("TET2",   "chr4",  106164914, "G", "A",  4.896, 2267, 111),
    ("ASXL1",  "chr20",  31022903, "G", "T",  4.824, 3151, 152),
    ("TP53",   "chr17",   7577559, "G", "A",  4.822, 3028, 146),
    ("CBL",    "chr11", 119148988, "C", "T",  4.693, 2280, 107),
    ("SF3B1",  "chr2",  198266713, "C", "T",  4.605, 1976,  91),
    ("FLT3",   "chr13",  28592642, "C", "A",  4.360, 2752, 120),
    ("IDH2",   "chr15",  90631838, "C", "T",  4.356, 2617, 114),
    ("JAK2",   "chr9",    5073770, "G", "T",  4.359, 2019,  88),
    ("EZH2",   "chr7",  148514471, "C", "T",  4.017, 2141,  86),
]


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    calls = {}
    with open(sys.argv[1]) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["chromosome"], int(row["start"]), row["REF"], row["ALT"])
            calls[key] = row

    print(f"{'gene':8} {'position':20} {'change':6} "
          f"{'depth':>7} {'var':>6} {'VAF %':>8} | {'paper %':>8} {'diff':>7}")
    print("-" * 82)

    found = 0
    for gene, chrom, pos, ref, alt, paper_vaf, paper_dp, paper_var in CERTIFIED:
        row = calls.get((chrom, pos, ref, alt))
        if row is None:
            print(f"{gene:8} {chrom + ':' + str(pos):20} {ref + '>' + alt:6} "
                  f"{'-':>7} {'-':>6} {'NOT CALLED':>8} | {paper_vaf:8.3f} {'':>7}")
            continue
        found += 1
        dp = int(row["total_depth"])
        var = int(row["variant_depth"])
        vaf = float(row["VAF"]) * 100
        print(f"{gene:8} {chrom + ':' + str(pos):20} {ref + '>' + alt:6} "
              f"{dp:7d} {var:6d} {vaf:8.3f} | {paper_vaf:8.3f} {vaf - paper_vaf:+7.3f}")

    print("-" * 82)
    print(f"{found}/{len(CERTIFIED)} certified variants present in the demo output")
    if found == len(CERTIFIED):
        print("All 14 detected.")
    else:
        print("Missing variants above are a problem: check the demo was built "
              "from the loci BED, and that the run completed.")


if __name__ == "__main__":
    main()
