# Demo

A small simulated dataset for checking that an installation works, and that the
chromosomal rearrangement caller reproduces a known answer.

## What this is

`BCR_ABL_demo.bam` is a simulated sample carrying a **t(9;22) BCR::ABL1**
rearrangement at **40% VAF in 100 cells**, generated as described in
Supplementary Note Box 3 of the manuscript: the BCR and ABL1 sequences were
taken from Ensembl (GRCh37 release 104), a breakpoint chosen at random from the
commonly affected introns (ABL1 intron 1–2, BCR intron 14–15), fusion sequences
built by joining the ends, and copies made in the proportions implied by the
cell number and VAF. Reads were then simulated to match the real panel —
fragmented to ~200 bp, "captured" with the actual panel probe sequences, and
written as paired 150 bp reads carrying 3 bp duplex UMIs — and processed to a
mapped SSCS BAM.

Because it is simulated, the truth is known exactly
(`BCR_ABL_demo_simulation_information.txt`):

```
number_cells = 100,  VAF = 0.4
BCR  breakpoint = chr22 23632754
ABL1 breakpoint = chr9  133717747
```

## Running it

Needs the reference genome and ANNOVAR (see [../docs/INSTALL.md](../docs/INSTALL.md)),
and both environments active:

```bash
source scripts/activate.sh
mkdir -p demo_output
python $TETRIS_SEQ/chromosomal_rearrangement_caller/Watson_code_translocation_calling_all_types_2025_v1_targeted.py \
    --infile demo/BCR_ABL_demo.bam \
    --sample-name BCR_ABL_demo \
    --min-mapping-quality 20 --read-length 150 \
    --bed $TETRIS_SEQ/pipeline_tools/TWIST_CNV_panel_TE-95031423_h19.bed \
    --targeted_bed $TETRIS_SEQ/chromosomal_rearrangement_caller/Translocation_regions_of_interest.bed \
    --min-reads 5 --min-softclip-length 10 \
    --chromosomal_ideogram $TETRIS_SEQ/pipeline_tools/chromosome_ideogram_hg19.txt \
    --ref $EXTERNAL_TOOLS/Homo_sapiens_assembly19.fasta \
    --out-directory demo_output
```

## Expected output

`demo_output/BCR_ABL_demo_translocations_found_just_those_specifically_targeted_both_sides_panel_grouped_and_filtered.csv`
should contain the reciprocal pair, notated **t(9;22)(q34.12;q11.23)**:

| left breakpoint | right breakpoint | supporting reads | VAF |
|---|---|---|---|
| chr9 133717746 (ABL1, intronic) | chr22 23632752 (BCR, intronic) | 120 | 0.32 / 0.30 |
| chr22 23632753 (BCR, intronic) | chr9 133717745 (ABL1, intronic) | 96 | 0.28 / 0.27 |

`BCR_ABL_demo_expected_output.csv` holds the full reference table.

Two things to note when comparing: the called breakpoints sit 1–2 bp from the
simulated ones, which is the caller's resolution; and the measured VAFs run a
little below the simulated 40%, which is the known slight underestimation for
these rearrangements described in the Supplementary Note.

## Expected run time

Under a minute on a desktop or a compute node, excluding the one-off reference
genome and ANNOVAR installation.
