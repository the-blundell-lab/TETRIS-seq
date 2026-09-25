# Demo

A small dataset for checking an installation against a known answer: a
commercial reference standard with certified variants at known allele fractions,
sequenced on the TETRIS-seq SNV/indel panel. Because the certified values come
from the manufacturer rather than from this pipeline, it tests the whole chain
from FASTQ to annotated variant calls.

## Getting the data

The FASTQ pair is 311 MB, too large for this repository, and is archived on
Zenodo:

**TETRIS-seq demo data**, <https://doi.org/10.5281/zenodo.XXXXXXX>

## What it contains

Horizon Discovery Myeloid DNA Reference Standard **HD829**, sequenced as
`Myeloid_100A` (SLX-20548, xGenUDI1). The full capture is 12.4 GB; this is a
locus-restricted subset built by [`make_horizon_demo.sh`](make_horizon_demo.sh),
which keeps every read aligning within 150 bp of the 14 certified SNV positions
or across the FLT3-ITD region, and discards the rest.

Selection is by locus, not by downsampling. Every read of a UMI family shares a
start position, so a family is either kept whole or dropped entirely, and duplex
consensus calling behaves at these positions exactly as it does on the full
library. Random downsampling would fragment families and destroy the duplex
yield.

## Running it

```bash
source scripts/activate.sh

run_sample.sh -o SLX_20548 -t 16 \
    --snv-csv /path/to/demo/horizon_samples.csv \
    --snv-r1  /path/to/demo/SLX-20548.xGenUDI1.DEMO.s_2.r_1.fq.gz \
    --snv-r2  /path/to/demo/SLX-20548.xGenUDI1.DEMO.s_2.r_2.fq.gz \
    --no-cnv
```

Use absolute paths, or paths relative to the working directory: the sample name
is read from the second dot-separated field of the FASTQ filename, so a leading
`../` prevents the sample being identified.

`--no-cnv` because HD829 was sequenced on the SNV/indel panel only.

## Expected output

In `Horizon_HD829_demo_xGenUDI1/DCS_files_MUFS3/Horizon_HD829_demo_SNV_watson_code_DCS_variants_MUFs_3_annotated.txt`,
the 14 certified SNVs, with the VAFs measured from the full (non-subsetted)
capture for comparison:

| gene | variant | position (hg19) | depth | variant | VAF % | full capture |
|---|---|---|---|---|---|---|
| KRAS | p.G13D | 12:25398281 C>T | 2399 | 894 | 37.266 | 37.266 |
| RUNX1 | p.M240I | 21:36206711 C>T | 2407 | 819 | 34.026 | 34.026 |
| NRAS | p.Q61L | 1:115256529 T>A | 2230 | 191 | 8.565 | 8.565 |
| DNMT3A | p.R882C | 2:25457243 G>A | 3014 | 160 | 5.309 | 5.309 |
| IDH1 | p.R132C | 2:209113113 G>A | 2154 | 109 | 5.060 | 5.063 |
| TET2 | p.R1261H | 4:106164914 G>A | 2267 | 111 | 4.896 | 4.896 |
| ASXL1 | p.W796C | 20:31022903 G>T | 3151 | 152 | 4.824 | 4.824 |
| TP53 | p.S241F | 17:7577559 G>A | 3028 | 146 | 4.822 | 4.822 |
| CBL | p.S403F | 11:119148988 C>T | 2280 | 107 | 4.693 | 4.693 |
| SF3B1 | p.G740E | 2:198266713 C>T | 1976 | 91 | 4.605 | 4.605 |
| FLT3 | p.D835Y | 13:28592642 C>A | 2752 | 120 | 4.360 | 4.360 |
| JAK2 | p.V617F | 9:5073770 G>T | 2019 | 88 | 4.359 | 4.359 |
| IDH2 | p.R172K | 15:90631838 C>T | 2617 | 114 | 4.356 | 4.356 |
| EZH2 | p.R413Q | 7:148514471 C>T | 2141 | 86 | 4.017 | 4.017 |

Thirteen reproduce the full capture exactly; IDH1 differs by one reference read
(2154 against 2153), which shifts its VAF by 0.002 percentage points.

[`check_horizon_demo.py`](check_horizon_demo.py) prints this comparison
automatically:

```bash
python demo/check_horizon_demo.py \
    Horizon_HD829_demo_xGenUDI1/DCS_files_MUFS3/Horizon_HD829_demo_SNV_watson_code_DCS_variants_MUFs_3_annotated.txt
```

And in `Horizon_HD829_demo_xGenUDI1/FLT3_calling/Pindel/Horizon_HD829_demo_FLT3_ITD_calls_DCS.csv`,
the certified FLT3 internal tandem duplication:

```
chr13:28,608,046-28,608,330    DUP:TANDEM
300 bp (284 bp duplication + 16 bp inserted sequence)
depth 2946   variant 120   VAF 4.07%   allelic ratio 0.042
```

## Run time and resources

About **40 minutes** using 16 threads (measured: 42 minutes), with a peak of
about **12.5 GB** of memory during consensus calling and indel realignment. A
machine with 16 GB will be tight; 32 GB is comfortable. Most of the time is in
SSCS and DCS consensus calling.

## What this demo does not cover

- **The *in silico* noise-correction model.** It fits a per-position error
  distribution across roughly 40 samples from a sequencing lane, so it cannot
  run on a single sample. The demo therefore reproduces depths and VAFs, not
  the model's p-values or final call set.
- **Panel-wide metrics.** Coverage and depth summaries across the panel will
  look degenerate, because reads exist only at the demo's target regions.
- **The mCA/chromosomal rearrangement panel.** That is a separate capture with
  its own FASTQ files, and HD829 carries no mosaic chromosomal alteration or
  rearrangement to detect. See
  [`../mCA_caller/`](../mCA_caller) and
  [`../chromosomal_rearrangement_caller/`](../chromosomal_rearrangement_caller).
