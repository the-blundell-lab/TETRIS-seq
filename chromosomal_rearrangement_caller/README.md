# Chromosomal rearrangement caller

Custom caller for **chromosomal rearrangements / translocations** from TETRIS-seq
data (the mCA/translocation panel). It identifies split/discordant read support
for rearrangements at targeted regions of interest, then filters and groups the
calls.

**Workflow (see `../docs/pipeline_overview.md`):** rearrangement calling is first
run on the **mapped merged BAM** (pre-SSCS). If a rearrangement is called as
real, the caller is **re-run on the SSCS BAM** to obtain a more accurate VAF
estimate.

Runs in the analysis environment (`environment_analysis.yml`, Python 3.11).
Additional dependencies: `pysam`, `pyfaidx`, `networkx`, `fuzzywuzzy`
(+ `python-Levenshtein`), `xlsxwriter` and `openpyxl` for the grouped call
tables, and `ZODB`/`BTrees`/`transaction` for the `*_shelve_edit.py` variant
only.

ANNOVAR is used to annotate which genes the called breakpoints fall in. The
caller finds it via `$ANNOVAR_HOME`, or `$EXTERNAL_TOOLS/annovar` — so either
export one, or source the pipeline's config first:

```bash
source ../config/config.sh
```

## Files

| File | Role |
|---|---|
| `Watson_code_translocation_calling_all_types_2025_v1_targeted.py` | **Main caller** (v1.0, July 2025). In: mapped/SSCS BAM, sample name, min mapping quality, output folder, panel-probe BED, read length. Out: CSV of rearrangements + read-type support, metrics file, filtered/grouped CSV + XLSX |
| `Watson_code_translocation_calling_all_types_2025_v1_targeted_shelve_edit.py` | As above, using on-disk ZODB storage (lower memory for large inputs) |
| `Watson_code_filter_and_group_chromosomal_rearrangements.py` | Filters and groups raw rearrangement calls (uses `networkx`) |
| `Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1_specific_regions.py` | SSCS consensus calling restricted to the rearrangement regions of interest |
| `Filtering translocation calls.ipynb` | Interactive filtering / review of translocation calls |
| `Translocation_regions_of_interest.bed` | Target regions of interest for rearrangement detection |

## Worked example

Rearrangement calling is a **two-pass** process: call on the raw mapped BAM, then
re-call over just the breakpoint regions of any hit, using consensus reads, to
get an accurate VAF.

Set `ANNOVAR_HOME` first (or `source ../config/config.sh`), and run from a
working directory outside the repository.

### Pass 1 — call on the mapped merged BAM

```bash
P=/path/to/TETRIS-seq
mkdir -p TEMP <sample>_output

python $P/chromosomal_rearrangement_caller/Watson_code_translocation_calling_all_types_2025_v1_targeted.py \
    --infile <sample>_mapped_merged_bam.bam \
    --sample-name <sample> \
    --min-mapping-quality 20 --read-length 146 \
    --bed $P/pipeline_tools/TWIST_CNV_panel_TE-95031423_h19.bed \
    --targeted_bed $P/chromosomal_rearrangement_caller/Translocation_regions_of_interest.bed \
    --min-reads 5 --min-softclip-length 10 \
    --chromosomal_ideogram $P/pipeline_tools/chromosome_ideogram_hg19.txt \
    --ref /path/to/Homo_sapiens_assembly19.fasta \
    --out-directory <sample>_output
```

This annotates the breakpoints with ANNOVAR and then runs
`Watson_code_filter_and_group_chromosomal_rearrangements.py` itself. The table to
read is:

```
<sample>_translocations_found_just_those_specifically_targeted_both_sides_panel_grouped_and_filtered.csv
```

Each row is one breakpoint-supporting event, grouped where breakpoints fall
within 500 bp of each other. Rearrangements within a single gene are excluded,
so a sample with none between different genes produces an empty table.

### Pass 2 — re-call over the breakpoint regions with consensus reads

Pass 1 detects rearrangements in the raw reads; the VAF is then measured on
**consensus** reads over the breakpoint regions only. One command does it:

```bash
$P/scripts/rearrangement_SSCS_recall.sh \
    -i <sample>_mapped_merged_bam.bam \
    -s <sample> \
    -o <sample>_output \
    -r "21 36210000 36215000 8 93078000 93080000"
```

`-r` takes `chrom start end` triples — a window around the `LEFT COORDINATE` and
another around the `RIGHT COORDINATE` of the call in pass 1's
`*_grouped_and_filtered.csv`. The example above is the RUNX1 and RUNX1T1 pair of
a t(8;21) sample.

Note `-i` is the **same raw mapped merged BAM as pass 1**, not any output of it:
the consensus step regroups reads by UMI, so it needs the pre-consensus reads.

The wrapper runs three steps:

1. `Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1_specific_regions.py`
   over those regions, at the script's default consensus settings (minimum family
   size 1, 0.9 agreement threshold, base quality 10) — as used in the manuscript
2. `samtools sort` and `index` on the resulting SSCS BAM
3. the translocation caller again, with `--infile` pointing at the sorted SSCS BAM

Reference genome, panel BED, regions BED and ideogram come from
`config/config.sh` exactly as for the panel scripts, and `ANNOVAR_HOME` is
exported for the annotation step. Run the three steps by hand if you need to
vary the consensus settings — `-h` on the wrapper prints what it runs.

The VAFs from this consensus pass are the ones to report; pass 1 is for
detection only.

Review of the filtered call set is manual — see
[`Filtering translocation calls.ipynb`](Filtering%20translocation%20calls.ipynb).

`Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1_specific_regions.py`
is the region-restricted form of the consensus caller that the FLT3-ITD arm also
uses — see [`../docs/pipeline_overview.md`](../docs/pipeline_overview.md) for how
the two arms relate.
