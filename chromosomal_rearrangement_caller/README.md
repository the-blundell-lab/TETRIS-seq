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
(+ `python-Levenshtein`), and `ZODB`/`BTrees`/`transaction` for the
`*_shelve_edit.py` variant only.

## Files

| File | Role |
|---|---|
| `Watson_code_translocation_calling_all_types_2025_v1_targeted.py` | **Main caller** (v1.0, July 2025). In: mapped/SSCS BAM, sample name, min mapping quality, output folder, panel-probe BED, read length. Out: CSV of rearrangements + read-type support, metrics file, filtered/grouped CSV + XLSX |
| `Watson_code_translocation_calling_all_types_2025_v1_targeted_shelve_edit.py` | As above, using on-disk ZODB storage (lower memory for large inputs) |
| `Watson_code_filter_and_group_chromosomal_rearrangements.py` | Filters and groups raw rearrangement calls (uses `networkx`) |
| `Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1_specific_regions.py` | SSCS consensus calling restricted to the rearrangement regions of interest |
| `Filtering translocation calls.ipynb` | Interactive filtering / review of translocation calls |
| `Translocation_regions_of_interest.bed` | Target regions of interest for rearrangement detection |

The `2025_v1` scripts are the current caller used in the manuscript. They
supersede an earlier Manta-based flow and the translocations-only consensus
caller (`Watson_code_SSCS_calling_for_translocations_1.1.py`); both have been
removed from this release, as neither was used for the manuscript and nothing
in the pipeline called them.

`Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1_specific_regions.py`
is the region-restricted form of the consensus caller that the FLT3-ITD arm also
uses — see [`../docs/pipeline_overview.md`](../docs/pipeline_overview.md) for how
the two arms relate.
