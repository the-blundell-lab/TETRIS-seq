# mCA caller

Detection of **mosaic chromosomal alterations (mCAs)** from TETRIS-seq CNV-panel
data. Uses change-point detection (`ruptures`) on B-allele frequencies and
read-depth (LRR) to segment chromosomes and classify gains, losses and
copy-neutral LOH.

Runs in the analysis environment (`environment_analysis.yml`, Python 3.11);
Activate it with `conda activate tetris-seq-analysis`.
additional dependencies: `ruptures`, `seaborn`, `statsmodels`.

## Files

| File | Role |
|---|---|
| `watson_code_create_PON.ipynb` | **Builds the panel of normals (PON)** from QC-passing final-timepoint control samples and writes per-sample PON-normalised LRR (log2-ratio) files, with **leave-one-out** normalisation for samples that are themselves in the PON. Run this before mCA calling. |
| `watson_code_mCA_caller_v13.ipynb` | **Unphased mCA caller** (five-detector consensus) |
| `watson_code_mCA_caller_phased_v5.ipynb` | **Phased mCA caller** — uses the haplotype phasing of an index sample to call the same mCA at earlier timepoints |

## Running it on your own data

`watson_code_create_PON.ipynb` is written around the UKCTOCS cohort, and its
configuration cells name that cohort explicitly. To build a PON from your own
controls you must edit three things:

| setting | what it is | change it to |
|---|---|---|
| `libraries` | the sequencing library folders to look in, e.g. `['SLX_19285', 'SLX_20125', 'SLX_20127']` | your own library folder names |
| `final_timepoint_controls` | an explicit list of the control samples that go into the PON | your own control sample names |
| `if not sample_name.startswith("CNTRL")` | assumes control samples are named `CNTRL_*` | your own naming convention, or remove the filter |

Create its output directory before running — the notebook writes into it but
does not create it:

```bash
mkdir -p <data_directory>/CNV_panel_final_timepoint_read_depths/PON_normalised_read_depths
```

It expects the per-sample depth files laid out as

```
<library>/<sample>/CNV_read_depths/<sample>_sample_normalised_read_depths.txt
```

relative to the working directory — note the folder name and the file prefix
must match, which is not how the CNV panel writes them, so the per-sample output
needs reorganising first (see [`../docs/data_layout.md`](../docs/data_layout.md)).

The notebook does two things in sequence, each with its own configuration cell:

1. **Builds the PON** from the control libraries — the settings in the table above.
2. **Applies it**, writing per-sample PON-normalised log-R ratios for the samples
   you want to call mCAs in. The second configuration cell has its own
   `libraries` list (the samples to normalise, not the controls), and it only
   processes sample folders whose names begin `CNTRL` or `C92`:

   ```python
   samples_to_process = [s for s in all_subfolders if s.startswith('CNTRL') or s.startswith('C92')]
   ```

   Anything named otherwise is skipped without a message. It also reads the panel
   BED and `chromosome_ideogram_hg19.txt` by bare filename;
   `scripts/run_notebook.sh` links those in from `pipeline_tools/` for you.

If none of the configured samples are found the notebook does not fail: it
reports `Skipping missing library` and produces an empty panel of normals. If
your PON comes out empty, check these settings first.

Run the notebook against a data directory without moving it:

```bash
conda activate tetris-seq-analysis
scripts/run_notebook.sh watson_code_create_PON.ipynb /path/to/your/data
```

Run the PON notebook first, then the unphased caller; the phased caller is used
where an mCA is present at the final timepoint and you want to look for it in
that individual's earlier samples.

The simulation framework used to benchmark the callers' sensitivity and
specificity is not here — it belongs to the manuscript analysis and lives in the
[preAML_evolutionary_dynamics](https://github.com/the-blundell-lab/preAML_evolutionary_dynamics)
repository (Supplementary Figs. 26–30).
