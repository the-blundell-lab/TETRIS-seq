# mCA caller

Detection of **mosaic chromosomal alterations (mCAs)** from TETRIS-seq CNV-panel
data. Uses change-point detection (`ruptures`) on B-allele frequencies and
read-depth (LRR) to segment chromosomes and classify gains, losses and
copy-neutral LOH.

Runs in the analysis environment (`environment_analysis.yml`, Python 3.11);
Activate it with `conda activate tetris-seq-analysis`.
additional dependencies: `ruptures`, `seaborn`, `statsmodels`.

## Files

mCA calling is done by running three notebooks in order, not by a pipeline
script. Each is run manually and its output is the next one's input.

### Step 1 — build the panel of normals

[`watson_code_create_PON.ipynb`](watson_code_create_PON.ipynb)

Builds the PON from QC-passing final-timepoint control samples, then applies it
to produce per-sample PON-normalised read depths and log-R ratios, with
**leave-one-out** normalisation for samples that are themselves in the PON.

### Step 2 — call mCAs at the index timepoint

[`watson_code_mCA_caller_unphased.ipynb`](watson_code_mCA_caller_unphased.ipynb)

The unphased caller. For each sample it merges the per-SNP B-allele frequencies
with the PON-normalised log-R ratio, scans each chromosome with a sliding window
of heterozygous SNPs, and calls gains, losses and CN-LOH from the combined
BAF-deviation and LRR evidence. Run this at each individual's index (final)
timepoint.

### Step 3 — track mCAs back through earlier timepoints

[`watson_code_mCA_caller_phased.ipynb`](watson_code_mCA_caller_phased.ipynb)

The longitudinal phased caller. Once step 2 has called an mCA in an index
sample, the haplotype carrying it is phased from that sample and the same
phasing applied to every earlier timepoint from the same participant. Summing
the BAF deviation along a known haplotype is far more sensitive than looking for
an unphased BAF split, so an mCA can be traced back to timepoints where the
unphased caller sees nothing. The mean phased deviation is tested with a
one-sided t-test and converted to a cell fraction with a 95% confidence
interval.

### Getting from pipeline output to the callers' input

The pipeline leaves its output nested per sample, but steps 2 and 3 read
everything from **one flat directory**, set at the top of each notebook:

```python
CNV_DEPOSIT_DIR = 'Data_files/mCA_calling/Real_data/EGA_deposit_CNV_BAF_LRR'
```

They need two files per sample, which the pipeline writes in different places:

| file | written by |
|---|---|
| `<sample>_PON_normalised_read_depths_and_LRR.txt` | step 1 (the PON notebook), under `<library>/<sample>/PON_normalised_log2ratios_Feb2026/` |
| `<sample>_..._variant_calling_only_SNPs_annovar_annotated.txt` | the CNV panel, in the sample's output folder |

Collect them with:

```bash
scripts/collect_mCA_inputs.sh <your results root> <flat directory>
```

which searches for both kinds and symlinks them into one place, then prints the
`CNV_DEPOSIT_DIR` line to paste into the notebooks. It warns if either kind is
missing, since the callers need both. The same two notebooks generate Supplementary Figs. 31–38 and are
also published in the
[preAML_evolutionary_dynamics](https://github.com/the-blundell-lab/preAML_evolutionary_dynamics)
repository, where that default path points at the deposited data; keep the two
copies in step if you change either.

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
