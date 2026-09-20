# mCA caller

Detection of **mosaic chromosomal alterations (mCAs)** from TETRIS-seq CNV-panel
data. Uses change-point detection (`ruptures`) on B-allele frequencies and
read-depth (LRR) to segment chromosomes and classify gains, losses and
copy-neutral LOH.

Runs in the analysis environment (`environment_analysis.yml`, Python 3.11);
additional dependencies: `ruptures`, `seaborn`, `statsmodels`.

## Files

| File | Role |
|---|---|
| `watson_code_create_PON.ipynb` | **Builds the panel of normals (PON)** from QC-passing final-timepoint control samples and writes per-sample PON-normalised LRR (log2-ratio) files, with **leave-one-out** normalisation for samples that are themselves in the PON. Run this before mCA calling. |
| `watson_code_mCA_caller_v13.ipynb` | **Unphased mCA caller** (five-detector consensus) |
| `watson_code_mCA_caller_phased_v5.ipynb` | **Phased mCA caller** — uses the haplotype phasing of an index sample to call the same mCA at earlier timepoints |

Run the PON notebook first, then the unphased caller; the phased caller is used
where an mCA is present at the final timepoint and you want to look for it in
that individual's earlier samples.

The simulation framework used to benchmark the callers' sensitivity and
specificity is not here — it belongs to the manuscript analysis and lives in the
[preAML_evolutionary_dynamics](https://github.com/the-blundell-lab/preAML_evolutionary_dynamics)
repository (Supplementary Figs. 26–30).
