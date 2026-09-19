# mCA caller

Detection of **mosaic chromosomal alterations (mCAs)** from TETRIS-seq CNV-panel
data, plus the simulation framework used to benchmark the caller's sensitivity
and specificity. Uses change-point detection (`ruptures`) on B-allele
frequencies and read-depth (LRR) to segment chromosomes and classify gains,
losses and copy-neutral LOH.

Runs in the analysis environment (`environment_analysis.yml`, Python 3.11);
additional dependencies: `ruptures`, `seaborn`, `statsmodels`.

## Files

| File | Role |
|---|---|
| `watson_code_create_PON.ipynb` | **Builds the panel of normals (PON)** from QC-passing final-timepoint control samples and writes per-sample PON-normalised LRR (log2-ratio) files, with **leave-one-out** normalisation for samples that are themselves in the PON. Run this before mCA calling. |
| `watson_code_mCA_caller_v13_use_this_one.ipynb` | **Main mCA caller** (current version, v13) |
| `watson_code_mCA_caller_simulated_samples_v13_use_this_one.py` | Runs the caller on simulated samples and compares to ground truth (`all_detected_mCAs.csv`, `truth_comparison.csv`) |
| `watson_code_mCA_simulations_v7_training_and_test_sets.ipynb` | Generates the simulated training/test sample sets |
| `watson_code_mCA_caller_simulated_sample_sensitivity_analysis_v7.py` | Sensitivity / specificity / precision metrics and figures |
| `watson_code_mCA_caller_simulated_sample_sensitivity_analysis_phased_v7.py` | As above, for phased (haplotype-resolved) calling |
| `watson_code_mCA_caller_phased_v5.ipynb` | Phased mCA caller (earlier v5) |

Files marked **"use_this_one"** are the current/canonical versions. The `v5`
phased notebook and `v7` analysis scripts are retained for the phased-calling
and benchmarking workflows; confirm which you wish to publish before release.
