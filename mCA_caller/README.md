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

### A known quirk in the control list

`watson_code_create_PON.ipynb` defines `final_timepoint_controls` as a literal
list. A missing comma after the first entry makes Python concatenate it with the
second:

```python
final_timepoint_controls = ['CNTRL_169_s7'      # <- no comma
                            'CNTRL_188_s7',
```

so that element becomes `'CNTRL_169_s7CNTRL_188_s7'`, matches no sample, and
both of those controls are excluded. The panel of normals behind the published
mCA calls was built this way.

**This is left as it is deliberately.** Adding the comma would change the panel
of normals, and therefore potentially the mCA calls, so the code no longer
reproduces the published results. If you are building a PON for your own data
you will be replacing this list anyway.

Run the PON notebook first, then the unphased caller; the phased caller is used
where an mCA is present at the final timepoint and you want to look for it in
that individual's earlier samples.

The simulation framework used to benchmark the callers' sensitivity and
specificity is not here — it belongs to the manuscript analysis and lives in the
[preAML_evolutionary_dynamics](https://github.com/the-blundell-lab/preAML_evolutionary_dynamics)
repository (Supplementary Figs. 26–30).
