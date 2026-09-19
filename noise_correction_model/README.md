# In silico noise correction model

The position-specific beta-binomial error model that separates true
low-frequency variants from sequencing error, fitted per sequencing lane.

| File | Role |
|---|---|
| `Duplex_Error_Model_initial_variant_calling_v5.py` | Fits the model across a lane and writes the per-timepoint beta-binomial call files |
| `Checking_ASXL1_indels.ipynb` | QC check on *ASXL1* indel calls, a known artefact-prone region |
| `data_files/` | Panel site lists used by the model |

Runs in the analysis environment (`environment_analysis.yml`, Python 3.11).

## What happens next

Post-processing of these calls — rescuing variants called as errors at one
timepoint, germline filtering, removal of error trajectories, and the variant
call tables and trajectory figures built from them — is part of the manuscript
analysis rather than the pipeline, and is published in the
[preAML_evolutionary_dynamics](https://github.com/the-blundell-lab/preAML_evolutionary_dynamics)
repository (`Post_processing_variant_calls.ipynb` and the figure notebooks).
