# Data layout and end-to-end flow

This pipeline has **two stages** that use **different directory layouts**. Stage 2
(the noise model + notebooks) does **not** read directly from Stage 1's output;
the per-sample results are first reorganised into a curated "results root". This
page documents both layouts and the reorganisation step that links them, so the
workflow can be reproduced from raw FASTQ to final figures.

```
 FASTQ ──[Stage 1: scripts/*.sh]──► per-sample pipeline output
                                          │
                                          │  reorganise (see §3)
                                          ▼
                                    results root  ──[Stage 2: noise_correction_model/]──► calls + figures
```

---

## 1. Stage 1: pipeline output (what the `.sh` scripts write)

Each run of `Watson_code_SNV_panel_v1.7.sh` creates, **in the current working
directory**, one folder per sample named `<sample_name_UDI>` (sample name + UDI
index). `<m>` / `<d>` are the SSCS / DCS minimum family sizes (`-m` / `-a`,
default 3); `<SAMPLE>` is the sample abbreviation.

```
<sample_name_UDI>/
├── Metrics_and_images/
├── Gene_read_depths/
├── SSCS_files_MUFS<m>/
│   ├── <SAMPLE>_SNV_watson_code_SSCS_variants_MUFs_<m>_all_positions.vcf
│   └── <SAMPLE>_SNV_watson_code_SSCS_variants_MUFs_<m>_annotated.txt
├── DCS_files_MUFS<d>/
│   ├── <SAMPLE>_SNV_watson_code_DCS_variants_MUFs_<d>_all_positions.vcf
│   └── <SAMPLE>_SNV_watson_code_DCS_variants_MUFs_<d>_annotated.txt
├── FLT3_calling/                       ← written by the two FLT3 scripts
│   ├── SSCS_files_MUFS<m>/<SAMPLE>_SNV_watson_code_SSCS_mapped_merged.bam
│   ├── DCS_files_MUFS<m>/<SAMPLE>_SNV_watson_code_DCS_mapped_merged.bam
│   └── Pindel/
│       ├── <SAMPLE>_FLT3_ITD_calls_SSCS.csv
│       └── <SAMPLE>_FLT3_ITD_calls_DCS.csv
└── VarDictJava_MUFS<m>/
    ├── <SAMPLE>_SNV_watson_code_DCS_VarDictJava.txt
    └── <SAMPLE>_SNV_watson_code_DCS_VarDictJava_annotated.txt
```

Note: the `-o` option (`OUTFOLDER`) is a run/library label only. It is **not** used
to place the output; files are written under `<sample_name_UDI>/` in the
working directory.

## 2. Stage 2: analysis input (what the noise model + notebooks read)

`Duplex_Error_Model_initial_variant_calling_v5.py` and the notebooks use
**relative paths** and must be run from inside a "results root" directory
(in the manuscript this was `UKCTOCS_sample_level_results_V2/`) that is laid out
**per sample, by variant class and consensus type**:

```
<results root>/                         ← run the notebooks from here
├── <sample>/                           ← <sample> = the sample abbreviation (<SAMPLE> above)
│   └── SNV/
│       ├── DCS/
│       │   ├── <sample>_SNV_watson_code_DCS_variants_MUFs_3_all_positions.vcf
│       │   ├── <sample>_SNV_watson_code_DCS_variants_MUFs_3_annotated.txt
│       │   └── VarDictJava/<sample>_SNV_watson_code_DCS_VarDictJava_indel_variant_calls.txt
│       └── SSCS/ …
└── Data_files/                         ← reference tables + collected VCFs (sibling of the samples)
    ├── TWIST_SNV_panel_all_possible_changes_panel_sites.csv
    ├── All_SNV_panel_positions_seen_in_COSMIC_haem_lymphoid.csv
    └── VCF_files/DCS/<sample>_SNV_watson_code_DCS_variants_MUFs_3_all_positions.vcf
```

The model then writes its results (e.g.
`<sample>_SNV_watson_code_DCS_MUFs_3_beta_binomial_SNV_all_variant_calls_*.txt`)
back into `<results root>/<sample>/SNV/DCS/`, which the trajectory notebook reads.

In this repository the reference tables live in
[`../noise_correction_model/data_files/`](../noise_correction_model). Per-sample
model inputs and outputs are individual-level participant data and are not
distributed here; the columns are described above, and the calls themselves are
available via the EGA and in the paper's supplementary tables.

## 3. The reorganisation step (bridge between Stage 1 and Stage 2)

To go from §1 to §2, each sample's pipeline output is remapped as follows.
File names are unchanged; only the folders move or rename.

| Stage 1 (pipeline output)                | Stage 2 (results root)            |
|------------------------------------------|-----------------------------------|
| `<sample_name_UDI>/`                     | `<sample>/`                       |
| `<sample_name_UDI>/DCS_files_MUFS3/`     | `<sample>/SNV/DCS/`               |
| `<sample_name_UDI>/SSCS_files_MUFS3/`    | `<sample>/SNV/SSCS/`             |
| `<sample_name_UDI>/VarDictJava_MUFS3/`   | `<sample>/SNV/DCS/VarDictJava/`   |

Notes:
- `<sample>` (the abbreviation) is the filename prefix before `_SNV_watson_code`.
- Some collected VCFs are also gathered under `Data_files/VCF_files/DCS/`.
- The model tolerates an alternative `…_SNV_SNV_watson_code_…` double-prefix
  filename via a fallback, so both naming variants are accepted.

Reproducibility note: this reorganisation was performed as a curation step in
the original workflow and is **not** scripted in this repository. To reproduce
the analysis end to end, arrange the pipeline output into the Stage 2 layout
shown above before running the notebooks.
