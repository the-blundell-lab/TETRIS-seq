# Running the pipeline, in order

This is the run order from raw FASTQ to final calls. It separates the
automated, per-sample stages (green in the pipeline figure in the main
README) from the cross-sample and
longitudinal stages, which you trigger separately once you have processed as
many samples as you want. Most of those are scripts and notebooks in this
repository; only a handful of steps are genuinely manual (marked below).

Two things to keep in mind throughout:

- The SNV panel and the CNV/mCA panel are separate captures with separate
  FASTQ files. They are run independently and never share inputs.
- The wrapper processes one sample per call. Consensus calling is memory-hungry,
  so how many samples you can run concurrently depends on the memory available;
  a cluster job array is a convenient way to do a whole lane.

---

## Disk space and temporary files

The panel scripts create a `TEMP/` directory **in the directory you run them
from** and point everything heavy at it, rather than at the system `/tmp`:

| consumer | setting |
|---|---|
| Picard (sorting, merging, adapter marking) | `TMP_DIR=TEMP` |
| fgbio `ExtractUmisFromBam` | `--tmp-dir ./TEMP` |
| SSCS / DCS consensus calling | `shelve` database of UMI-grouped reads in `TEMP/` |

Consensus calling holds the UMI-grouped reads on disk rather than in memory, so
`TEMP/` can become large for a deeply sequenced sample. Run from a filesystem
with room for it (scratch space rather than a quota-limited home directory), and
keep the output directory there too.

Cleanup is only partial. The consensus caller removes its own `shelve` files
when it finishes, and Picard removes its spill files on success, but nothing
deletes `TEMP/` itself — and **a run that fails part-way leaves its temporary
files behind**. If you have had failures, check and clear it:

```bash
du -sh TEMP
rm -rf TEMP
```

Temporary files are named per sample, so several samples can safely share one
working directory.

---

## Stage A: automated, per sample

Run these for each sample. You can either call the three scripts yourself,
or use the wrapper that chains them in order.

### Option 1: the wrapper

The SNV/indel panel and the mCA/rearrangement panel are **separate captures,
sequenced as separate libraries**. Each has its own FASTQ pair, its own library
name and its own sample-index CSV, so each is a separate run:

The panel scripts write their output into `<sample>_<UDI>/` **relative to the
directory you run them from**, so work from a data directory outside the
repository and call the scripts by their path:

```bash
source /path/to/TETRIS-seq/scripts/activate.sh   # both environments, in order
cd /path/to/your/working/directory

# SNV / indel panel, and FLT3-ITD (which re-uses the SNV panel's mapped BAM)
/path/to/TETRIS-seq/scripts/run_sample.sh -o SLX_20125 -t 8 --no-cnv \
    --snv-csv snv_samples.csv --snv-r1 sampleA_SNV_R1.fq.gz --snv-r2 sampleA_SNV_R2.fq.gz

# mCA / rearrangement panel - different library, different FASTQ
/path/to/TETRIS-seq/scripts/run_sample.sh -o SLX_20124 -t 8 --no-snv --no-flt3 \
    --cnv-csv cnv_samples.csv --cnv-r1 sampleA_CNV_R1.fq.gz --cnv-r2 sampleA_CNV_R2.fq.gz
```

If both FASTQ pairs are available at once, a single call with all six arguments
and no `--no-*` flags runs all three stages in order:

```bash
/path/to/TETRIS-seq/scripts/run_sample.sh -o SLX_20125 -t 8 \
    --snv-csv snv_samples.csv --snv-r1 sampleA_SNV_R1.fq.gz --snv-r2 sampleA_SNV_R2.fq.gz \
    --cnv-csv cnv_samples.csv --cnv-r1 sampleA_CNV_R1.fq.gz --cnv-r2 sampleA_CNV_R2.fq.gz
```

`run_sample.sh` runs three steps in order: the SNV panel, then FLT3-ITD (which
re-calls its own consensus reads from the SNV panel's mapped merged BAM), then
the CNV panel. It prints the manual next steps at the end.
Use `--no-snv` / `--no-cnv` / `--no-flt3` to run a subset, and `-h` for full
help. To do a lane, loop the wrapper over samples (see its `-h` example).

### Option 2: the individual scripts

| Step | Command | Notes |
|---|---|---|
| 1. SNV panel | `./Watson_code_SNV_panel_v1.7.sh -s <snv_csv> -o <lib> -f <snv_R1> -g <snv_R2> -t 8` | Builds both SSCS and DCS; runs the custom caller + VarDictJava, then ANNOVAR. |
| 2a. FLT3 consensus | `./Watson_code_SNV_panel_for_FLT3_calling_v1.1.sh -s <SAMPLE> -u <sample_name_UDI> -o <lib> -t 8` | Calls SSCS and DCS again from the SNV panel's mapped merged BAM, without SAM flag filtering, into `<sample_name_UDI>/FLT3_calling/`. Required: the SNV panel's own SSCS/DCS BAMs discard the discordant and soft-clipped reads that carry an ITD. |
| 2b. FLT3-ITD | `./Watson_code_Pindel_FLT3_caller.sh -s <SAMPLE> -u <sample_name_UDI> -o <lib> -t 8` | Runs Pindel on the BAMs from 2a. `<SAMPLE>`/`<sample_name_UDI>` come from the SNV FASTQ + CSV. Use the same `-m` in 2a and 2b. |
| 3. CNV panel | `./Watson_code_CNV_panel_v2.3.sh -s <cnv_csv> -o <lib> -f <cnv_R1> -g <cnv_R2> -t 8` | Builds SSCS only; computes read depths, BAFs, KMT2A-PTD depth ratio, breakpoint coverage. |

Between Stage A and the SNV cross-sample steps, the per-sample output is
reorganised into the "results root" layout the notebooks expect. See
[`data_layout.md`](data_layout.md).

---

## Stage B: cross-sample and longitudinal (triggered separately)

Run these once you have processed the samples you want to analyse together.
Steps 1, 2 and 4-7 are scripts and notebooks in this repository; they are not
part of the per-sample driver, but they are not done by hand either. The steps that are
actually manual are marked **manual** below.

### SNV panel

1. Position-specific noise (error) model, fitted across the batch of samples
   (in the paper, the ~40 samples on a lane):
   `noise_correction_model/Duplex_Error_Model_initial_variant_calling_v5.py`
2. Post-model processing of the error-corrected SNVs and VarDictJava indels,
   and the trajectory figures, are part of the manuscript analysis rather than
   this pipeline: see [preAML_evolutionary_dynamics](https://github.com/the-blundell-lab/preAML_evolutionary_dynamics)
   (`Post_processing_variant_calls.ipynb` and the figure notebooks).
3. FLT3-ITD curation (**manual**): call an ITD real only if it is detected in
   both the SSCS and DCS Pindel outputs.

### CNV / mCA panel

4. Build the Panel of Normals (once), from QC-passing control samples (36 were used here):
   `mCA_caller/watson_code_create_PON.ipynb`
5. Call mCAs against the PON:
   `mCA_caller/watson_code_mCA_caller_v13.ipynb`
6. Longitudinal mCA calling. For any individual with an mCA, sequence all
   timepoints, phase the index sample, then call at earlier timepoints. All
   of this is in `mCA_caller/` (the phased notebook,
   `watson_code_mCA_caller_phased_v5.ipynb`).
7. Chromosomal rearrangements: call on the raw BAM; for any hit, build an SSCS
   BAM over the breakpoint regions only
   (`Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1_specific_regions.py`,
   `--regions` as chrom/start/end triples) and re-call on it for accurate VAF
   (`chromosomal_rearrangement_caller/Watson_code_translocation_calling_all_types_2025_v1_targeted.py`,
   then `Watson_code_filter_and_group_chromosomal_rearrangements.py`). Review of
   the filtered call set is **manual**
   (`chromosomal_rearrangement_caller/Filtering translocation calls.ipynb`).
8. Curation (**manual**): KMT2A-PTD cell fraction `2 × (R − 1)`; mCA germline
   exclusion. See [`pipeline_overview.md`](pipeline_overview.md#manual-curation-steps-not-automated-in-this-repository).

---

See the pipeline figure in the main README for the visual version of this flow,
and [`pipeline_overview.md`](pipeline_overview.md) for per-step script detail.
