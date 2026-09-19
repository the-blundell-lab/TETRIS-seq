# Pipeline schematic

This diagram shows the TETRIS-seq workflow across both capture panels.
Green boxes are automated (run per sequencing lane / library by the `scripts/*.sh`
entry points). Amber boxes are downstream steps that you trigger manually,
either because they operate across a group of samples or because they require
curation. It renders natively on GitHub.

```mermaid
flowchart TD
    FASTQ["Raw paired-end FASTQ — SNV library<br/>(~40 samples per lane)"]
    FASTQ_C["Raw paired-end FASTQ — mCA / rearrangement library<br/>(separate capture, separate library)"]

    %% ---------------- SNV PANEL ----------------
    subgraph SNVP["SNV panel — Watson_code_SNV_panel_v1.7.sh"]
        direction TB
        S_NAME["Resolve sample + library name<br/>from the UDI index in the FASTQ name"]
        S_F2B["FASTQ &rarr; unmapped BAM<br/>Picard FastqToSam"]
        S_UMI["Move the 3 bp inline UMI into the RX tag<br/>fgbio ExtractUmisFromBam, read structure 3M2S146T<br/>(3 bp UMI, 2 bp spacer dropped, 146 bp template)"]
        S_ADPT["Mark adapter read-through in the XT tag<br/>Picard MarkIlluminaAdapters"]
        S_ALIGN["Align to hg19, keeping the tags<br/>SamToFastq (adapters masked) | bwa mem -p | MergeBamAlignment<br/>alignments merged back onto the unmapped BAM, so RX survives"]
        S_QC["QC: insert-size metrics + raw-BAM read coverage"]
        S_SSCS["Group reads into UMI families<br/>single-strand consensus (SSCS)<br/>&ge;3 reads, &ge;90% base agreement"]
        S_DCS["Pair complementary strands<br/>duplex consensus (DCS)<br/>disagreement between strands = N"]
        S_POST["Re-process each consensus level separately (SSCS and DCS)<br/>unmap &rarr; SortSam queryname &rarr; MarkIlluminaAdapters &rarr;<br/>SamToFastq | bwa mem | MergeBamAlignment &rarr;<br/>fgbio ClipBam (hard, 3 bp ends + overlaps) &rarr; GATK indel realignment"]
        S_CALL["SNV/indel calling<br/>samtools mpileup (-BOa -d 1000000 -Q0) &rarr; custom caller<br/>+ VarDictJava, run twice (with and without -k 0 local realignment)<br/>+ ANNOVAR annotation of each output"]
        S_NAME --> S_F2B --> S_UMI --> S_ADPT --> S_ALIGN
        S_ALIGN --> S_QC
        S_ALIGN --> S_SSCS
        S_SSCS --> S_DCS
        S_SSCS --> S_POST
        S_DCS --> S_POST
        S_POST --> S_CALL
    end

    subgraph FLT3["FLT3-ITD — Watson_code_SNV_panel_for_FLT3_calling_v1.1.sh, then Watson_code_Pindel_FLT3_caller.sh"]
        direction TB
        F_CONS["Re-call SSCS + DCS from the mapped merged BAM<br/>same grouping key, no SAM flag filter<br/>re-aligned, not clipped or indel-realigned"]
        F_PINDEL["Pindel on SSCS + DCS"]
        F_CONS --> F_PINDEL
    end

    %% ---------------- CNV PANEL ----------------
    subgraph CNVP["mCA / rearrangement panel — Watson_code_CNV_panel_v2.3.sh"]
        direction TB
        C_PREP["Same front end as the SNV panel, on its own library:<br/>UDI name lookup &rarr; FastqToSam &rarr; ExtractUmisFromBam (3M2S146T)<br/>&rarr; MarkIlluminaAdapters"]
        C_ALIGN["Align to hg19, keeping the tags<br/>SamToFastq | bwa mem -p | MergeBamAlignment"]
        C_QC["QC: insert-size metrics"]
        C_CONS["Group reads into UMI families<br/>single-strand consensus (SSCS)<br/>&ge;1 read (default -m 1), &ge;90% base agreement<br/>no duplex step on this panel"]
        C_POST["Re-process the consensus BAM<br/>unmap &rarr; SortSam &rarr; MarkIlluminaAdapters &rarr; re-align &rarr;<br/>ClipBam + GATK indel realignment"]
        C_KMT2A["KMT2A-PTD: mean depth exon 3:27"]
        C_BAF["Per-sample read depths + BAFs<br/>(mpileup, SNP calling)"]
        C_BRK["Translocation breakpoint<br/>read coverage"]
        C_PREP --> C_ALIGN --> C_CONS --> C_POST
        C_ALIGN --> C_QC
        C_POST --> C_KMT2A
        C_POST --> C_BAF
        C_POST --> C_BRK
    end

    %% ---------------- MANUAL / CROSS-SAMPLE ----------------
    BRIDGE["Reorganise per-sample output<br/>into the results root<br/>docs/data_layout.md section 3"]
    NOISE["Position-specific noise model<br/>across the ~40 samples<br/>noise_correction_model/"]
    POSTMODEL["Post-model variant processing<br/>+ germline filtering<br/>noise_correction_model/"]
    DYNAMICS["Acquisition age + fitness inference,<br/>phylogenies, Muller plots<br/>manuscript analysis code"]
    SNV_FINAL["Final SNV + indel calls<br/>(post-model processing)"]

    PON["Build Panel of Normals<br/>(36 control samples)<br/>mCA_caller/create_PON"]
    MCA_INDEX["Call mCA on index sample<br/>+ haplotype phasing<br/>mCA_caller/"]
    MCA_EARLY["Call mCA at earlier timepoints (phased)<br/>mCA_caller/ (phased notebook)"]

    REARR["Chromosomal rearrangement calling<br/>on raw BAM (pre-SSCS)<br/>chromosomal_rearrangement_caller/"]
    REARR_SSCS["Build SSCS over the breakpoint regions<br/>(no SAM flag filter, --regions)<br/>re-call for accurate VAF (only if hit)"]
    CURATE["Manual curation:<br/>FLT3-ITD (SSCS & DCS), KMT2A-PTD cell<br/>fraction, mCA germline exclusion"]

    %% ---------------- EDGES ----------------
    FASTQ --> SNVP
    FASTQ_C --> CNVP
    S_ALIGN --> F_CONS
    S_CALL --> BRIDGE --> NOISE --> POSTMODEL --> SNV_FINAL
    POSTMODEL --> DYNAMICS
    F_PINDEL --> CURATE

    C_KMT2A --> CURATE
    C_BAF --> PON --> MCA_INDEX --> MCA_EARLY --> CURATE
    C_BRK --> REARR
    CNVP --> REARR --> REARR_SSCS --> CURATE

    %% ---------------- STYLES ----------------
    classDef auto fill:#d5efd5,stroke:#2e7d32,color:#1b3d1b;
    classDef manual fill:#ffe6b3,stroke:#cc8400,color:#4d3300;
    classDef script fill:#f5e9f6,stroke:#8a4b93,color:#5a2a61;
    classDef io fill:#e3eefc,stroke:#2c5aa0,color:#12203a;
    classDef qc fill:#f1f1ef,stroke:#8c8f93,color:#4d5155;

    class FASTQ,FASTQ_C,SNV_FINAL io;
    class S_QC,C_QC qc;
    class S_NAME,S_F2B,S_UMI,S_ADPT,S_ALIGN,S_SSCS,S_DCS,S_POST,S_CALL,F_CONS,F_PINDEL,C_PREP,C_ALIGN,C_CONS,C_POST,C_KMT2A,C_BAF,C_BRK auto;
    class NOISE,POSTMODEL,DYNAMICS,PON,MCA_INDEX,MCA_EARLY,REARR,REARR_SSCS script;
    class BRIDGE,CURATE manual;
```

## How to read it

1. The two panels are **separate captures sequenced as separate libraries**, so each
   starts from its own FASTQ pair, run per lane by the `scripts/*.sh` entry points (green).
2. Both begin with the same front end: the sample and library names are read from the
   UDI index in the FASTQ filename, the reads are converted to an unaligned BAM, the 3 bp
   inline UMI on each read is moved into the RX tag, adapter read-through is marked, and
   the reads are aligned in a single `SamToFastq | bwa mem | MergeBamAlignment` pipe — the
   alignments are merged back onto the unaligned BAM so the UMI tags survive alignment.
   The same round trip is used again on each consensus BAM, since consensus reads are new
   sequences and cannot inherit their input's coordinates.
3. The green boxes are fully scripted: one command per panel produces
   consensus reads and the first-pass calls for every sample on the lane.
4. The amber boxes are triggered by you, after the automated run:
   - the noise model needs the whole group of ~40 samples together;
   - the PON is built once from the 36 controls, then mCAs are called and,
     where present, phased on the index sample and re-called at earlier
     timepoints, all within `mCA_caller/`;
   - chromosomal rearrangements are called from the raw BAM in
     `chromosomal_rearrangement_caller/` (the panel script only computes
     breakpoint read coverage), then re-quantified on the SSCS BAM if a hit was
     found;
   - final curation (FLT3-ITD requiring both SSCS & DCS, KMT2A-PTD cell
     fraction, mCA germline exclusion) is done by hand, as described in the manuscript.

See [`pipeline_overview.md`](pipeline_overview.md) for the per-step scripts and
[`data_layout.md`](data_layout.md) for input/output directory structure.
