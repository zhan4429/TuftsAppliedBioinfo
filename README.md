# Applied Bioinformatics

Tufts University - Department of Biology

This repository supports the Tufts Department of Biology course **Applied Bioinformatics**. The course emphasizes practical, real-world data analysis, reproducible workflows, and professional scientific communication. Students will work with real datasets, run bioinformatics tools and pipelines, organize projects with Git and GitHub, use notebooks for reporting, and practice critical evaluation of computational results.

Course policies, grading, due dates, and official requirements are governed by the current official syllabus. This repository is a working companion for labs; it is not a replacement for the syllabus.

## Course Website

This repository includes a static GitHub Pages site in `docs/`. After GitHub Pages is enabled for the repository, the default course site URL will be:

```text
https://zhan4429.github.io/TuftsAppliedBioinfo/
```

The deployment workflow is in `.github/workflows/pages.yml` and publishes the contents of `docs/` on pushes to `main`.

## What This Repository Is For

This repository is intended for in-class labs and related course materials. Weekly lab Markdown files should make it easy to:

- open the lab instructions in GitHub or a local editor
- copy commands during class
- follow the logic of each workflow step by step
- find sample sheets, small configuration files, and environment files
- keep a reproducible record of the tools, parameters, and outputs used in the course

Large datasets, large results directories, credentials, and private data should not be committed here unless the instructor explicitly says otherwise. Lab instructions should link to public data sources or explain how to obtain course-provided data.

## Getting Started

Clone the repository once:

```bash
git clone <repository-url>
cd <repository-name>
```

Before each class, pull the latest lab materials:

```bash
git pull
```

During lab:

1. Open the relevant weekly Markdown file in GitHub, a text editor, or your preferred notebook/workflow environment.
2. Read each section before copying commands.
3. Replace placeholders such as `<project-directory>`, `<sample-id>`, `<account>`, or `<partition>` with values provided in class.
4. Save your own notes, outputs, and reports in the location requested by the instructor. Do not overwrite shared course files unless asked.

## Suggested Repository Structure

The course repository may evolve, but a useful structure for weekly labs is:

```text
.
|-- README.md
|-- labs/
|   |-- week-00-prep/
|   |   `-- README.md
|   |-- week-01-command-line-hpc/
|   |   `-- README.md
|   |-- week-02-sequencing-qc/
|   |   `-- README.md
|   `-- ...
|-- envs/
|   `-- example-environment.yml
|-- containers/
|   `-- README.md
|-- configs/
|   `-- README.md
|-- samplesheets/
|   `-- README.md
`-- notebooks/
    `-- README.md
```

Each weekly lab should include the purpose of the lab, setup instructions, commands to run, expected outputs, interpretation questions, troubleshooting notes, and any homework or notebook/report expectations.

## Topic Overview

The syllabus organizes the course around applied workflows and reproducible computational practice.

| Week | Topics |
| --- | --- |
| Week 0 | Optional self-paced preparation in shell navigation and basic Python or R syntax. |
| Week 1 | Linux command line, text processing, Git/GitHub basics, HPC job submission with SLURM, conda/mamba environments, containers, public databases, and dataset inspection. |
| Week 2 | Sequencing technologies, experimental design, replication, confounding, batch effects, FASTQ format, FastQC, fastp, MultiQC, and QC interpretation. |
| Week 3 | Genome assembly and annotation, assembly QC, contamination screening, GFF3/GTF structure, functional annotation, annotation version issues, and IGV inspection. |
| Week 4 | Read alignment, SAM/BAM processing, alignment QC, variant calling with GATK HaplotypeCaller and DeepVariant, filtering, benchmarking, annotation, and VCF interpretation. |
| Week 5 | Bulk RNA-seq design, pseudoalignment, alignment-based quantification, library type and strandedness, nf-core/rnaseq, count matrices, metadata, and exploratory analysis. |
| Week 6 | Differential expression with DESeq2, design formulas, covariates and batch effects, multiple testing, visualization, enrichment analysis, and notebook/report preparation. |
| Week 7 | ATAC-seq, chromatin accessibility, peak calling with MACS3, nf-core/atacseq, ChIP-seq and CUT&RUN/CUT&Tag concepts, QC metrics, peak annotation, motif enrichment, and integration with expression results. |
| Week 8 | Single-cell RNA-seq preprocessing, 10x barcodes and UMIs, Cell Ranger, STARsolo, alevin-fry, nf-core/scrnaseq, AnnData and Seurat objects, QC, filtering, normalization, and feature selection. |
| Week 9 | Single-cell dimensionality reduction, clustering, marker genes, cell type annotation, batch integration, pseudobulk differential expression, and spatial transcriptomics. |
| Week 10 | Microbiome and metagenomics, 16S amplicon workflows, QIIME2/DADA2, compositional data, shotgun taxonomic and functional profiling, metagenome-assembled genomes, diversity metrics, and differential abundance. |
| Week 11 | Mini-project preparation, including dataset acquisition and QC as assigned. |
| Weeks 12-14 | Mini-project work: reproducing key findings from a published paper, maintaining a reproducible GitHub repository, preparing notebooks and reports, and presenting results. |

## Reproducibility Expectations

Bioinformatics work should be understandable and rerunnable by someone else, including your future self. For labs and projects:

- keep commands, parameters, sample IDs, and software versions visible
- use relative paths when possible
- separate raw data, intermediate files, final outputs, and reports
- document where input data came from
- preserve small configuration files and sample sheets
- avoid committing large generated files unless instructed
- record failed or unexpected results when they affect interpretation

The course emphasizes critical evaluation. A pipeline result is not automatically correct because it ran successfully. Compare tools when asked, check simple baselines, inspect QC reports, and explain limitations clearly.

## Tooling Expectations

**Git and GitHub:** Use Git from the beginning of the course. Pull updates before class, make clear commits for your own work, and use GitHub repositories to organize final project materials.

**HPC:** Many analyses will run on a high-performance computing system. Use the login node for light work only, submit compute jobs through the scheduler when needed, and follow course-specific guidance for accounts, partitions, modules, storage, and job scripts.

**conda and mamba:** Use environment files when provided. Prefer reproducible environment creation over ad hoc installation. If you add packages for an assignment or project, record the change.

**Containers:** Containers such as Docker images and Apptainer/Singularity images help make workflows portable. Record image names, tags, and any bind mounts or runtime options needed to rerun an analysis.

**Nextflow and nf-core:** Workflow managers are a reproducibility standard in research and industry. For Nextflow and nf-core labs, keep sample sheets, parameter files, configuration files, logs, and reports organized so that another student could understand how the run was produced.

**Notebooks:** Jupyter and R Markdown notebooks should communicate both computation and reasoning. A good notebook runs from top to bottom, explains major choices, labels outputs clearly, and includes interpretation rather than only code.

## Responsible AI Use

AI coding assistants are permitted and expected as part of the course, but they do not replace scientific judgment.

- Disclose AI use in submitted notebooks or reports, including the tool used and what it helped with.
- Verify generated commands and code against documentation, error messages, and expected outputs.
- You are responsible for every line you submit.
- Written interpretation, discussion, and limitations sections should be your own scientific reasoning.
- Watch for common AI failure modes: hallucinated arguments, deprecated APIs, plausible but inappropriate statistical choices, and code that runs while producing the wrong analysis.

## Final Project Reminder

The final project asks groups to reproduce key findings from a published paper using the original dataset. A strong project repository should include code, environment or container files, sample sheets, notebooks, and a README sufficient for another group to rerun the analysis. The final report should honestly describe what could and could not be reproduced and why.
