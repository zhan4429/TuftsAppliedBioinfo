# Reading Materials

Applied Bioinformatics - Tufts University Department of Biology

These readings support the current labs and homework for Weeks 1–4 and 10.
Start with the named excerpts under **Before class**, aiming for about 20–30 minutes
per session. Use **For homework** as a reference while working on an assignment.
**Optional background** provides more detail when you want to explore a topic.

For research papers, read the specified sections and figures first. You do not need
to understand every algorithm or run a second workflow from a tutorial. Use the
commands and cluster settings in the course labs for your analysis.

## The Study Behind Homeworks 2 and 3

Jia et al. (2023), [A high-quality reference genome for the fission yeast Schizosaccharomyces osmophilus](https://pmc.ncbi.nlm.nih.gov/articles/PMC10085805/).

Read the abstract and the sequencing, genome assembly, and gene annotation methods.
This is the study behind the data you download and clean in Homework 2, then assemble
and annotate in Homework 3. Use the published genome as biological context for your
results. The classroom analysis uses a different workflow, so identical results are
not expected.

## Week 1

### Week 1, Session 1

**Command line, text processing, and Git**

- **Before class:** Software Carpentry's [Pipes and Filters](https://swcarpentry.github.io/shell-novice/04-pipefilter.html). Focus on redirecting output and connecting commands with pipes. These ideas explain how the lab counts and summarizes annotation features.
- **Before class:** Software Carpentry's [Tracking Changes with Git](https://swcarpentry.github.io/git-novice/04-changes.html). Follow the modify–stage–commit cycle and the use of `git status` and `git diff`.
- **For homework:** Keep the course [AWK cheat sheet](../tutorials/awk-cheatsheet.html) and [sed cheat sheet](../tutorials/sed-cheatsheet.html) open while working on Homework 1's feature queries and text replacements. Use the Git reading when creating your course-script repository and recording commits.

### Week 1, Session 2

**HPC, environments, and fetching public data**

- **Before class:** Tufts Research Technology's [Running your software](https://rtguides.it.tufts.edu/hpc/access/60-brief-slurm.html). Focus on compute resources, interactive sessions, and batch scripts.
- **For homework:** Nextflow training's [Manage workflow executions](https://training.nextflow.io/latest/nextflow_run/03_manage_executions/). Read the sections on execution history and cleaning work directories. They explain `nextflow log`, run names, and the relationship between cleanup and the ability to resume a run.

Connect these readings to Homework 1's job-script troubleshooting and after-class
cleanup. A submitted job, a completed job, and verified output are different checks.

## Week 2

### Week 2, Session 1

**Read quality across sequencing platforms**

- **Before class:** Data Carpentry's [Assessing Read Quality](https://datacarpentry.github.io/wrangling-genomics/02-quality-control.html). Read the FASTQ and quality-control explanations, especially how Phred scores relate to the FastQC plots.
- **For homework:** The [NanoPlot documentation](https://github.com/wdecoster/NanoPlot) provides a reference for long-read plots and report outputs. Compare read-length distributions, base yield, and quality rather than judging a long-read dataset from a short-read QC flag.
- **For homework:** In the [SeqKit statistics documentation](https://bioinf.shenwei.me/seqkit/usage/#stats), look up `num_seqs`, `sum_len`, N50, and `AvgQual`. Check what each number measures before using it in Homework 2's raw-read summaries and coverage calculations.

### Week 2, Session 2

**Trimming, filtering, and MultiQC**

- **Before class:** In the [fastp documentation](https://github.com/OpenGene/fastp#adapters), read the adapter-trimming and quality-filtering sections. Distinguish removing sequence from a read from discarding the whole read.
- **Before class:** MultiQC's [Using reports](https://docs.seqera.io/multiqc/reports/) explains how to navigate tables and plots and compare samples. Use the raw-versus-trimmed examples in the lab as your concrete case.
- **For homework:** Read the filtering section of the [fastplong documentation](https://github.com/OpenGene/fastplong#filtering) alongside the lab's explanation of mean quality and length thresholds.

For Homework 2, explain what preprocessing changed using retained reads and bases,
quality, and adapter content together. A higher retained percentage or more green
QC flags does not by itself establish a better result.

## Week 3

### Week 3, Session 1

**Genome assembly and evaluation**

- **Before class:** Li and Durbin (2024), [Genome assembly in the telomere-to-telomere era](https://www.genome.gov/sites/default/files/media/files/2024-10/Genome-assembly-in-the-telomere-to-telomere-era.pdf). Read the introduction, sequencing-data overview, and assembly-evaluation section. Focus on coverage, repeats, and the distinction between contiguity, completeness, and correctness.
- **For homework:** Rhie et al. (2020), [Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies](https://link.springer.com/article/10.1186/s13059-020-02134-9). Start with the abstract, Figure 1, and the explanations of QV and k-mer completeness. These support Homework 3's accuracy assessment.

Be ready to explain why a large N50 is useful but cannot establish base accuracy,
and why an accuracy estimate depends on the evidence used to calculate it.

### Week 3, Session 2

**Annotation and evidence**

- **Before class:** Gabriel et al. (2024), [BRAKER3: Fully automated genome annotation using RNA-seq and protein evidence with GeneMark-ETP, AUGUSTUS, and TSEBRA](https://pmc.ncbi.nlm.nih.gov/articles/PMC11216308/). Read the abstract, Figure 1, and its workflow explanation. Identify the inputs and what RNA and protein evidence contribute.
- **For homework:** Use the [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) as a lookup reference for feature types and the `ID` and `Parent` attributes. Focus on how genes, transcripts, exons, and CDS features are related.

Connect the readings to Homework 3's feature counts, intron-containing transcripts,
and RNA-seq support. An annotation feature count and the number of genes expressed
in one sampled condition are different quantities.

## Week 4

### Week 4, Session 1

**Alignment and alignment QC**

- **Before class:** In the official [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf), read the example in Section 1.1 and the mandatory fields in Section 1.4. Use the FLAG table and CIGAR definitions as references; the binary-format details are not needed for this lab.
- **For homework:** In the [mosdepth documentation](https://github.com/brentp/mosdepth), look up the output summaries, coverage thresholds, and mapping-quality filter. Relate these to the lab's comparison of all-read coverage with coverage supported by higher-MAPQ alignments.

For Homework 4, explain what a record and its coverage evidence mean, including
ambiguously placed reads, clipping, and supplementary alignments.

### Week 4, Session 2

**Variant calling and interpretation**

- **Before class:** GATK's [VCF—Variant Call Format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format). Focus on the header, REF/ALT, FILTER, INFO versus FORMAT, and genotype fields such as GT, AD, and DP.
- **For homework:** Krusche et al. (2019), [Best practices for benchmarking germline small-variant calls in human genomes](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699627/). Read the benchmarking recommendations and discussion of variant representations. Use the concepts to interpret normalization, caller agreement, and the limits of comparison without a truth set.
- **Optional background:** Morales et al. (2022), [A joint NCBI and EMBL-EBI transcript set for clinical genomics and research](https://pmc.ncbi.nlm.nih.gov/articles/PMC9007741/). Read the introduction and MANE rationale for the human transcript-choice example in Homework 4.

Keep the identity of the biological samples in mind. Differences between reads and
an assembly can have biological as well as technical explanations. Caller agreement
alone does not identify an assembly error.

## Week 10

### Week 10, Session 1

**Amplicon reads to an ASV table**

- **Before class:** In QIIME 2's [Gut-to-soil axis tutorial](https://amplicon-docs.qiime2.org/en/stable/tutorials/gut-to-soil/), read the sequence-quality-control and feature-table-construction section, including the paired-end DADA2 example and denoising statistics.
- **For the lab:** Relate truncation choices to the overlap needed for paired reads to merge. Follow the read counts through filtering, denoising, merging, and chimera removal, then connect the resulting ASVs to the feature table.
- **Optional background:** Gloor et al. (2017), [Microbiome Datasets Are Compositional: And This Is Not Optional](https://pmc.ncbi.nlm.nih.gov/articles/PMC5695134/). Read the introduction and Figure 1 for the lab's final relative-abundance exercise. Reuse this reading in Session 2.

### Week 10, Session 2

**Shotgun profiles, community statistics, and MAG quality**

- **Before class:** Read or revisit the introduction and Figure 1 of [Gloor et al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5695134/). Focus on why a change in relative abundance does not establish the same change in absolute abundance.
- **For the lab:** Read the dispersion discussion in the [`adonis2` documentation](https://vegandevs.github.io/vegan/reference/adonis.html), alongside the description of [`betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.html). Connect the two outputs to differences between groups and variation within groups.
- **For MAG assessment:** Bowers et al. (2017), [Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea](https://pmc.ncbi.nlm.nih.gov/articles/PMC6436528/). Read the quality categories and their requirements, including the RNA-gene requirements for high-quality drafts.
- **Optional background:** Nearing et al. (2022), [Microbiome differential abundance methods produce different results across 38 datasets](https://www.nature.com/articles/s41467-022-28034-z). Read the abstract and discussion for the motivation behind comparing methods and reporting sensitivity to analysis choices.

Use these readings to explain the statistical assumptions and limitations of a
result, including why completeness and contamination estimates alone do not establish
every requirement for a high-quality MAG.
