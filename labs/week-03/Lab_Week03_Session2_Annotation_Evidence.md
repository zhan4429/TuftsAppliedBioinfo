# Lab 3B — Annotation, and What Counts as Evidence

**Applied Bioinformatics · Week 3, Session 2**
**Time: 30 minutes · Runs alongside the lecture**

> You assembled a genome on Tuesday. Today you put genes on it — and, more importantly, you find out how much you should believe them.

---

## Before you start

```bash
srun --pty -p interactive -n 8 --mem=32g --time=0-03:00:00 bash
hostname

export COURSE=/cluster/tufts/bio_appbio/shared/wk3
export ASM=~/appbio/wk3/assembly/asm.p_ctg.fa
export RNA=~/appbio/wk3/fetch_rna/results/fastq

mkdir -p ~/appbio/wk3/annot && cd ~/appbio/wk3/annot
ls -lh $ASM
ls -lh $RNA
```

Fallbacks if either is missing:

```bash
export ASM=$COURSE/fallback_assembly/asm.p_ctg.fa
export RNA=$COURSE/fallback_rna
```

```bash
module load anaconda
source activate appbio
conda install -y -c conda-forge -c bioconda hisat2 samtools agat
```

---

## Part 1 — Submit the annotation, then stop waiting for it (5 minutes)

Structural annotation of a eukaryote is a **project, not a step**. A full BRAKER3 run on this genome takes hours. Submit it now; it will be there for your homework.

```bash
mkdir -p logs

cat > annotate.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=annotate
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load singularity
# Instructor will give you the exact command for the tool we are using.
# It will be roughly:
#   singularity exec <braker3.sif> braker.pl \
#       --genome=asm.p_ctg.fa --bam=rna.sorted.bam \
#       --prot_seq=proteins.fa --softmasking --threads=16
EOF

echo "Annotation script written — the instructor will complete it with the site-specific command."
```

**Question 1.** Why is eukaryotic annotation hours of work when the bacterial equivalent is four minutes?

::: Answer
Introns. A bacterial coding sequence is a continuous open reading frame, so you find ORFs and you are largely done. In a eukaryote the coding sequence is scattered across exons separated by tens of kilobases, so there *are* no long ORFs to find.

Add alternative splicing, untranslated regions with no coding signal at all, repeats that contain their own ORFs and get annotated as genes unless masked first, and pseudogenes that look exactly like genes — and you have a problem that needs evidence, training, and reconciliation rather than a scan.
:::

**So what can you do in thirty minutes?** Produce the evidence yourself, and use it to judge gene models. That is today's lab, and it is the part that actually teaches you something.

---

## Part 2 — Make the evidence (11 minutes)

### 2.1 Index the assembly

```bash
hisat2-build -p 8 $ASM asm_idx 2> logs/hisat2_build.log
ls asm_idx*
```

### 2.2 Align the RNA-seq

```bash
R1=$(ls $RNA/*_1.fastq.gz | head -1)
R2=${R1/_1.fastq.gz/_2.fastq.gz}
echo "Aligning: $(basename $R1)"

hisat2 -p 8 -x asm_idx -1 $R1 -2 $R2 --summary-file rna_summary.txt \
  | samtools sort -@ 4 -o rna.sorted.bam -
samtools index rna.sorted.bam

cat rna_summary.txt
```

**Question 2.** What overall alignment rate did you get? What would a low rate mean here, and what does a *spliced* aligner do that `bwa` would not?

::: Answer
For RNA-seq against its own genome, expect 80–95%. Substantially lower means the RNA came from a different organism or strain, the assembly is poor, or there is heavy contamination.

HISAT2 is **splice-aware**: it can place one read across an intron, producing a CIGAR string containing `N`. A DNA aligner like `bwa` would either soft-clip the read at the junction or fail to place it — and you would lose exactly the reads that tell you where introns are.
:::

### 2.3 Find the introns

```bash
samtools view rna.sorted.bam | awk '$6 ~ /N/' | wc -l
samtools view rna.sorted.bam | grep -o '[0-9]*N' | sort -n | uniq -c | sort -rn | head
```

**Question 3.** Those `N` operations are introns your data observed directly. Why is that different in kind from a gene finder predicting an intron?

::: Answer
A prediction says *this looks like a splice site, based on a model of what splice sites look like in this kind of organism.*

A spliced read says *a transcript was here, and it was spliced at this exact position.* One is a hypothesis; the other is an observation.

This is the whole distinction between *ab initio* and evidence-based annotation. Neither is sufficient alone — prediction finds genes that were not expressed, evidence confirms genes that were.
:::

### 2.4 Coverage over the assembly

```bash
samtools depth -a rna.sorted.bam | \
  awk '{t++; if($3>0) c++} END {printf "bases with any RNA coverage: %.1f%%\n", 100*c/t}'
```

**Question 4.** Only a fraction of the genome is covered. Is that a problem?

::: Answer
No — it is the expected result and an important one. Only transcribed regions get reads. Intergenic sequence, introns (mostly), and genes not expressed in this condition all sit at zero.

The corollary matters: **absence of RNA-seq coverage is not evidence that there is no gene there.** It may mean the gene was silent in the tissue and condition that was sequenced. This is why evidence alone cannot annotate a genome.
:::

> **CHECKPOINT 1** — Everyone has a sorted, indexed RNA-seq BAM against their own assembly.

---

## Part 3 — Judge gene models against evidence (10 minutes)

A reference annotation for this organism is staged. Load it alongside your evidence and see how they compare.

### 3.1 Look at the annotation

```bash
cp $COURSE/reference.gff3 .
grep -v '^#' reference.gff3 | cut -f3 | sort | uniq -c | sort -rn
```

That is the Week 1 counting idiom, on a file you now understand.

```bash
awk -F'\t' '$3=="gene"' reference.gff3 | wc -l
awk -F'\t' '$3=="mRNA"' reference.gff3 | wc -l
awk -F'\t' '$3=="exon" {s+=$5-$4+1} END {print "exonic bases:", s}' reference.gff3
```

**Question 5.** There are more mRNA features than gene features. Why, and what does that mean for counting "how many genes are expressed"?

::: Answer
Alternative splicing — one gene locus produces several transcripts. The GFF3 hierarchy is `gene → mRNA → exon/CDS`, and a gene can parent many mRNAs.

For counting, it means you must decide what you are counting *before* you count. "2,000 transcripts detected" and "1,400 genes detected" can describe the same result. This is the transcript-choice problem from the lecture, and it is why MANE Select exists for human.
:::

### 3.2 How many genes have evidence?

```bash
awk -F'\t' '$3=="gene"' reference.gff3 | \
  awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5}' | sort -k1,1 -k2,2n > genes.bed

samtools bedcov genes.bed rna.sorted.bam | \
  awk '{n++; if($4>0) e++} END {printf "genes with RNA-seq support: %d / %d (%.1f%%)\n", e, n, 100*e/n}'
```

**Question 6.** You get some percentage well below 100. Give two different explanations, and say how you would distinguish them.

::: Answer
1. **The gene is real but not expressed** in the condition sequenced. Most likely, and expected — a typical single-condition RNA-seq experiment detects perhaps 50–70% of annotated genes.
2. **The gene model is wrong** — an over-prediction from an *ab initio* method, a repeat annotated as a gene, or a pseudogene.

To distinguish them: sequence more conditions or tissues, which raises the detection fraction for explanation 1 and leaves explanation 2 untouched. Also check whether the unsupported models cluster in repetitive regions, have unusual codon usage, or lack homologues — all signatures of over-prediction.

Note what you cannot do: conclude from one RNA-seq sample that a gene does not exist.
:::

### 3.3 Look at it

Open IGV through OnDemand (**Interactive Apps → IGV**, or Files → download and use IGV on your laptop). Load:

1. **Genome:** `asm.p_ctg.fa` (Genomes → Load Genome from File)
2. **Annotation:** `reference.gff3`
3. **Evidence:** `rna.sorted.bam`

Find a well-covered gene and zoom in.

**Question 7.** Looking at one gene with good coverage: do the RNA-seq reads agree with the exon boundaries? Find one place where they do not and describe it.

::: Answer
Common disagreements worth spotting:

- **Reads extending past the annotated gene end** — the UTR is longer than annotated. UTRs have no coding signal and are the first thing *ab initio* methods get wrong
- **A spliced read crossing where the model shows an exon** — an unannotated intron, or a different isoform
- **Coverage in an annotated intron** — retained intron, an unannotated exon, or genomic DNA contamination in the library
- **A gap in coverage inside an exon** — usually just low expression, sometimes a mis-assembly

Any of these, described accurately, is a full answer. The skill is looking rather than trusting.
:::

> **CHECKPOINT 2** — Everyone has IGV open with three tracks and has found one disagreement.

---

## Part 4 — Functional annotation, and what it is worth (4 minutes)

Structural annotation says *where*. Functional annotation says *what it does* — and it is inference all the way down.

```bash
head -3 $COURSE/reference_functional.tsv
awk -F'\t' 'NR>1 && $2 ~ /hypothetical|unknown|uncharacterized/' $COURSE/reference_functional.tsv | wc -l
awk -F'\t' 'NR>1' $COURSE/reference_functional.tsv | wc -l
```

**Question 8.** A large fraction come back "hypothetical protein". Is that a failure of the annotation?

::: Answer
No. It means nobody has characterised anything similar enough to transfer a name from. It is an honest statement about the reference databases, not about your assembly or the organism.

For a well-studied fungus expect 20–35% hypothetical; for an under-studied lineage, far more. Lowering thresholds to reduce that number does not create knowledge — it creates confident wrong names, which are worse than honest silence.
:::

**Question 9.** Your structural annotation truncated a gene because of a frameshift in the assembly. What does the functional step do?

::: Answer
It faithfully assigns a product name to the truncated protein, with no warning. You now have a confident functional annotation for a protein that does not exist.

Functional annotation **inherits structural errors**. This is why Tuesday's correctness question matters even when all you want is a gene list — and why an assembly with residual indel errors produces systematically shorter predicted proteins. Check mean CDS length against expectation; it is the cheapest diagnostic you have.
:::

> **CHECKPOINT 3** — Everyone can state why functional annotation cannot be more reliable than the structural annotation beneath it.

---

## Takeaways

| Task | Command |
|---|---|
| Index a genome for RNA | `hisat2-build -p 8 asm.fa asm_idx` |
| Spliced alignment | `hisat2 -x idx -1 R1 -2 R2 \| samtools sort -o out.bam` |
| Find observed introns | filter the BAM for CIGAR strings containing N |
| Coverage per feature | `samtools bedcov genes.bed rna.sorted.bam` |
| Annotation summary | `agat_sp_statistics.pl --gff file.gff3` |

**Four things to carry forward:**

1. A predicted intron is a hypothesis; a spliced read is an observation. Good annotation uses both, because each covers the other's blind spot.
2. Absence of RNA-seq coverage is not evidence of absence of a gene. It may just not have been expressed.
3. "Hypothetical protein" is an honest answer about the databases, not a failure.
4. Functional annotation inherits every structural error, silently. Correctness upstream is not optional.

