# Lab 3B — Annotation, and What Counts as Evidence

**Applied Bioinformatics · Week 3, Session 2**
**Runs alongside the lecture**

> You assembled a genome on Tuesday. Today you put genes on it — and, more importantly, you find out how much you should believe them.

---

## Before you start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED="$COURSE/shared/wk3"
export ASM=$MYWORK/week-03/session-01/asm.p_ctg.fa
export RNA="$MYWORK/week-01/session-02/fetch/results/fastq"

mkdir -p $MYWORK/week-03/session-02 && cd $MYWORK/week-03/session-02
ls -lh $ASM
ls -lh $RNA
```

Fallbacks if either is missing:

```bash
export ASM=$SHARED/fallback_assembly/asm.p_ctg.fa
export RNA="$SHARED/fallback_rna"
```

```bash
module load hisat2 samtools agat
module list
```

---

## Part 1 — Make the evidence, then submit the annotation

Structural annotation of a eukaryote is a **project, not a step**. A full BRAKER3 run on this genome takes hours, so it has to start early. It uses RNA-seq alignments as evidence, though, so align first; on a 12 Mb genome that is quick.

### 1.1 Index the assembly and align the RNA-seq

Use the Illumina RNA-seq run **SRR27956204**, already downloaded in Week 1 and
examined in Week 2. The same folder also contains paired-end BGISEQ genomic reads,
so a general `*_1.fastq.gz` search could select the wrong experiment. Use the exact
RNA-seq filenames below. If you use staged copies, check their names and adjust
`R1` and `R2` if needed.

```bash
mkdir -p logs
hisat2-build -p 8 $ASM asm_idx 2> logs/hisat2_build.log

R1="$RNA/SRX23612877_SRR27956204_1.fastq.gz"
R2="$RNA/SRX23612877_SRR27956204_2.fastq.gz"
ls -lh "$R1" "$R2"

# the first 4 million pairs are plenty of evidence for a 12 Mb genome
hisat2 -p 8 -u 4000000 -x asm_idx -1 $R1 -2 $R2 --summary-file rna_summary.txt \
  | samtools sort -@ 4 -o rna.sorted.bam -
samtools index rna.sorted.bam

cat rna_summary.txt
```

**Question 1.** What overall alignment rate did you get? What would a low rate mean here, and what does a *spliced* aligner do that `bwa` would not?

::: Answer
For RNA-seq against its own genome, expect 80–95%. Substantially lower means the RNA came from a different organism or strain, the assembly is poor, or there is heavy contamination.

HISAT2 is **splice-aware**: it can place one read across an intron, producing a CIGAR string containing `N`. A DNA aligner like `bwa` would either soft-clip the read at the junction or fail to place it — and you would lose exactly the reads that tell you where introns are.
:::

### 1.2 Submit BRAKER3

BRAKER3 annotates the same sequence the BAM was aligned to, so give it a copy of that assembly:

```bash
cp $ASM asm.fa
ls -lh $SHARED/orthodb_fungi.fa      # protein evidence, staged

cat > annotate.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=braker3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
module load braker3

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED="$COURSE/shared/wk3"
cd $MYWORK/week-03/session-02

# AUGUSTUS saves its trained model in its config directory, and the module's copy is read-only
[ -d augustus_config ] || cp -r "${AUGUSTUS_CONFIG_PATH:?the braker3 module did not set AUGUSTUS_CONFIG_PATH}" augustus_config

braker.pl \
    --genome=asm.fa \
    --bam=rna.sorted.bam \
    --prot_seq=$SHARED/orthodb_fungi.fa \
    --fungus \
    --species=wk3_${USER}_${SLURM_JOB_ID} \
    --AUGUSTUS_CONFIG_PATH=$PWD/augustus_config \
    --workingdir=$PWD/braker \
    --threads=${SLURM_CPUS_PER_TASK} \
    --gff3
EOF

bash -n annotate.sh
sbatch annotate.sh
squeue -u $USER
```

**Job ID: ______________**

Two choices in that script are deliberate. There is no repeat-masking step: yeast genomes carry few repeats, so we accept that risk here, but for a plant or animal you would mask first or BRAKER would annotate transposons as genes. And the species name includes your job ID, because AUGUSTUS refuses to overwrite a model that already exists. The result will be in `braker/braker.gff3` for your homework.

**Question 2.** Why is eukaryotic annotation a long job when the bacterial equivalent is quick?

::: Answer
Introns. A bacterial coding sequence is a continuous open reading frame, so you find ORFs and you are largely done. In a eukaryote the coding sequence is scattered across exons separated by tens of kilobases, so there *are* no long ORFs to find.

Add alternative splicing, untranslated regions with no coding signal at all, repeats that contain their own ORFs and get annotated as genes unless masked first, and pseudogenes that look exactly like genes — and you have a problem that needs evidence, training, and reconciliation rather than a scan.
:::

> **CHECKPOINT 1** — Everyone has a sorted, indexed RNA-seq BAM and a BRAKER3 job in the queue.

---

## Part 2 — Read the evidence

### 2.1 Find the introns

```bash
samtools view rna.sorted.bam | awk '$6 ~ /N/' | wc -l
samtools view rna.sorted.bam | awk '$6 ~ /N/ {print $6}' | grep -o '[0-9]*N' | sort | uniq -c | sort -rn | head
```

The first command counts spliced reads. The second lists the most common intron lengths. Both look only at column 6, the CIGAR string: searching the whole line would also count `N` bases in the read sequences.

**Question 3.** Those `N` operations are introns your data observed directly. Why is that different in kind from a gene finder predicting an intron?

::: Answer
A prediction says *this looks like a splice site, based on a model of what splice sites look like in this kind of organism.*

A spliced read says *a transcript was here, and it was spliced at this exact position.* One is a hypothesis; the other is an observation.

This is the whole distinction between *ab initio* and evidence-based annotation. Neither is sufficient alone — prediction finds genes that were not expressed, evidence confirms genes that were.
:::

### 2.2 Coverage over the assembly

```bash
samtools depth -a rna.sorted.bam | \
  awk '{t++; if($3>0) c++} END {printf "bases with any RNA coverage: %.1f%%\n", 100*c/t}'
```

**Question 4.** Only a fraction of the genome is covered. Is that a problem?

::: Answer
No — it is the expected result and an important one. Only transcribed regions get reads. Intergenic sequence, introns (mostly), and genes not expressed in this condition all sit at zero.

The corollary matters: **absence of RNA-seq coverage is not evidence that there is no gene there.** It may mean the gene was silent in the condition that was sequenced. This is why evidence alone cannot annotate a genome.
:::

---

## Part 3 — Judge gene models against evidence

A reference annotation for this organism is staged at `$SHARED/reference.gff3`. Before you compare it with your evidence, check that the two describe the same sequence.

### 3.1 Does this annotation describe your assembly?

`reference.gff3` was made on the staged assembly:

```bash
export REF_ASM=$SHARED/fallback_assembly/asm.p_ctg.fa
grep -m 3 '^>' $ASM $REF_ASM
md5sum $ASM $REF_ASM
```

**Question 5.** The contig names may look identical. Do the checksums agree? What would go wrong if you compared `reference.gff3` with a BAM made against your own assembly?

::: Answer
hifiasm names contigs by order (`ptg000001l`, `ptg000002l`, …), so two different assemblies share names. Unless you used the fallback assembly, the checksums differ: these are different sequences, and position 50,000 on one `ptg000001l` is not position 50,000 on the other.

Nothing would fail. `samtools bedcov` would happily report coverage over regions that do not correspond to the genes, giving plausible numbers and a wrong answer. This is the lecture's "seqid must match the FASTA" problem in its sharpest form: matching names are not enough; the sequence has to be the same.
:::

If the checksums differ, align the same RNA-seq to the staged assembly too:

```bash
if ! cmp -s $ASM $REF_ASM; then
  hisat2-build -p 8 $REF_ASM ref_idx 2> logs/hisat2_ref_build.log
  hisat2 -p 8 -u 4000000 -x ref_idx -1 $R1 -2 $R2 2> logs/hisat2_ref_align.log \
    | samtools sort -@ 4 -o ref.sorted.bam -
else
  ln -sf rna.sorted.bam ref.sorted.bam
fi
samtools index ref.sorted.bam
```

From here on, `ref.sorted.bam` and `reference.gff3` describe the same sequence.

### 3.2 Look at the annotation

```bash
cp $SHARED/reference.gff3 .
grep -v '^#' reference.gff3 | cut -f3 | sort | uniq -c | sort -rn
```

That is the Week 1 counting idiom, on a file you now understand.

```bash
awk -F'\t' '$3=="gene"' reference.gff3 | wc -l
awk -F'\t' '$3=="mRNA"' reference.gff3 | wc -l
awk -F'\t' '$3=="exon" {s+=$5-$4+1} END {print "exonic bases:", s}' reference.gff3
```

**Question 6.** There are more mRNA features than gene features. Why, and what does that mean for counting "how many genes are expressed"?

::: Answer
Alternative splicing — one gene locus produces several transcripts. The GFF3 hierarchy is `gene → mRNA → exon/CDS`, and a gene can parent many mRNAs.

For counting, it means you must decide what you are counting *before* you count. "2,000 transcripts detected" and "1,400 genes detected" can describe the same result. This is the transcript-choice problem from the lecture, and it is why MANE Select exists for human.
:::

### 3.3 How many genes have evidence?

```bash
awk -F'\t' '$3=="gene"' reference.gff3 | \
  awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5}' | sort -k1,1 -k2,2n > genes.bed

samtools bedcov genes.bed ref.sorted.bam | \
  awk '{n++; if($4>0) e++} END {printf "genes with RNA-seq support: %d / %d (%.1f%%)\n", e, n, 100*e/n}'
```

**Question 7.** You get some percentage below 100. Give two different explanations, and say how you would distinguish them.

::: Answer
1. **The gene is real but not expressed** in the condition sequenced. How common this is depends on the organism: in a compact yeast genome grown in one condition most genes have at least some reads, while one tissue of a multicellular organism detects far fewer.
2. **The gene model is wrong** — an over-prediction from an *ab initio* method, a repeat annotated as a gene, or a pseudogene.

To distinguish them: sequence more conditions, which raises the detection fraction for explanation 1 and leaves explanation 2 untouched. Also check whether the unsupported models cluster in repetitive regions, have unusual codon usage, or lack homologues — all signatures of over-prediction.

Note what you cannot do: conclude from one RNA-seq sample that a gene does not exist.
:::

### 3.4 One gene, read from the command line

Take the best-covered gene, list its annotated exons, and list the introns that reads actually crossed:

```bash
samtools bedcov genes.bed ref.sorted.bam | sort -k4,4nr | head -1 > top_gene.bed
read CHR START END COV < top_gene.bed
REGION="$CHR:$((START+1))-$END"
echo "$REGION"

# annotated exons, as start-end
awk -F'\t' -v c=$CHR -v s=$START -v e=$END \
  '$1==c && $3=="exon" && $4>s && $5<=e {print $4 "-" $5}' reference.gff3 | sort -u

# introns observed in reads, as start-end, with read counts
samtools view ref.sorted.bam "$REGION" | awk '{
  pos = $4; cig = $6
  while (match(cig, /[0-9]+[MIDNSHP=X]/)) {
    len = substr(cig, RSTART, RLENGTH-1); op = substr(cig, RSTART+RLENGTH-1, 1)
    if (op == "N") print pos "-" pos+len-1
    if (op ~ /[MDN=X]/) pos += len
    cig = substr(cig, RSTART+RLENGTH)
  }}' | sort | uniq -c | sort -rn | head
```

An annotated intron runs from one base after an exon ends to one base before the next exon starts. Most yeast genes have no intron, so if the junction list is empty, that agrees with a one-exon model: change `head -1` to `sed -n 2p`, `3p` and so on to try the next genes.

**Question 8.** Find a gene where reads cross at least one intron. Do the observed introns match the gaps between annotated exons? Describe one agreement or disagreement.

::: Answer
An intron supported by many reads at exactly the annotated coordinates is the model confirmed by observation. Disagreements worth describing:

- **A frequent junction where no intron is annotated** — an unannotated intron, or a different isoform
- **An annotated intron that no read crosses** — low expression, or a wrong model
- **A junction a few bases from the annotated boundary** — an alternative splice site, or the model is slightly off
- **A junction seen in one or two reads only** — weak evidence; read counts matter

Any of these, described accurately with the read counts, is a full answer. The skill is checking rather than trusting.
:::

> **CHECKPOINT 2** — Everyone has compared observed introns with annotated exon boundaries for at least one gene.

---

## Part 4 — Functional annotation, and what it is worth

Structural annotation says *where*. Functional annotation says *what it does* — and it is inference all the way down.

```bash
head -3 $SHARED/reference_functional.tsv
awk -F'\t' 'NR>1 && $2 ~ /hypothetical|unknown|uncharacterized/' $SHARED/reference_functional.tsv | wc -l
awk -F'\t' 'NR>1' $SHARED/reference_functional.tsv | wc -l
```

**Question 9.** A large fraction come back "hypothetical protein". Is that a failure of the annotation?

::: Answer
No. It means nobody has characterised anything similar enough to transfer a name from. It is an honest statement about the reference databases, not about your assembly or the organism.

For a well-studied fungus expect 20–35% hypothetical; for an under-studied lineage, far more. Lowering thresholds to reduce that number does not create knowledge — it creates confident wrong names, which are worse than honest silence.
:::

**Question 10.** Your structural annotation truncated a gene because of a frameshift in the assembly. What does the functional step do?

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
| Same sequence? | `md5sum assembly.fa annotated_assembly.fa` |
| Find observed introns | filter the BAM for CIGAR strings containing N |
| Coverage per feature | `samtools bedcov genes.bed rna.sorted.bam` |
| Annotation summary | `agat_sp_statistics.pl --gff file.gff3` |

**Five things to carry forward:**

1. A predicted intron is a hypothesis; a spliced read is an observation. Good annotation uses both, because each covers the other's blind spot.
2. Absence of RNA-seq coverage is not evidence of absence of a gene. It may just not have been expressed.
3. An annotation describes one exact sequence. Matching contig names are not enough; check that the assembly is the same.
4. "Hypothetical protein" is an honest answer about the databases, not a failure.
5. Functional annotation inherits every structural error, silently. Correctness upstream is not optional.
