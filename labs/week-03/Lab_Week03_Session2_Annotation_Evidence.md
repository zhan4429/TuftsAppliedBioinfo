# Week 3 Session 2: Annotation and Evidence

Applied Bioinformatics - Tufts University Department of Biology

> You assembled a genome on Tuesday. Today you put genes on it — and, more importantly, you find out how much you should believe them.

> **In class:** submit one job that aligns the RNA-seq and then runs BRAKER3 (Part 1). It runs for hours, so everything after that is done after class: Part 2 needs only the alignment, which the job finishes first, and Part 3 needs BRAKER3 to have finished.

---

## Before You Start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname
```

When `hostname` shows a compute node, not a login node, continue in that shell:

```bash
export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED="$COURSE/shared/wk3"
export ASM=$MYWORK/week-03/session-01/hifiasm/asm.nuc.fa
export RNA="$MYWORK/week-01/session-02/fetch/results/fastq"

mkdir -p $MYWORK/week-03/session-02 && cd $MYWORK/week-03/session-02
ls -lh $ASM
ls -lh $RNA
```

Fallbacks if either is missing:

```bash
export ASM=$SHARED/fallback_assembly/asm.nuc.fa
export RNA="$SHARED/fallback_rna"
```

The two inputs:

- `$ASM` is the nuclear assembly you made in Lab 3A Section 3.2: the chromosomes, without
  the extra mitochondrial and rDNA copies. Annotating the full assembly would only add
  duplicate rDNA and mitochondrial genes.
- `$RNA` holds the Illumina RNA-seq run **SRR27956204** that you downloaded in Week 1 and
  inspected in Week 2. [Its ENA record](https://www.ebi.ac.uk/ena/browser/view/SRR27956204)
  identifies it as paired-end RNA-seq from *Saccharomyces cerevisiae*. Reuse those files;
  you do not need another download.

Keep the Week 1 `results/fastq` files. Cleaning a finished Nextflow run's `work/`
directory, as taught in Week 1, is different from deleting its published `results/`;
this lab still needs these reads.

```bash
module load samtools
module list
```

---

## Part 1: Submit The Evidence And The Annotation As One Job

Structural annotation of a eukaryote is a **project, not a step**. A full BRAKER3 run on this genome takes hours, so it has to start early. It uses RNA-seq alignments as evidence, so one job does both: it aligns the RNA-seq first, then runs BRAKER3.

**Question 1.** Why does RNA-seq help annotate a genome, and what can it not tell you?

<details>
<summary>Answer</summary>

RNA-seq reads provide evidence for transcribed regions. Reads that span splice
junctions help locate exon boundaries and introns.

A gene that was not expressed in the sampled condition may have no RNA-seq support.
Absence of coverage does not establish absence of a gene, so annotation also uses
sequence-based predictions and protein evidence.

</details>

**Why RNA-seq and proteins?** A gene finder working from DNA alone has only sequence
patterns to go on: start and stop codons, the short motifs at splice sites, and the codon
usage typical of genes. That is enough to guess where genes are, not to be sure. BRAKER3
therefore adds two independent kinds of evidence, and this lab gives it one of each
(Figure 1):

- **RNA-seq from this yeast** (`rna.sorted.bam`). Cells copy their active genes into
  RNA, so reads pile up on the exons of the genes that were switched on when the cells
  were harvested, and a read split across an intron marks both ends of that intron to
  the exact base. RNA-seq is the most precise evidence for a gene's structure, but it
  says nothing about genes that were switched off.
- **Proteins from many fungi** (`$SHARED/orthodb_fungi.fa`, from the OrthoDB
  database). A gene's protein usually stays similar across related species even after
  its DNA has changed, because changes that break a protein tend to be lost. So when a
  fungal protein aligns to your genome, it points to a gene there, whether or not that
  gene was expressed in your sample. Protein evidence is less exact at the edges:
  related proteins differ in length, and the alignment has to estimate where each exon
  starts and stops.

Each covers the other's blind spot. RNA-seq gives exact structure for the genes it sees;
proteins find conserved genes that RNA-seq missed. Genes supported by both are the ones
BRAKER3 trusts most.

![Two kinds of gene evidence on the same stretch of genome. Gene A is switched on: RNA-seq reads stack on its three exons and split reads mark its two introns. Gene B is switched off: it has no reads, but fungal protein alignments cover it. Both genes appear in the final gene models.](assets/annotation-evidence/rna-and-protein-evidence.png)

*Figure 1. Gene A was switched on, so RNA-seq reads cover its exons and split reads
mark its introns exactly. Gene B was switched off and has no reads, but proteins from
many fungi still align to it. The gene models combine both. Most* S. cerevisiae *genes
have a single exon; the figure shows introns because that is where the two kinds of
evidence differ most.*

### 1.1 Check The Inputs

Use the Illumina RNA-seq run **SRR27956204**, already downloaded in Week 1 and
examined in Week 2. The same folder also contains paired-end BGISEQ genomic reads,
so a general `*_1.fastq.gz` search could select the wrong experiment. Use the exact
RNA-seq filenames below.

```bash
mkdir -p logs
cp $ASM asm.fa

R1="$RNA/SRX23612877_SRR27956204_1.fastq.gz"
R2="$RNA/SRX23612877_SRR27956204_2.fastq.gz"
ls -lh asm.fa "$R1" "$R2"
ls -lh $SHARED/orthodb_fungi.fa      # protein evidence, staged
```

All four files should exist before you submit. The job in Section 1.2 finds the same two RNA-seq files itself, and uses the staged copies in `$SHARED/fallback_rna` if yours are missing. If the staged files have different names, change the `R1` and `R2` lines in the script to match.

### 1.2 Submit One Job: Align, Then Annotate

BRAKER3 annotates the same sequence the RNA-seq was aligned to, so the job aligns to `asm.fa`, the copy you just made, and then gives BRAKER3 that same file:

```bash
cat > annotate.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=braker3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load hisat2 samtools braker3/3.1.3
set -euo pipefail

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED="$COURSE/shared/wk3"
cd $MYWORK/week-03/session-02

# the Week 1 RNA-seq, or the staged copies if yours are missing
RNA=$MYWORK/week-01/session-02/fetch/results/fastq
[ -s "$RNA/SRX23612877_SRR27956204_1.fastq.gz" ] || RNA=$SHARED/fallback_rna
R1=$RNA/SRX23612877_SRR27956204_1.fastq.gz
R2=$RNA/SRX23612877_SRR27956204_2.fastq.gz

# Step 1: align every RNA-seq pair to the assembly
hisat2-build -p ${SLURM_CPUS_PER_TASK} asm.fa asm_idx \
  > logs/hisat2_build.log 2>&1
hisat2 -p ${SLURM_CPUS_PER_TASK} -x asm_idx -1 $R1 -2 $R2 \
  --summary-file rna_summary.txt \
  | samtools sort -@ 4 -T sort_tmp -o rna.sorted.bam -
samtools index rna.sorted.bam

# Step 2: annotate, with the RNA-seq and protein evidence
braker.pl \
    --genome=asm.fa \
    --bam=rna.sorted.bam \
    --prot_seq=$SHARED/orthodb_fungi.fa \
    --fungus \
    --species=wk3_${USER}_${SLURM_JOB_ID} \
    --workingdir=$PWD/braker \
    --threads=${SLURM_CPUS_PER_TASK} \
    --gff3
EOF

bash -n annotate.sh
sbatch annotate.sh
squeue --me
```

**Step 1** is the alignment. `hisat2-build` indexes `asm.fa`, `hisat2` aligns every read pair, `samtools sort` writes them as a sorted BAM, and `samtools index` indexes it. `rna_summary.txt` records how well the reads aligned. All the reads are used: a 12 Mb genome needs far fewer for its well-expressed genes, but more reads give weakly expressed genes evidence too, and inside the job they cost only time.

**Step 2** is BRAKER3, given the BAM that step 1 just made. `set -euo pipefail` stops the job at the first error, so BRAKER3 never runs on a missing or half-written BAM.

**What each option tells BRAKER3:**

| Option | Meaning |
|---|---|
| `--genome=asm.fa` | The sequence to annotate. It must be the same assembly the RNA-seq was aligned to, which is why you copied `$ASM`. |
| `--bam=rna.sorted.bam` | The RNA-seq evidence: the sorted alignments from step 1 of the same job. BRAKER3 builds transcripts from them and reads the introns from the split reads. |
| `--prot_seq=$SHARED/orthodb_fungi.fa` | The protein evidence: fungal proteins from OrthoDB, which BRAKER3 aligns to the genome. |
| `--fungus` | Use GeneMark's model for fungal introns, which looks for the branch-point signal that fungal introns carry near their end. |
| `--species=wk3_${USER}_${SLURM_JOB_ID}` | The name for the gene model AUGUSTUS trains on this genome. It must be new, because AUGUSTUS will not overwrite an existing model; the job ID makes every run's name unique. |
| `--workingdir=$PWD/braker` | The folder for all results and intermediate files. |
| `--threads=${SLURM_CPUS_PER_TASK}` | How many CPUs to use: the 16 requested from SLURM. |
| `--gff3` | Also write the final genes as GFF3 (`braker.gff3`), alongside the default GTF. |

Figure 2 shows where each option acts. BRAKER3 runs two gene finders, GeneMark-ETP and
AUGUSTUS, which both use your evidence; a third tool, TSEBRA, then keeps the
best-supported version of each gene.

![BRAKER3 workflow. Inputs on the left: the genome (--genome), RNA-seq alignments (--bam) and fungal proteins (--prot_seq). Step 1, GeneMark-ETP (--fungus), uses the evidence to find strongly supported genes and trains step 2, AUGUSTUS (--species). Step 3, TSEBRA, combines both gene sets. Outputs in the braker folder (--workingdir, --gff3): braker.gff3, braker.gtf, braker.aa, braker.codingseq and braker.log.](assets/annotation-evidence/braker3-workflow.png)

*Figure 2. What BRAKER3 does with your three inputs, with each `annotate.sh` option
shown at the step it controls.*

One choice in that script is deliberate: there is no repeat-masking step. Yeast genomes carry few repeats, so we accept that risk here, but for a plant or animal you would mask repeats first, or BRAKER3 would annotate transposons as genes. The results will be in `braker/`; the gene models are in `braker/braker.gff3`.

**Question 2.** Most *S. cerevisiae* genes have no intron. Why can this BRAKER3 annotation still take much longer than a simple bacterial gene prediction?

<details>
<summary>Answer</summary>

For this yeast, long introns and extensive alternative splicing are not the main explanation. Most genes are intronless, and its genome is compact.

BRAKER3 combines RNA and protein evidence with gene-finding models. It must align protein evidence, train gene predictors such as AUGUSTUS, and combine candidate models with the evidence. Those computations can take substantial time even on a small genome. Runtime also depends on the protein database, software settings and available CPUs.

Intron-rich plant and animal genomes add further complexity, but that architecture should not be used to describe this yeast.

</details>

**Checkpoint 1.** Everyone has one job in the queue that aligns the RNA-seq and then runs BRAKER3. That is all for class time.

---

## After Class

Everything from here on is done after class. Work on a compute node. If this is a new login, repeat the setup in Before You Start, including `module load samtools`.

Part 2 needs only the RNA-seq alignment, which the job finishes first: you can start it as soon as `rna.sorted.bam.bai` exists, even while BRAKER3 is still running. Part 3 needs BRAKER3 to have finished. Check the job's state with `sacct`, replacing the example ID with your `annotate.sh` job ID:

```bash
cd $MYWORK/week-03/session-02
sacct -j 12345 --format=JobID,JobName,State,ExitCode,Elapsed
```

`COMPLETED` with exit code `0:0` means `braker/braker.gff3` is ready. If the job failed, read the end of `logs/braker3_<job ID>.err`; if you cannot fix it, use the staged annotation described at the start of Part 3.

---

## Part 2: Read The Evidence

### 2.1 Alignment Rate And Coverage

Step 1 of your job wrote the alignment and its summary. Check that both are there:

```bash
ls -lh rna.sorted.bam rna.sorted.bam.bai
cat rna_summary.txt
```

**Question 3.** What overall alignment rate did you get? What would a low rate mean here, and what does a *spliced* aligner do that `bwa` would not?

<details>
<summary>Answer</summary>

For RNA-seq against its own genome, expect 80–95%. Substantially lower means the RNA came from a different organism or strain, the assembly is poor, or there is heavy contamination.

HISAT2 is **splice-aware**: it can place one read across an intron, producing a CIGAR string containing `N`. A DNA aligner like `bwa` would either soft-clip the read at the junction or fail to place it — and you would lose exactly the reads that tell you where introns are.

</details>

Now measure how much of the assembly the RNA-seq covers, at three depths: at least one read, at least 10, and at least 100:

```bash
samtools depth -a rna.sorted.bam | \
  awk '{t++; if ($3 > 0) a++; if ($3 >= 10) b++; if ($3 >= 100) c++}
       END {printf "1+ reads: %.1f%%   10+ reads: %.1f%%   100+ reads: %.1f%%\n", 100*a/t, 100*b/t, 100*c/t}'
```

**Question 4.** Almost all of the genome has at least one read. Why so much, in a genome where annotated features cover only about three-quarters of the sequence? And what does the fall at 10 and 100 reads tell you?

<details>
<summary>Answer</summary>

Three things add up:

- **Yeast genes sit almost end to end.** Annotated features cover about 74% of the reference genome, the gaps between genes are usually a few hundred bases, and few genes have introns, so there is little silent sequence to begin with.
- **Yeast transcribes more than its annotated genes.** The reference annotation has no untranslated regions, so the ends of every mRNA fall outside annotated features. The genome is also pervasively transcribed: noncoding, antisense and short unstable transcripts.
- **One read is a very low bar.** The job aligned every read pair, so even rare transcripts and background reach one read: unspliced pre-mRNA, traces of genomic DNA in the RNA preparation, and reads placed at random among copies of a repeat.

The fall at higher thresholds separates real expression from that background: well-expressed genes keep their coverage, while trace coverage drops away. So "any coverage" says little on its own.

The lesson that still holds is at the level of genes: **a gene with no reads may simply be silent in the condition that was sequenced**, so absence of RNA-seq coverage is not evidence that there is no gene there. This is why evidence alone cannot annotate a genome.

</details>

---

## Part 3: Your Gene Models And Their Proteins

Your BRAKER3 job must have finished before this part. Point one variable at the gene models and another at the assembly they were made on:

```bash
cd $MYWORK/week-03/session-02
GFF=braker/braker.gff3
FA=asm.fa
ls -lh $GFF $FA
```

If your BRAKER3 job failed, use the staged annotation instead. It was made on the staged assembly, so use that assembly with it; the copy lets the tools write an index file beside it:

```bash
GFF=$SHARED/fallback_annotation/braker.gff3
cp $SHARED/fallback_assembly/asm.nuc.fa fallback_asm.fa
FA=fallback_asm.fa
```

### 3.1 Look At Your Annotation

```bash
grep -v '^#' $GFF | cut -f3 | sort | uniq -c | sort -rn
```

That is the Week 1 counting idiom, on a file you made yourself.

```bash
awk -F'\t' '$3=="gene"' $GFF | wc -l
awk -F'\t' '$3=="mRNA" || $3=="transcript"' $GFF | wc -l
awk -F'\t' '$3=="CDS" {s+=$5-$4+1} END {print "coding bases:", s}' $GFF
```

**Question 5.** *S. cerevisiae* has about 6,000 protein-coding genes. How close is your gene count? Are there more transcripts than genes, and what would that mean?

<details>
<summary>Answer</summary>

A count in the same range as 6,000 is a good result for an automatic annotation. Many more genes suggests real genes split into pieces, or genes predicted in repeats; many fewer suggests genes that were missed or merged.

BRAKER3 can give one gene more than one transcript when the evidence supports alternative versions; their IDs then end in `.t1`, `.t2` and so on. Equal counts mean one transcript per gene, which is typical for yeast. Before calling extra transcripts alternative splicing, look at what actually differs between them.

</details>

**Checkpoint 2.** Everyone has gene and transcript counts for their annotation.

### 3.2 Extract The CDS And Protein Sequences

The GFF3 records *where* each gene's coding exons are; the sequences themselves are in the assembly. `gffread` joins each transcript's coding exons, reverse-complements genes on the minus strand, and translates them:

```bash
module load gffread
gffread -g $FA -x cds.fasta -y protein.fasta $GFF
grep -c '>' cds.fasta protein.fasta
head -4 protein.fasta
```

`-g` names the assembly the gene models were made on, `-x` writes the coding sequences, and `-y` writes their translations. Both files hold one sequence per transcript, named by its transcript ID.

Now look for proteins with a stop codon in the middle. `gffread` leaves out each protein's final stop codon and writes any other stop as `.`, so a `.` anywhere in a protein is an internal stop:

```bash
awk '/^>/ {if (s ~ /\./) n++; s = ""; next} {s = s $0}
     END {if (s ~ /\./) n++; print n + 0, "proteins with an internal stop codon"}' protein.fasta
```

**Question 6.** How many sequences does each file hold, and how does that compare with your transcript count from Section 3.1? How many proteins have an internal stop codon, and what could cause one?

<details>
<summary>Answer</summary>

Each file should hold one sequence per transcript, so both counts should equal the transcript count from Section 3.1. A different number means some models could not be read, or the GFF3 and the assembly do not belong together.

For BRAKER3's own models, expect no internal stops: its gene finders only predict open reading frames in the sequence they were given. So an internal stop here usually means `gffread` read a different assembly from the one BRAKER3 used, for example the staged annotation with your own assembly.

For gene models made another way, such as an annotation transferred from a related genome, an internal stop also reveals a problem in the sequence: a small insertion or deletion inside a gene shifts the reading frame of everything after it. Translating gene models is the cheapest check for these frameshift errors.

BRAKER3 also writes these sequences itself, as `braker/braker.codingseq` and `braker/braker.aa`. `gffread` works for any GFF3, and making them yourself checks that the annotation and the assembly fit together.

</details>

**Checkpoint 3.** Everyone has `cds.fasta` and `protein.fasta`, with one sequence per transcript.

---

## Takeaways

| Task | Command |
|---|---|
| Index a genome for RNA | `hisat2-build -p 8 asm.fa asm_idx` |
| Spliced alignment | `hisat2 -x idx -1 R1 -2 R2 \| samtools sort -o out.bam` |
| Coverage across the genome | `samtools depth -a rna.sorted.bam` |
| Annotation feature counts | `grep -v '^#' braker/braker.gff3` followed by `cut -f3`, `sort` and `uniq -c` |
| CDS and protein FASTA from a GFF3 | `gffread -g asm.fa -x cds.fasta -y protein.fasta braker.gff3` |

**Three things to carry forward:**

1. A predicted intron is a hypothesis; a spliced read is an observation. Good annotation uses both, because each covers the other's blind spot.
2. Absence of RNA-seq coverage is not evidence of absence of a gene. It may just not have been expressed.
3. A gene model is only useful with the assembly it was made on. Translating the models is the cheapest check that the two belong together.
