# Week 4 Session 1: Alignment QC

Applied Bioinformatics - Tufts University Department of Biology

> Last week you built a genome. This week it is your reference. You will align genomic reads labeled as the same yeast strain to it, inspect the alignments, and investigate differences and mapping problems.

> **Read Part 1 first and submit the alignment job immediately.** It runs while you do the rest of the lab.

---

## Before You Start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname                      # must NOT say login-p0x

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED=$COURSE/shared/wk4
export FQ=$MYWORK/week-01/session-02/fetch/results/fastq
export GENOME=12100000        # expected genome size in bases

mkdir -p $MYWORK/week-04/session-01/logs && cd $MYWORK/week-04/session-01

module load bwa-mem2 minimap2 samtools mosdepth seqkit
module list
```

Every tool is an environment module; `module avail <tool>` shows the versions installed.

### Dataset Identity Matters

The [BGISEQ BioSample](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN24477750) and
[HiFi BioSample](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN24477612) both identify
**S288C**, with collection dates 2021-01-21 and 2021-12-29 respectively. These are
separate BioSamples; the records do not establish one shared DNA preparation or
prove identical genomes. We use a haploid nuclear-genome model as a teaching
assumption and interpret disagreements as candidates to investigate, rather than
automatically correcting the assembly.

Gather the three inputs: your Week 3 assembly, the BGISEQ genomic reads from Week 1, and the full HiFi input linked in Lab 3A.

```bash
cp $MYWORK/week-03/session-01/asm.p_ctg.fa asm.fa
md5sum asm.fa > asm.md5

DNB1=$(ls $FQ/*SRR17374239*_1.fastq.gz 2>/dev/null)
DNB2=$(ls $FQ/*SRR17374239*_2.fastq.gz 2>/dev/null)
ln -sf $MYWORK/week-03/session-01/hifi.fastq.gz hifi.fq.gz

ls -lhL asm.fa "$DNB1" "$DNB2" hifi.fq.gz      # -L follows the link, so a missing file shows up
```

Fallbacks for anything missing:

```bash
cp $COURSE/shared/wk3/fallback_assembly/asm.p_ctg.fa asm.fa && md5sum asm.fa > asm.md5
DNB1=$(ls $COURSE/shared/wk2/fallback_fastq/*SRR17374239*_1.fastq.gz)
DNB2=$(ls $COURSE/shared/wk2/fallback_fastq/*SRR17374239*_2.fastq.gz)
rm -f hifi.fq.gz && seqkit head -n 30000 $(ls $COURSE/shared/wk3/fallback_hifi/*.fastq.gz | head -1) -o hifi.fq.gz
```

---

## Part 1: Submit The Alignment, Now

### 1.1 Take Two Million Read Pairs

```bash
seqkit head -n 2000000 "$DNB1" -o dnb_1.fq.gz
seqkit head -n 2000000 "$DNB2" -o dnb_2.fq.gz
seqkit stats -T dnb_1.fq.gz dnb_2.fq.gz | \
  awk -v g=$GENOME 'NR>1 {s+=$5} END {printf "coverage: %.0fx\n", s/g}'
```

**Question 1.** Week 1 fetched two short-read runs from this project. Why align the BGISEQ run and not the Illumina one?

<details>
<summary>Answer</summary>

The Illumina run is RNA-seq. Its reads come from transcripts, so their depth measures expression rather than copy number, intergenic sequence gets almost nothing, and reads spanning introns would be clipped by a DNA aligner. Every coverage-based diagnosis in this lab would be wrong.

The BGISEQ run is genomic DNA: every part of the genome is sampled roughly evenly. Choosing reads for an analysis starts with what the library measured, not with which platform made it. This is the Week 1 warning, arriving as a consequence.

</details>

### 1.2 Index The Reference Once

```bash
bwa-mem2 index asm.fa 2> logs/bwa_index.log
samtools faidx asm.fa
ls asm.fa*
```

On a 12 Mb genome this takes seconds. A human genome takes about an hour and tens of gigabytes, which is why shared references come with their indexes already built.

### 1.3 Submit Both Alignments

```bash
cat > align.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
module load bwa-mem2 samtools minimap2
set -euo pipefail
cd $MYWORK/week-04/session-01
T=${SLURM_CPUS_PER_TASK}
mkdir -p "tmp/$SLURM_JOB_ID"

# short reads: align, add mate tags, sort, mark duplicates, index
bwa-mem2 mem -t $T -R '@RG\tID:SRR17374239\tSM:SRR17374239\tLB:lib1\tPL:DNBSEQ' \
    asm.fa dnb_1.fq.gz dnb_2.fq.gz \
  | samtools fixmate -m -@ 2 - - \
  | samtools sort -@ 2 -m 1G -T "tmp/$SLURM_JOB_ID/sort" - \
  | samtools markdup -@ 2 - sample.bam
samtools index sample.bam

# HiFi reads from Week 3, same genome
minimap2 -ax map-hifi -t $T -R '@RG\tID:SRR18210286\tSM:SRR18210286\tPL:PACBIO' asm.fa hifi.fq.gz \
  | samtools sort -@ 2 -T "tmp/$SLURM_JOB_ID/hifi_sort" -o hifi.bam -
samtools index hifi.bam
EOF

bash -n align.sh
ALIGN_JOB_ID=$(sbatch --parsable align.sh)
printf '%s\n' "$ALIGN_JOB_ID" > align_job_id.txt
squeue -j "$ALIGN_JOB_ID"
```

> Check back at Part 3.

**Question 2.** The read groups record separate BGISEQ and HiFi sample identifiers as well as run and platform information. Why preserve these identifiers when both samples are labeled S288C?

<details>
<summary>Answer</summary>

The strain name is not a unique physical-sample identifier. The reads came from different BioSamples, so their provenance must remain distinct. Read groups preserve sample and run information inside the BAM and are required by tools such as GATK.

The `ID` identifies the read group; `SM` identifies the sample. `DNBSEQ` is the SAM platform label for BGISEQ/MGI instruments. Do not use one shared `SM` merely because two libraries have the same strain label.

</details>

**Checkpoint 1.** Everyone's alignment job is submitted and its ID saved.

---

## Part 2: Read A SAM File While It Runs

You do not need the full BAM to learn to read one. Align ten thousand pairs interactively:

```bash
seqkit head -n 10000 dnb_1.fq.gz -o peek_1.fq.gz
seqkit head -n 10000 dnb_2.fq.gz -o peek_2.fq.gz
bwa-mem2 mem -t 8 -R '@RG\tID:SRR17374239\tSM:SRR17374239\tLB:lib1\tPL:DNBSEQ' \
  asm.fa peek_1.fq.gz peek_2.fq.gz > peek.sam 2> logs/peek.log
```

### 2.1 The Header

```bash
samtools view -H peek.sam | cut -c1-80 | head -5
samtools view -H peek.sam | grep -E '^@(RG|PG)' | cut -c1-80
```

The `@SQ` lines carry your contig names and lengths. A BAM made against the published S288C reference would name chromosomes `chrI`, `chrII` and so on, and the two files could never be compared directly.

### 2.2 One Record, Field By Field

```bash
samtools view peek.sam | head -2 | cut -f1-9
```

**Question 3.** Pick the first record. What are its flag, contig, position, MAPQ and CIGAR, and where is its mate?

<details>
<summary>Answer</summary>

Most first records are a properly paired read: flag 99 or 83 (forward or reverse first-in-pair, properly paired), MAPQ 60, a CIGAR such as `100M` or `150M` (whatever the read length is), `=` in field 7 (mate on the same contig) and a mate position a few hundred bases away. `samtools flags` decodes any flag you are unsure of.

Reading one record completely is the skill; everything else in this lab is that, repeated with `awk`.

</details>

### 2.3 The Flags, Counted

```bash
samtools view peek.sam | cut -f2 | sort | uniq -c | sort -rn | head
samtools flags 99
samtools flags 2145
```

**Question 4.** Most reads carry flags 83, 99, 147 or 163. What do the rarer flags above 2048 mean here, and why does a fragmented assembly produce more of them than a complete one?

<details>
<summary>Answer</summary>

Flags with 2048 set are supplementary alignments: part of a read aligned in one place and the rest somewhere else. In a reference genome they usually mark structural variants.

Against your own assembly they also appear where a read crosses the end of one contig and continues into a sequence that sits at the start of another. Every contig break is a place the genome continues but the assembly does not, so more breaks mean more split reads.

</details>

### 2.4 Mapping Quality And Clipping

```bash
samtools view peek.sam | \
  awk '{ if ($5 == 0) z++; else if ($5 < 30) l++; else h++ } END {print "MAPQ 0:", z, "  1-29:", l, "  30+:", h}'
samtools view peek.sam | awk '$6 ~ /S/' | wc -l
```

**Question 5.** What can produce MAPQ 0 reads when short genomic reads are aligned to a yeast assembly? Does the S288C label prove that every ambiguous alignment is a repeat problem?

<details>
<summary>Answer</summary>

Repeated sequences are a major source: a short read can fit several Ty elements, rDNA copies or subtelomeric repeats similarly well. A read must reach discriminating sequence to be placed confidently.

Low-quality or low-complexity reads, assembly duplications and other ambiguous alignments can also lower MAPQ. A strain label does not prove ploidy, identical samples or perfect sequence quality. Use the alignment and coverage evidence to test an explanation.

</details>

**Checkpoint 2.** Everyone can read a SAM record and say what its flag, MAPQ and CIGAR mean.

---

## Part 3: QC The Real Alignment

```bash
squeue -j "$(cat align_job_id.txt)"
ls -lh sample.bam hifi.bam
```

### 3.1 The Headline Numbers

```bash
samtools flagstat -@ 4 sample.bam
samtools stats -@ 4 sample.bam | grep ^SN | cut -f2-3 | \
  grep -E 'reads mapped:|reads properly paired:|reads duplicated:|insert size average:|error rate:'
```

**Question 6.** What mapping rate did you observe, and what would you investigate if it were only 90%?

<details>
<summary>Answer</summary>

A high mapping rate is expected for related S288C genomic libraries, but there is no guaranteed percentage. The samples were not verified as one DNA preparation, and the assembly may omit some sequences.

Investigate adapter sequence, contamination, read quality, missing organellar or plasmid sequence, and differences between the cultures. A 90% mapping rate does not establish that 10% of the genome is missing: the denominator is reads, and genomic regions are not sampled equally. Inspect unmapped reads and their QC before assigning a cause.

</details>

### 3.2 Coverage, Contig By Contig

```bash
mosdepth -t 4 --by 2000 --fast-mode all sample.bam
mosdepth -t 4 --by 2000 --fast-mode --mapq 20 q20 sample.bam

# fraction of the genome at 10x or more, all reads and MAPQ 20+ only
awk '$1=="total" && $2==10 {print "all reads:", $3}' all.mosdepth.global.dist.txt
awk '$1=="total" && $2==10 {print "MAPQ 20+: ", $3}' q20.mosdepth.global.dist.txt

# mean depth per contig, all reads against MAPQ 20+
paste all.mosdepth.summary.txt q20.mosdepth.summary.txt | grep -v _region | \
  awk 'NR>1 {printf "%-14s %9d %9.1f %9.1f\n", $1, $2, $4, $10}' | sort -k3,3nr | head -20
```

**Question 7.** A few contigs sit at many times the typical depth. Name the likely candidates in a yeast genome, and say why each is over-represented.

<details>
<summary>Answer</summary>

- **The mitochondrial genome** (about 86 kb in S288C): many copies per cell.
- **The 2-micron plasmid** (about 6 kb), if hifiasm assembled it: dozens of copies per cell.
- **rDNA**: about 150 tandem copies of a 9 kb unit, which the assembly collapses into one or a few. Every copy's reads pile onto the few that were assembled.

The first two are real copy number; the third is a collapsed repeat. Depth alone does not tell them apart. Contig length, GC (the mitochondrial genome is very AT-rich) and Homework 3's GC-against-depth table do.

</details>

### 3.3 Where Short Reads Cannot Be Placed

```bash
zcat all.regions.bed.gz > all.bed
zcat q20.regions.bed.gz > q20.bed
paste all.bed q20.bed | awk '$4 > 10 && $8 < 0.5*$4' > mapq0_islands.bed
wc -l mapq0_islands.bed
head mapq0_islands.bed
```

Each line is a 2 kb window where more than half the depth comes from reads with MAPQ below 20.

**Question 8.** These windows are covered, but mostly by reads that cannot be trusted to be there. Why does a mean-depth summary hide them, and what does that mean for Thursday?

<details>
<summary>Answer</summary>

Mean depth counts every read wherever the aligner put it, including the ones it placed at random. These windows look adequately covered and are not: no variant caller should, or by default will, trust MAPQ 0 reads.

For Session 2 it means variant calls are blind in these windows, and any calls that do appear there are suspect. Report breadth at a stated depth and mapping quality, not mean depth alone.

</details>

**Checkpoint 3.** Everyone has a flagstat, a per-contig depth table and a list of MAPQ 0 islands.

---

## Part 4: The Same Genome, Seen By Long Reads

```bash
mosdepth -t 4 --by 2000 --fast-mode --mapq 20 hq20 hifi.bam
zcat hq20.regions.bed.gz > hq20.bed

# in each MAPQ 0 island: short-read depth, short MAPQ 20+ depth, HiFi MAPQ 20+ depth
paste all.bed q20.bed hq20.bed | awk '$4 > 10 && $8 < 0.5*$4 {printf "%s:%d-%d  %6.1f %6.1f %6.1f\n", $1, $2+1, $3, $4, $8, $12}' | head
```

**Question 9.** In most of these windows, HiFi reads with MAPQ 20 or more cover what short reads could not place. Why, and what does it mean for how this assembly was built?

<details>
<summary>Answer</summary>

A HiFi read is longer than a Ty element, so it reaches unique sequence on both sides and can only fit one copy. A 150 bp read inside the element fits every copy.

That is also why the assembly resolved these copies at all: hifiasm built it from the same long reads. The repeats that defeat short-read alignment are the same ones that defeat short-read assembly. The lecture's threshold argument, from Week 3, seen from the alignment side.

</details>

**Checkpoint 4.** `sample.bam` and its index exist. Session 2 calls variants on it.

---

## Takeaways

| Task | Command |
|---|---|
| Index a reference | `bwa-mem2 index asm.fa` and `samtools faidx asm.fa` |
| Align, mark duplicates, index | `bwa-mem2 mem -R … \| samtools fixmate -m \| samtools sort \| samtools markdup` |
| Long reads | `minimap2 -ax map-hifi` |
| Decode a flag | `samtools flags 99` |
| Headline QC | `samtools flagstat`, `samtools stats` |
| Coverage, all reads and trusted reads | `mosdepth --by 2000` with and without `--mapq 20` |

**Three things to carry forward:**

1. Choose reads by what the library measured. The Illumina run in this project is RNA-seq and cannot answer questions about genome coverage.
2. Mean depth includes reads placed at random. Report breadth at a stated depth and mapping quality.
3. Repeats often make short reads ambiguous; longer reads can reach unique flanks. Check sample provenance and alignment evidence before assigning a cause.
