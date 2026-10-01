# Week 2 Session 1: QC Across Four Sequencing Platforms

Applied Bioinformatics - Tufts University Department of Biology

Last week you fetched the same yeast genome sequenced four different ways: Illumina,
BGISEQ, PacBio HiFi and Oxford Nanopore. In this lab you subsample all four, measure
them, and run the right QC tool for each.

Because the sample is identical across all four, every difference you see is the
platform. That is unusual and it is the point of the exercise.

## Before You Start

Get a compute node. This lab runs real tools and does not belong on a login node.

```bash
srun --pty -p interactive -n 4 --mem=16g --time=0-03:00:00 bash
hostname
```

Set your paths:

```bash
export COURSE=/cluster/tufts/bio_appbio
export MYWORK=$COURSE/$USER
export READS=$MYWORK/week-01/session-02/fetch/results/fastq
export GENOME=12100000

mkdir -p "$MYWORK/week-02/session-01"
cd "$MYWORK/week-02/session-01"

ls -lh "$READS"
```

You should see six files: two paired-end runs with `_1` and `_2`, and two single-file
long-read runs.

If that directory is empty, your download did not finish. Use the staged copies and sort
yours out afterwards:

```bash
export READS=$COURSE/shared/wk2/fallback_fastq
ls -lh "$READS"
```

Load the software:

```bash
module load miniforge
conda activate appbio-week1
conda install -y -c conda-forge -c bioconda fastqc nanoplot
```

## Part 1: Subsample Before You Measure

### 1.1 Why This Comes First

Your Illumina run has about 45 million reads. FastQC on a file that size takes longer
than this lab slot, and it would tell you nothing that 400,000 reads will not.

Remind yourself of the scale:

```bash
cd "$MYWORK/week-01/session-02/fetch"
awk -F'\t' 'NR>1 && !seen[$2]++ {printf "%-13s %-16s %6.0fx coverage\n", $1, $2, $5/12100000}' runs.tsv
cd "$MYWORK/week-02/session-01"
```

### 1.2 Take a Slice of Each

```bash
mkdir -p sub

# short-read runs: 400,000 reads from each mate, order preserved so pairs stay matched
for f in "$READS"/*_1.fastq.gz; do
  base=$(basename "$f" _1.fastq.gz)
  seqkit head -n 400000 "$READS/${base}_1.fastq.gz" -o "sub/${base}_1.fastq.gz"
  seqkit head -n 400000 "$READS/${base}_2.fastq.gz" -o "sub/${base}_2.fastq.gz"
  echo "subsampled $base"
done

# long-read runs: fewer reads, because each one is far longer
for f in "$READS"/*.fastq.gz; do
  base=$(basename "$f" .fastq.gz)
  case "$base" in *_1|*_2) continue ;; esac
  seqkit head -n 20000 "$f" -o "sub/${base}.fastq.gz"
  echo "subsampled $base"
done

ls -lh sub/
```

**Question 1.** `seqkit head` takes the first N reads rather than a random sample. What
bias does that introduce, and why is it acceptable here?

<details>
<summary>Answer</summary>

Reads are stored roughly in the order they came off the instrument, so the first N come
disproportionately from one region of the flow cell. If that region ran badly, your slice
looks worse than the run really is.

It is acceptable here because you are comparing platforms, not certifying a dataset, and
the bias applies equally to all four. For a real QC report you would use `seqkit sample`
to draw randomly, which is slower because it has to read the whole file.

The habit worth keeping: know which shortcut you took, and say so.

</details>

**Checkpoint 1.** `ls sub/` shows six files, and they are far smaller than the originals.

## Part 2: Measure All Four Side by Side

### 2.1 One Command, Every Platform

```bash
seqkit stats -a sub/*.fastq.gz
```

Fill this in from the output. One row per file.

```text
file                 num_seqs    sum_len    avg_len    N50     Q20(%)   Q30(%)
```

**Question 2.** Compare `avg_len` across the four platforms. Which two are nearly
constant, and which two vary? Why?

<details>
<summary>Answer</summary>

Illumina and BGISEQ are nearly constant, because read length is set by the number of
sequencing cycles the operator chose. Every cluster is read the same fixed number of
times.

PacBio and Nanopore vary enormously, because read length is set by the length of the DNA
molecule that went through the pore or polymerase. That is determined by your extraction,
not by the instrument — which is why high molecular weight extraction is the part of a
long-read project that actually decides the result.

</details>

### 2.2 Quality Is Not One Number Either

Look at the `Q20(%)` and `Q30(%)` columns.

**Question 3.** Rank the four platforms by Q30. Does that ranking tell you which platform
is best?

<details>
<summary>Answer</summary>

Expect the short-read platforms highest, PacBio HiFi close behind, and Nanopore lowest.

It does not tell you which is best, because accuracy is only one axis. The nanopore run
has far lower Q30 and reads that are a hundred times longer, and length is what resolves
repeats. A platform is only "best" relative to a question.

Week 3 makes this concrete: you will assemble the HiFi reads, because they are the only
ones that are both long and accurate.

</details>

### 2.3 Read Length as a Distribution

An average hides the shape. Look at the real distribution for one long-read file:

```bash
LONG=$(ls sub/*.fastq.gz | grep -v '_[12].fastq.gz' | head -1)
echo "$LONG"

seqkit fx2tab -nl "$LONG" | awk '{print $NF}' | sort -n | \
  awk '{a[NR]=$1} END {
    printf "min    %8d\n", a[1]
    printf "Q1     %8d\n", a[int(NR*0.25)]
    printf "median %8d\n", a[int(NR*0.50)]
    printf "Q3     %8d\n", a[int(NR*0.75)]
    printf "max    %8d\n", a[NR] }'
```

Now the same for a short-read file:

```bash
seqkit fx2tab -nl sub/*_1.fastq.gz | awk '{print $NF}' | sort -n | uniq -c | head
```

**Question 4.** The short-read file gives one length repeated hundreds of thousands of
times. Why does that make the mean a useful summary there and a misleading one for long
reads?

<details>
<summary>Answer</summary>

When every value is the same, the mean describes the data perfectly.

For long reads the distribution spans orders of magnitude and is heavily skewed, so the
mean sits far from the typical read and far from where the bases actually are. A run with
a mean of 4 kb and an N50 of 22 kb is a good run with a lot of short junk in it. A mean of
4 kb and an N50 of 6 kb has no long reads at all. The mean cannot tell those apart.

This is why long-read reports lead with N50, and it comes back in Week 3 when contiguity
depends on whether reads span repeats.

</details>

**Checkpoint 2.** You have a filled-in stats table and can name one thing each platform
does that the others cannot.

## Part 3: FastQC, and Where It Applies

### 3.1 Run It on the Short Reads Only

```bash
mkdir -p fastqc
fastqc -t 4 -o fastqc sub/*_1.fastq.gz sub/*_2.fastq.gz
ls fastqc/
```

View a report through OnDemand: **Files > Home Directory**, navigate to your
`week-02/session-01/fastqc/` directory, and click an `.html` file.

### 3.2 Read It Properly

Work through these panels for one Illumina sample and write down what you see.

```text
Basic Statistics          total sequences, length, %GC
Per base sequence quality where does it start dropping?
Per base sequence content skewed at the start, or all the way through?
Per sequence GC content   one peak, or two?
Adapter content           does it rise, and from which position?
Overrepresented sequences anything named?
```

**Question 5.** Compare the Illumina and BGISEQ reports. Both are short-read platforms.
Name one panel where they differ.

<details>
<summary>Answer</summary>

Any accurate observation earns this. Common differences: the per-base quality profile has
a different shape, because the chemistry and error model differ; adapter content differs
because the library kits use different adapter sequences; %GC may differ slightly from
coverage bias.

The general point is that "short read" is not one thing. FastQC's thresholds were set
with Illumina in mind, so another short-read platform can trip flags that mean nothing.

</details>

### 3.3 Try It Where It Does Not Belong

```bash
fastqc -t 4 -o fastqc "$LONG"
```

Open that report.

**Question 6.** Look at Basic Statistics and Per base sequence quality. Why is this
report close to useless?

<details>
<summary>Answer</summary>

Sequence length is reported as a range spanning orders of magnitude, so "per base"
position has no consistent meaning — position 5,000 exists in a handful of reads and not
in most of them. The per-base quality plot bins positions that are not comparable, and the
tail is computed from a shrinking and unrepresentative set of reads.

FastQC also assumes a fixed length for several modules, and silently produces something
anyway rather than refusing.

A tool producing output is not the same as a tool being appropriate. That judgement is
yours.

</details>

**Checkpoint 3.** You have FastQC reports for both short-read platforms and have seen one
fail on long reads.

## Part 4: The Right Tool for Long Reads

### 4.1 NanoPlot

```bash
mkdir -p nanoplot
NanoPlot --fastq "$LONG" -o nanoplot -t 4 --N50 --loglength 2>/dev/null
cat nanoplot/NanoStats.txt
```

### 4.2 Four Numbers That Matter

From `NanoStats.txt`, record:

```text
Mean read length
Median read length
Read length N50
Mean read quality
Number of reads
```

**Question 7.** NanoPlot reports both mean and median read length, and they differ
substantially. Which is larger, and what does that tell you about the run?

<details>
<summary>Answer</summary>

The mean is larger than the median, because the distribution has a long right tail — a
small number of very long reads pull the average up while most reads are shorter.

N50 is larger still, because it weights by bases rather than by reads. The gap between
median and N50 is a quick read on how much of your sequence sits in long molecules, which
is what matters for assembly.

</details>

### 4.3 Compare the Two Long-Read Platforms

Run NanoPlot on the PacBio file as well, then compare:

```bash
PB=$(ls sub/*.fastq.gz | grep -v '_[12].fastq.gz' | sed -n 2p)
mkdir -p nanoplot_pb
NanoPlot --fastq "$PB" -o nanoplot_pb -t 4 --N50 2>/dev/null
grep -E 'Mean read length|Median read|N50|Mean read quality' nanoplot/NanoStats.txt nanoplot_pb/NanoStats.txt
```

**Question 8.** Both are long-read platforms. Which has higher read quality, and which
has longer reads? What does that trade-off mean for choosing between them?

<details>
<summary>Answer</summary>

PacBio HiFi has substantially higher quality; nanopore typically has longer reads.

HiFi earns its accuracy by reading the same molecule repeatedly and taking a consensus,
which costs length. Nanopore reads the molecule once, so length is limited only by your
DNA.

Choose HiFi when you need per-base accuracy — assembly you intend to annotate, variant
calling. Choose nanopore when you need to span something very long, or need the answer
today. Many projects use both, which is why this study has both.

</details>

**Checkpoint 4.** You have NanoPlot stats for both long-read platforms and can state the
trade-off between them.

## Takeaways

Measure every platform the same way first:

```bash
seqkit stats -a sub/*.fastq.gz
```

Subsample before QC, and know which shortcut you took:

```bash
seqkit head -n 400000 in.fastq.gz -o out.fastq.gz
```

Use FastQC for short reads and NanoPlot for long ones:

```bash
fastqc -t 4 -o fastqc sub/*_1.fastq.gz
NanoPlot --fastq long.fastq.gz -o nanoplot --N50
```

Three things to carry forward:

1. Read length is fixed by the operator on short-read platforms and by your DNA extraction
   on long-read platforms. That single difference explains most of what you measured today.
2. A tool that produces output is not a tool that applies. FastQC runs happily on nanopore
   data and tells you almost nothing true.
3. No platform is best. Accuracy, length, throughput and turnaround trade against each
   other, and which one you need is a property of your question.
