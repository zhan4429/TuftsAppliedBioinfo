# Week 2 Session 1: QC Across Four Sequencing Platforms

Applied Bioinformatics - Tufts University Department of Biology

Last week you fetched four yeast sequencing datasets from Illumina, BGISEQ, PacBio HiFi
and Oxford Nanopore. In this lab you measure the downloaded reads and use QC tools
suited to short or long reads.

| Run | Platform | Library type |
| --- | --- | --- |
| SRR27956204 | Illumina NovaSeq 6000 | Paired-end RNA-seq |
| SRR17374239 | BGISEQ-500 | Paired-end whole-genome sequencing |
| SRR18210286 | PacBio HiFi | Single-end whole-genome sequencing |
| SRR17374240 | Oxford Nanopore | Single-end whole-genome sequencing |

**The Illumina run contains RNA-seq reads; the other three contain genomic reads.**
These are separate datasets with different BioSample records, rather than a controlled
experiment sequencing the same sample on four platforms. Differences in QC reports can
reflect platform, library type, sample preparation, and read processing. Use the
metadata when interpreting the plots; a difference between files does not establish a
platform effect.

## Before You Start

Get a compute node. This lab runs real tools and does not belong on a login node.

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=4 \
  --mem=16G --time=03:00:00 --pty bash
hostname
```

Set your paths:

```bash
export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK=$COURSE/$USER
export READS=$MYWORK/week-01/session-02/fetch/results/fastq

mkdir -p "$MYWORK/week-02/session-01"
cd "$MYWORK/week-02/session-01"

ls -lh "$READS"
```

You should see six files: two paired-end runs with `_1` and `_2`, and two single-file
long-read runs. Use the samplesheet in
`$MYWORK/week-01/session-02/fetch/results/samplesheet/samplesheet.csv` to match files to
runs; filenames may contain experiment identifiers rather than the run accessions above.

If that directory is empty, your download did not finish. Use the staged copies and sort
yours out afterwards:

```bash
export READS=$COURSE/shared/wk2/fallback_fastq
ls -lh "$READS"
```

Load the software:

```bash
module load seqkit fastqc nanoplot
module list
```

`module load` makes these programs available in your shell. SeqKit provides sequence
statistics, FastQC creates read-QC reports, and NanoPlot plots long-read data. The
NanoPlot command is written `NanoPlot` with capital letters, even though its module is
named `nanoplot`.

### Identify the Input Files

Use the fetchngs samplesheet to identify the two Illumina mate files, the two BGISEQ
mate files, and the single PacBio and Nanopore files. Replace each filename below with
the actual filename in your `$READS` directory. Keep the variable names on the left.

```bash
ILLUMINA_R1="$READS/REPLACE_WITH_ILLUMINA_R1.fastq.gz"
ILLUMINA_R2="$READS/REPLACE_WITH_ILLUMINA_R2.fastq.gz"
BGISEQ_R1="$READS/REPLACE_WITH_BGISEQ_R1.fastq.gz"
BGISEQ_R2="$READS/REPLACE_WITH_BGISEQ_R2.fastq.gz"
PB="$READS/REPLACE_WITH_PACBIO.fastq.gz"
LONG="$READS/REPLACE_WITH_NANOPORE.fastq.gz"

ls -lh "$ILLUMINA_R1" "$ILLUMINA_R2" "$BGISEQ_R1" "$BGISEQ_R2" "$PB" "$LONG"
```

All six paths should point to existing files before you continue. `$LONG` names the
Nanopore file and `$PB` names the PacBio file. These variables refer directly to the
downloaded reads; this lab does not create a separate subset.

## Part 1: Measure All Four Side by Side

### 1.1 One Command, Every Platform

Start with a numerical summary of every FASTQ file. SeqKit reads compressed
`.fastq.gz` files directly, so you do not need to decompress them first.

```bash
seqkit stats -a "$READS"/*.fastq.gz
```

| Part | Meaning |
| --- | --- |
| `seqkit stats` | Calculate sequence counts, lengths, and other statistics, with one output row per file. |
| `-a` | Show all available statistics, including length quartiles, N50, and FASTQ quality summaries. |
| `"$READS"/*.fastq.gz` | Use every filename ending in `.fastq.gz` in the input directory. The shell expands `*` to the matching filenames. |

`num_seqs` is the read count, `sum_len` is the total number of bases, and `avg_len`
is the average read length. The command prints a table in the terminal and leaves
the FASTQ files unchanged.

Fill this in from the output. One row per file.

```text
file                 num_seqs    sum_len    avg_len    N50     Q20(%)   Q30(%)
```

**Question 1.** Compare `avg_len` across the four platforms. Which two are nearly
constant, and which two vary? Why?

<details>
<summary>Answer</summary>

Illumina and BGISEQ reads often have a narrow length distribution because sequencing
uses a chosen number of cycles. Trimming and other processing can shorten some reads,
so check the actual distributions in these files.

PacBio and Nanopore reads typically have broader length distributions. DNA fragment
length, library preparation, size selection, instrument settings, and read processing
all contribute. The comparison illustrates the datasets you have, not a universal
length distribution for each platform.

</details>

### 1.2 Quality Is Not One Number Either

Look at the `Q20(%)` and `Q30(%)` columns.

**Question 2.** Rank the four datasets by their reported Q30. Does that ranking tell you
which platform is best?

<details>
<summary>Answer</summary>

Use the values from your own output rather than an expected ranking. Reported base
quality is one measure; read length, yield, and the biological question also matter.
These runs differ in library type and sample preparation, so the ranking cannot isolate
the effect of platform.

Week 3 makes the choice concrete: you will use the PacBio HiFi genomic reads for genome
assembly. The Illumina RNA-seq run measures expressed transcripts and serves a different
purpose.

</details>

### 1.3 Read Length as a Distribution

An average hides the shape. Look at the real distribution for one long-read file:

`seqkit fx2tab` converts sequence records into a tab-separated table. Here we request
only read names and lengths, then summarize those lengths with the rest of the pipeline.

```bash
echo "$LONG"

seqkit fx2tab -nl "$LONG" | awk '{print $NF}' | sort -n | \
  awk '{a[NR]=$1} END {
    printf "min    %8d\n", a[1]
    printf "Q1     %8d\n", a[int(NR*0.25)]
    printf "median %8d\n", a[int(NR*0.50)]
    printf "Q3     %8d\n", a[int(NR*0.75)]
    printf "max    %8d\n", a[NR] }'
```

| Part | Meaning |
| --- | --- |
| `seqkit fx2tab` | Write one table row for each read. |
| `-n` | Keep read names without printing their sequences or quality strings. |
| `-l` | Add each read's length in bases. |
| `-nl` | A shorter way to write `-n -l`; both options apply. |
| `"$LONG"` | Read the Nanopore FASTQ file selected in the setup. |

`echo` displays the chosen filename. Each pipe (`|`) passes output to the next command.
The first `awk` extracts the final field, the read length; `sort -n` sorts the lengths
numerically. The last `awk` reports the minimum, quartiles, median, and maximum.
The backslash (`\`) continues the command on the next line.

Use the same SeqKit options for one Illumina mate file. This time, count how often
each read length occurs:

```bash
seqkit fx2tab -nl "$ILLUMINA_R1" | awk '{print $NF}' | sort -n | uniq -c | head
```

`"$ILLUMINA_R1"` selects the first Illumina mate file. `-nl` again requests names and
lengths. After extracting and sorting the lengths, `uniq -c` counts equal adjacent
values, and `head` shows the first ten length groups. The two output columns are the
number of reads and their length in bases.

**Question 3.** The short-read file gives one length repeated hundreds of thousands of
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

**Checkpoint 1.** You have a filled-in stats table and can describe differences between
the datasets while considering both platform and library type.

## Part 2: FastQC, and Where It Applies

### 2.1 Run It on the Short Reads Only

FastQC checks each input file for patterns such as declining base quality, adapter
content, and unusual sequence composition. It creates a separate report for each
mate file and does not trim or change the reads.

```bash
mkdir -p fastqc
fastqc -t 4 --outdir fastqc \
  "$ILLUMINA_R1" "$ILLUMINA_R2" "$BGISEQ_R1" "$BGISEQ_R2"
ls fastqc/
```

| Part | Meaning |
| --- | --- |
| `mkdir -p fastqc` | Create the output directory if needed; FastQC expects it to exist. |
| `fastqc` | Run the quality-control checks on the listed input files. |
| `-t 4` | Set the processing thread limit to four, matching the four CPUs requested above. |
| `--outdir fastqc` | Save the reports in the `fastqc/` directory. `-o fastqc` means the same thing. |
| The four filename variables | Analyze the two Illumina mate files and the two BGISEQ mate files. Each file gets its own report. |
| `ls fastqc/` | List the reports that were created. |

For each input, expect an HTML report for viewing and a ZIP archive containing report
data. The backslash after `--outdir fastqc` continues the command onto the next line.

View a report through OnDemand: **Files > Home Directory**, navigate to your
`week-02/session-01/fastqc/` directory, and click an `.html` file.

### 2.2 Read It Properly

Work through these panels for one Illumina sample and write down what you see.

```text
Basic Statistics          total sequences, length, %GC
Per base sequence quality where does it start dropping?
Per base sequence content skewed at the start, or all the way through?
Per sequence GC content   one peak, or two?
Adapter content           does it rise, and from which position?
Overrepresented sequences anything named?
```

**Question 4.** Compare the Illumina and BGISEQ reports. Both are short-read platforms.
Name one panel where they differ. Why should you consider library type as well as
platform when interpreting that difference?

<details>
<summary>Answer</summary>

Any accurate observation supported by the reports earns credit. The Illumina run is
RNA-seq, while the BGISEQ run is whole-genome sequencing. Expressed transcripts are
present at different abundances, which can affect duplication and sequence composition
in RNA-seq. Higher duplication alone does not establish a library failure.

Quality profiles can also reflect sequencing chemistry and read processing; adapter
content depends on library preparation. Describe what differs, then explain plausible
causes without assigning every difference to the platform. FastQC warnings need to be
interpreted in the context of the library type.

</details>

### 2.3 Try It Where It Does Not Belong

Run the same FastQC checks on the Nanopore file to examine how its report behaves
with long reads:

```bash
fastqc -t 4 -o fastqc "$LONG"
```

Here `-t 4` keeps the same thread limit, `-o fastqc` saves the report alongside the
short-read reports, and `"$LONG"` selects the Nanopore input. This is a demonstration
of tool interpretation; it does not convert or filter the reads.

Open that report.

**Question 5.** Look at Basic Statistics and Per base sequence quality. Why is this
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

**Checkpoint 2.** You have FastQC reports for both short-read platforms and have seen one
fail on long reads.

## Part 3: The Right Tool for Long Reads

### 3.1 NanoPlot

NanoPlot generates **plots and an HTML report**, as well as a text summary. It can
summarize both Nanopore and PacBio HiFi FASTQ reads; its name does not limit it to
Nanopore data.

```bash
mkdir -p nanoplot
NanoPlot --fastq "$LONG" -o nanoplot -t 4 --N50 --loglength
ls nanoplot/
cat nanoplot/NanoStats.txt
```

| Option | Meaning |
| --- | --- |
| `--fastq "$LONG"` | Read the selected Nanopore FASTQ file, including its stored quality scores. Compressed `.gz` input is supported. |
| `-o nanoplot` | Write the report, statistics, and plots into `nanoplot/`. |
| `-t 4` | Allow up to four processing threads, matching the compute allocation. |
| `--N50` | Add an N50 marker to the read-length histogram. |
| `--loglength` | Also produce plots with read length on a logarithmic scale. This changes the visualization, not the input reads. |

`mkdir -p` prepares the output folder, `ls` lists its contents, and `cat` displays the
text summary. These commands inspect the output without changing the reads.

The main outputs are:

| Output | What to use it for |
| --- | --- |
| `nanoplot/NanoPlot-report.html` | Open in a browser to explore the plots and summary statistics. |
| `nanoplot/NanoStats.txt` | Read the numerical summary in the terminal. |
| Other plot files in `nanoplot/` | Open individual plots; filenames and available formats depend on the installed version. |

**Open the report.** In OnDemand, choose **Files** and navigate to
`$MYWORK/week-02/session-01/nanoplot/` in your course directory. Open
`NanoPlot-report.html` in a browser. If OnDemand downloads HTML instead of displaying
it, download the output folder and open the report locally, keeping its files together.

Look for these plots:

- **Read-length distribution:** shows how many reads fall in each length range. `--N50` adds an N50 marker; at least half of the sequenced bases are in reads at least this long.
- **Read-quality distribution:** shows the spread of average quality scores across reads, when available.
- **Read length versus average quality:** shows whether longer and shorter reads have different reported qualities. The report may show points or a density view.

`--loglength` adds plots with a logarithmic length axis, making broad length ranges
easier to compare. Some quality plots may be skipped when the scores have no variation.
These plots summarize the scores stored in the FASTQ file; they do not independently
measure sequencing accuracy.

### 3.2 Five Numbers That Matter

From `NanoStats.txt`, record:

```text
Mean read length
Median read length
Read length N50
Mean read quality
Number of reads
```

**Question 6.** NanoPlot reports both mean and median read length, and they differ
substantially. Which is larger, and what does that tell you about the run?

<details>
<summary>Answer</summary>

The mean is larger than the median, because the distribution has a long right tail — a
small number of very long reads pull the average up while most reads are shorter.

N50 is larger still, because it weights by bases rather than by reads. The gap between
median and N50 is a quick read on how much of your sequence sits in long molecules, which
is what matters for assembly.

</details>

### 3.3 Compare the Two Long-Read Platforms

Run NanoPlot on the PacBio file as well. Use a separate output directory so its report
is easy to distinguish from the Nanopore report:

```bash
mkdir -p nanoplot_pb
NanoPlot --fastq "$PB" -o nanoplot_pb -t 4 --N50 --loglength
grep -E 'Mean read length|Median read|N50|Mean read quality' nanoplot/NanoStats.txt nanoplot_pb/NanoStats.txt
```

This command uses the same NanoPlot options: `--fastq "$PB"` selects the PacBio
file, `-o nanoplot_pb` chooses its output folder, `-t 4` sets the thread limit,
`--N50` marks N50, and `--loglength` adds logarithmic length plots. Both runs therefore
use the same plot settings. The final `grep -E` command selects matching summary lines
from the two `NanoStats.txt` files; it does not run NanoPlot again.

Open **`nanoplot_pb/NanoPlot-report.html`** for the PacBio report, using the same
OnDemand or download steps as above. Compare the read-length and quality plots with
those in `nanoplot/NanoPlot-report.html`, as well as the numbers in `NanoStats.txt`.

**Question 7.** Both are long-read platforms. Which has higher read quality, and which
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

**Checkpoint 3.** You have opened both NanoPlot HTML reports, compared their plots and
summary statistics, and can describe the trade-off between the two datasets.

## Takeaways

Measure every dataset the same way first. `stats` prints one summary row per file,
`-a` includes the additional metrics, and the filename pattern selects all input FASTQs:

```bash
seqkit stats -a "$READS"/*.fastq.gz
```

Use FastQC for short reads. `-t 4` sets its thread limit, `-o fastqc` selects the
report folder, and the remaining arguments list the four short-read files:

```bash
fastqc -t 4 -o fastqc "$ILLUMINA_R1" "$ILLUMINA_R2" "$BGISEQ_R1" "$BGISEQ_R2"
```

Use NanoPlot for long reads. `--fastq` names the input, `-o` names the output folder,
`-t 4` sets the thread limit, and `--N50` adds the N50 marker:

```bash
NanoPlot --fastq "$LONG" -o nanoplot -t 4 --N50
```

Three things to carry forward:

1. Interpret QC with the metadata. The Illumina reads are RNA-seq; the BGISEQ, PacBio and
   Nanopore reads are genomic. Platform, library preparation and read processing all matter.
2. A tool that produces output is not a tool that applies. FastQC runs happily on nanopore
   data and tells you almost nothing true.
3. No platform is best. Accuracy, length, throughput and turnaround trade against each
   other, and which one you need is a property of your question.

## Dataset Records and QC Reference

- [Illumina RNA-seq: SRR27956204](https://www.ebi.ac.uk/ena/browser/view/SRR27956204)
- [BGISEQ genomic reads: SRR17374239](https://www.ebi.ac.uk/ena/browser/view/SRR17374239)
- [PacBio genomic reads: SRR18210286](https://www.ebi.ac.uk/ena/browser/view/SRR18210286)
- [Nanopore genomic reads: SRR17374240](https://www.ebi.ac.uk/ena/browser/view/SRR17374240)
- [FastQC: duplication and RNA-seq libraries](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)
- [NanoPlot documentation and output descriptions](https://github.com/wdecoster/NanoPlot)
- [SeqKit command and parameter reference](https://bioinf.shenwei.me/seqkit/usage/)
- [FastQC installation and command-line usage](https://github.com/s-andrews/FastQC/blob/master/INSTALL.md)
