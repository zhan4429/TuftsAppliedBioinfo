# Week 2 Session 2: Trimming, Filtering, and Reading QC Across Samples

Applied Bioinformatics - Tufts University Department of Biology

On Tuesday you measured four platforms. Today you act on what you found: trim the short
reads, clean and filter the long reads, aggregate everything into one report, and diagnose
four samples that are deliberately broken.

The last part is the point of the session. Running QC is easy. Deciding what a bad plot
means is the skill.

## Before You Start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=4 \
  --mem=16G --time=03:00:00 --pty bash
hostname

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK=$COURSE/$USER
export READS=$MYWORK/week-01/session-02/fetch/results/fastq

mkdir -p "$MYWORK/week-02/session-02"
cd "$MYWORK/week-02/session-02"

ls -lh "$READS"
```

Use the original downloaded FASTQ files, as in Session 1. If `$READS` is empty,
finish the Week 1 download or use the staged FASTQ files:

```bash
export READS=$COURSE/shared/wk2/sub
ls -lh "$READS"
```

```bash
module load seqkit fastp fastplong fastqc multiqc
module list
fastplong --version
```

`fastplong --version` should report 0.4.1, the version Week 3's assembly pipeline uses. A
different version can give different read counts from the same command.

Repeat the six filename assignments from
[Session 1: Identify the Input Files](week-02-01-qc-across-platforms.html#identify-the-input-files),
using the files in your current `$READS` directory. You will use `ILLUMINA_R1`,
`ILLUMINA_R2`, `BGISEQ_R1`, `BGISEQ_R2`, `PB`, and `LONG` below.

## Part 1: Trim the Short Reads

### 1.1 Run fastp on One Sample

```bash
mkdir -p trimmed reports

S=illumina
echo "working on $S"

fastp \
  --in1 "$ILLUMINA_R1" \
  --in2 "$ILLUMINA_R2" \
  --out1 "trimmed/${S}_1.trim.fq.gz" \
  --out2 "trimmed/${S}_2.trim.fq.gz" \
  --detect_adapter_for_pe \
  --json "reports/${S}.fastp.json" \
  --html "reports/${S}.fastp.html" \
  --thread 4
```

`--detect_adapter_for_pe` is not on by default for paired-end data, and it is almost
always what you want. That is the kind of default that quietly costs people a project.

### 1.2 How Much Survived

```bash
python3 -c "
import json, glob
d = json.load(open(sorted(glob.glob('reports/*.fastp.json'))[0]))
b = d['summary']['before_filtering']['total_reads']
a = d['summary']['after_filtering']['total_reads']
print(f'before: {b:,}')
print(f'after : {a:,}')
print(f'kept  : {100*a/b:.1f}%')
"
```

**Question 1.** You keep 97%. Is that good? What retention would worry you, and in which
direction?

<details>
<summary>Answer</summary>

97% is healthy; 95 to 99 is typical.

Below about 90% needs explaining, usually short inserts producing heavy adapter
read-through, or genuinely poor data.

But the direction people get wrong is the other one. Over-trimming is worse than
under-trimming: it removes real sequence, biases you toward reads that were easy to
sequence, and shortens reads so they map worse. Aligners already soft-clip bad ends. Trim
adapters; leave quality mostly alone.

</details>

### 1.3 Trim the Other Short-Read Sample

```bash
S=bgiseq
fastp \
  --in1 "$BGISEQ_R1" --in2 "$BGISEQ_R2" \
  --out1 "trimmed/${S}_1.trim.fq.gz" --out2 "trimmed/${S}_2.trim.fq.gz" \
  --detect_adapter_for_pe \
  --json "reports/${S}.fastp.json" --html "reports/${S}.fastp.html" \
  --thread 4
ls reports/
```

**Checkpoint 1.** Both short-read samples trimmed, with JSON and HTML reports.

## Part 2: Long Reads Are Filtered, Not Quality-Trimmed

### 2.1 Why the Difference

A short read has adapter at a predictable end and quality that decays toward the 3' end,
so trimming the ends is meaningful. A long read has errors distributed along its whole
length. Trimming the ends of a 20 kb read removes almost nothing useful and leaves the
problem untouched.

So for long reads you do two things. Remove any adapter left at the read ends, as you did
for short reads. Then make a keep-or-discard decision per read, on length and mean quality.

`fastplong` does both in one pass. It comes from the authors of `fastp` and writes the
same kind of JSON and HTML reports. It is also the first step of Week 3's assembly
pipeline, nf-core/genomeassembler, which runs it on every long-read file. The settings you
choose today are the ones you will give the pipeline.

### 2.2 Start With the Defaults

Run `fastplong` on the Nanopore reads with no filtering options. Without `--out` it writes
only its reports, not a new FASTQ file.

```bash
mkdir -p long_reports

base=$(basename "$LONG" .fastq.gz)
fastplong \
  --in "$LONG" \
  --json "long_reports/${base}.default.fastplong.json" \
  --html "long_reports/${base}.default.fastplong.html" \
  --thread 4
```

You asked for no filtering, but `fastplong` still filtered. Its defaults come from `fastp`:
it trims adapters, discards a read if more than 40% of its bases are below Q15, and
discards reads shorter than 20 bp. Nothing in the command shows any of that.

### 2.3 Set the Thresholds Yourself

```bash
for f in "$LONG" "$PB"; do
  base=$(basename "$f" .fastq.gz)
  fastplong \
    --in "$f" \
    --out "trimmed/${base}.filt.fq.gz" \
    --mean_qual 10 \
    --length_required 1000 \
    --unqualified_percent_limit 100 \
    --json "long_reports/${base}.fastplong.json" \
    --html "long_reports/${base}.fastplong.html" \
    --thread 4
  echo "filtered $base"
done
```

- `--mean_qual 10` discards reads whose average quality is below Q10.
- `--length_required 1000` discards reads shorter than 1,000 bp.
- `--unqualified_percent_limit 100` turns off the per-base rule from 2.2, so average
  quality and length are the only filters.
- Adapter trimming stays on, as it is by default.

The options are written out in full so that anyone reading the command knows what it did.
Week 3 gives this same set to the assembly pipeline.

### 2.4 What Was Removed

Summarize the three reports in one table:

```bash
python3 - <<'EOF'
import json, glob, os

print(f"{'report':<34}{'reads in':>12}{'reads out':>12}{'kept':>8}{'low qual':>10}{'too short':>11}{'adapter':>9}")
for path in sorted(glob.glob("long_reports/*.fastplong.json")):
    d = json.load(open(path))
    before = d["summary"]["before_filtering"]["total_reads"]
    after = d["summary"]["after_filtering"]["total_reads"]
    failed = d["filtering_result"]
    adapter = d.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)
    name = os.path.basename(path).replace(".fastplong.json", "")
    print(f"{name:<34}{before:>12,}{after:>12,}{100*after/before:>7.1f}%"
          f"{failed['low_quality_reads']:>10,}{failed['too_short_reads']:>11,}{adapter:>9,}")
EOF
```

`low qual` and `too short` count the reads removed by each kind of rule. `adapter` counts
reads that had adapter trimmed from an end. Then compare the lengths before and after:

```bash
seqkit stats -a -j 4 "$LONG" "$PB" trimmed/*.filt.fq.gz
```

Open the `.html` files in `long_reports/` through OnDemand Files to see the length and
quality plots.

**Question 2.** Compare the before and after rows. What happened to read count, total
bases, and N50 in each file? Why did the length filter remove nothing?

<details>
<summary>Answer</summary>

For the Nanopore file, any reads removed appear under `low qual`. Compare the share of
reads removed with the share of bases removed. Low-quality Nanopore reads are often
shorter than average, so total bases usually falls by proportionally less than read count,
and N50 rises or stays flat. If N50 had fallen, the filter would have been removing your
longest reads, which would mean the thresholds were wrong.

The PacBio HiFi file barely changes. A HiFi read is a consensus built from several passes
over the same molecule, and a read is only called HiFi at about Q20 or better, so a Q10
threshold has almost nothing to remove. Running it still checks the read ends for adapter.

`too short` is zero, or nearly so, for both files because nothing was short to begin
with. Session 1 showed that the shortest Nanopore read is 2,000 bp and the shortest HiFi
read is 1,332 bp: the submitters removed short reads before depositing the data. A
threshold below the data's minimum does nothing. Check the data before assuming a filter
did its job.

</details>

**Question 3.** Why `--mean_qual 10` rather than the Q20 you might use for short reads?

<details>
<summary>Answer</summary>

Q10 is one error in ten; Q20 is one in a hundred. Those would be alarming thresholds for
Illumina.

For nanopore they are not, because the whole platform operates at a lower per-base
accuracy and the value comes from length and from consensus across many reads. Filtering
at Q20 would discard most of a perfectly usable run.

Thresholds are platform-specific. Carrying a short-read instinct onto long-read data
throws away good data.

One caution about the number itself. `fastplong` averages the quality scores directly.
Quality scores are on a log scale, so a few very poor bases raise the read's true error
rate far more than they lower this average. A read that averages Q10 here has an error
rate of at least one in ten, usually worse. Treat the threshold as a cut-off, not as a
precise error rate.

</details>

**Question 4.** Compare the two Nanopore rows: the default run and yours. Which kept more
reads, and which rule removed reads in each?

<details>
<summary>Answer</summary>

Both runs list their quality failures under `low qual`, but different rules put them
there. In the default run it is the per-base rule: more than 40% of bases below Q15. In
your run it is the average: below Q10.

The two rules measure different things. The per-base rule counts how many bases are poor;
a read's average can stay above Q10 even when a large share of its bases is below Q15. In
Session 1 only about 61% of Nanopore bases reached Q20, so many reads have a sizeable
share below Q15, and the per-base rule is likely the stricter of the two for these data.
Your table shows whether it was.

The larger point is that the default run filtered reads although you asked for no
filtering. Week 3's pipeline runs `fastplong` with its defaults unless you give it
options, which is why you will give it the ones you chose here.

</details>

**Checkpoint 2.** Both long-read files filtered with thresholds you set, and you can say
which rule removed reads in each Nanopore run.

## Part 3: Aggregate Everything

### 3.1 One Report for All Four Platforms

```bash
mkdir -p fastqc_trimmed
fastqc -t 4 -o fastqc_trimmed trimmed/*.trim.fq.gz

multiqc --force --title "Week 2 all platforms" -o multiqc \
  "$MYWORK/week-02/session-01/fastqc" fastqc_trimmed reports
```

Open `multiqc/multiqc_report.html` through OnDemand Files. MultiQC has no module for
`fastplong`, which is why the long-read reports went to `long_reports/` rather than
`reports/`.

### 3.2 Start at the General Statistics Table

Click the column headers to sort: **M Seqs**, **% Dups**, **% GC**, **% Passed**.

**Question 5.** You are looking at raw and trimmed versions of two platforms in one table.
What does sorting by %GC show you, and why is that reassuring?

<details>
<summary>Answer</summary>

All samples should sit at roughly the same %GC, because they are the same yeast genome.
Yeast is around 38%.

That is reassuring because it is an independent check that you are looking at what you
think you are. If one platform's %GC were markedly different, it would point at
contamination, a sample swap, or a coverage bias in the library prep rather than at
anything the sequencer did.

Agreement across platforms on a biological property is evidence; disagreement is a lead
worth following.

</details>

### 3.3 What One Report Cannot Show

Scroll to the overlaid per-base quality and adapter content plots.

**Question 6.** Pick one thing visible here that you could not have seen in any single
FastQC report.

<details>
<summary>Answer</summary>

Anything that is a comparison rather than a measurement: that one platform's quality
decays faster than another's, that trimming removed adapter from one sample and there was
none in the other to remove, that read counts differ by an order of magnitude.

The general point is that a sample's own report shows its numbers without context. A
sample with a tenth of the reads of its siblings looks fine alone and is obviously wrong
in a table.

</details>

**Checkpoint 3.** One MultiQC report covering raw and trimmed data from both short-read
platforms.

## Part 4: Four Samples That Are Not Fine

Your own data was reasonably well behaved. Public data usually is, because it was
deposited by people who had already fixed the obvious problems.

These four have not been fixed. **Each has exactly one thing wrong with it.**

```bash
mkdir -p mystery_qc
fastqc -t 4 -o mystery_qc "$COURSE"/shared/wk2/mystery/*.fastq.gz
multiqc --force --title "Mystery samples" -o multiqc_mystery mystery_qc
```

Open `multiqc_mystery/multiqc_report.html`.

For each sample write one sentence: **what is wrong, what is the evidence, what would you
do?**

```text
Sample A
Sample B
Sample C
Sample D
```

Useful extra checks:

```bash
cd mystery_qc
for z in *_fastqc.zip; do unzip -o -q "$z"; done
grep -H '^FAIL' */summary.txt | sed 's|/summary.txt||' | sort
cd ..
```

<details>
<summary>Answers — do not open until you have written yours</summary>

**Sample A — adapter contamination.** Adapter content rises sharply from around position
100 and overrepresented sequences include adapter. The insert was shorter than the read
length, so the sequencer read through into the adapter. Trimming with
`--detect_adapter_for_pe` fixes it. Not a reason to discard the sample.

**Sample B — low yield.** An order of magnitude fewer reads than the others. Invisible in
its own report and obvious in the MultiQC table. Under-loading, failed quantification, or
a pooling error. This sample cannot support the same analysis as the others; re-sequence
or exclude, and say which in your methods.

**Sample C — quality crash partway through the run.** Per-base quality is fine and then
collapses at a specific cycle, the same cycle in both mates. A run-level problem — a
reagent issue or a bubble in the flow cell. Hard-trim to before the crash and accept
shorter reads, or re-sequence. Check whether other samples on the same run show it, which
tells you whether it is the sample or the run.

**Sample D — contamination.** The per-sequence GC plot has two peaks rather than one, and
overall %GC does not match the others. Another organism in the library, or a sample swap.
Trimming will not touch this. Screen taxonomically before doing anything else.

**The general lesson.** A, C and D would each have been missed by looking only at
pass/fail flags. B was invisible in its own report entirely. QC is reading numbers across
samples, not collecting green ticks.

</details>

**Checkpoint 4.** Four diagnoses written down before you opened the answers.

## Takeaways

Trim short reads, filter long reads:

```bash
fastp --in1 R1.fq.gz --in2 R2.fq.gz --out1 t1.fq.gz --out2 t2.fq.gz \
      --detect_adapter_for_pe --json report.json

fastplong --in long.fastq.gz --out filt.fq.gz \
          --mean_qual 10 --length_required 1000 --unqualified_percent_limit 100 \
          --json long.fastplong.json
```

Aggregate before you decide anything:

```bash
multiqc --force -o multiqc fastqc_dir reports_dir
```

Three things to carry forward:

1. Over-trimming is worse than under-trimming. Trim adapters; leave quality mostly alone.
2. Thresholds are platform-specific, and so are defaults. Q20 is routine for Illumina and
   would destroy a nanopore run, and a default carried over from short reads can filter
   long reads without telling you.
3. Run MultiQC even on one sample, and always before deciding a dataset is fine. The
   problems that matter most are only visible across samples.
