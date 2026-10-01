# Week 2 Session 2: Trimming, Filtering, and Reading QC Across Samples

Applied Bioinformatics - Tufts University Department of Biology

On Tuesday you measured four platforms. Today you act on what you found: trim the short
reads, filter the long reads, aggregate everything into one report, and diagnose four
samples that are deliberately broken.

The last part is the point of the session. Running QC is easy. Deciding what a bad plot
means is the skill.

## Before You Start

```bash
srun --pty -p interactive -n 4 --mem=16g --time=0-03:00:00 bash
hostname

export COURSE=/cluster/tufts/bio_appbio
export MYWORK=$COURSE/$USER
export SUB=$MYWORK/week-02/session-01/sub

mkdir -p "$MYWORK/week-02/session-02"
cd "$MYWORK/week-02/session-02"

ls -lh "$SUB"
```

If `$SUB` is empty, you did not finish Tuesday's lab. Use the staged subsamples:

```bash
export SUB=$COURSE/shared/wk2/sub
ls -lh "$SUB"
```

```bash
module load miniforge
conda activate appbio-week1
conda install -y -c conda-forge -c bioconda fastp multiqc chopper
```

## Part 1: Trim the Short Reads

### 1.1 Run fastp on One Sample

```bash
mkdir -p trimmed reports

S=$(basename $(ls "$SUB"/*_1.fastq.gz | head -1) _1.fastq.gz)
echo "working on $S"

fastp \
  --in1 "$SUB/${S}_1.fastq.gz" \
  --in2 "$SUB/${S}_2.fastq.gz" \
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
for f in "$SUB"/*_1.fastq.gz; do
  S=$(basename "$f" _1.fastq.gz)
  [ -f "reports/${S}.fastp.json" ] && continue
  fastp \
    --in1 "$SUB/${S}_1.fastq.gz" --in2 "$SUB/${S}_2.fastq.gz" \
    --out1 "trimmed/${S}_1.trim.fq.gz" --out2 "trimmed/${S}_2.trim.fq.gz" \
    --detect_adapter_for_pe \
    --json "reports/${S}.fastp.json" --html "reports/${S}.fastp.html" \
    --thread 4 2>/dev/null
  echo "trimmed $S"
done
ls reports/
```

**Checkpoint 1.** Both short-read samples trimmed, with JSON and HTML reports.

## Part 2: Long Reads Are Filtered, Not Trimmed

### 2.1 Why the Difference

A short read has adapter at a predictable end and quality that decays toward the 3' end,
so trimming the ends is meaningful. A long read has errors distributed along its whole
length. Trimming the ends of a 20 kb read removes almost nothing useful and leaves the
problem untouched.

So for long reads you make a keep-or-discard decision per read, on length and mean
quality.

### 2.2 Filter with chopper

```bash
for f in "$SUB"/*.fastq.gz; do
  base=$(basename "$f" .fastq.gz)
  case "$base" in *_1|*_2) continue ;; esac
  zcat "$f" | chopper -q 10 -l 1000 --threads 4 2>/dev/null | gzip > "trimmed/${base}.filt.fq.gz"
  echo "filtered $base"
done

seqkit stats -a "$SUB"/*.fastq.gz trimmed/*.filt.fq.gz 2>/dev/null | grep -v '_[12]'
```

**Question 2.** Compare the before and after rows. What happened to read count, total
bases, and N50?

<details>
<summary>Answer</summary>

Read count falls, total bases falls by less, and N50 goes up or stays flat.

That pattern is the filter working as intended: it removes short and low-quality reads,
which are numerous but contribute few bases. You lose a lot of reads and little sequence,
and what remains is longer on average.

If N50 had fallen you would have been discarding your best reads, which would mean the
thresholds were wrong.

</details>

**Question 3.** Why `-q 10` rather than the `-q 20` you might use for short reads?

<details>
<summary>Answer</summary>

Q10 is one error in ten; Q20 is one in a hundred. Those would be alarming thresholds for
Illumina.

For nanopore they are not, because the whole platform operates at a lower per-base
accuracy and the value comes from length and from consensus across many reads. Filtering
at Q20 would discard most of a perfectly usable run.

Thresholds are platform-specific. Carrying a short-read instinct onto long-read data
throws away good data.

</details>

**Checkpoint 2.** Both long-read files filtered, and you can explain what the filter
removed.

## Part 3: Aggregate Everything

### 3.1 One Report for All Four Platforms

```bash
mkdir -p fastqc_trimmed
fastqc -t 4 -o fastqc_trimmed trimmed/*.trim.fq.gz

multiqc --force --title "Week 2 all platforms" -o multiqc \
  "$MYWORK/week-02/session-01/fastqc" fastqc_trimmed reports
```

Open `multiqc/multiqc_report.html` through OnDemand Files.

### 3.2 Start at the General Statistics Table

Click the column headers to sort: **M Seqs**, **% Dups**, **% GC**, **% Passed**.

**Question 4.** You are looking at raw and trimmed versions of two platforms in one table.
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

**Question 5.** Pick one thing visible here that you could not have seen in any single
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

zcat long.fastq.gz | chopper -q 10 -l 1000 | gzip > filt.fq.gz
```

Aggregate before you decide anything:

```bash
multiqc --force -o multiqc fastqc_dir reports_dir
```

Three things to carry forward:

1. Over-trimming is worse than under-trimming. Trim adapters; leave quality mostly alone.
2. Thresholds are platform-specific. Q20 is routine for Illumina and would destroy a
   nanopore run.
3. Run MultiQC even on one sample, and always before deciding a dataset is fine. The
   problems that matter most are only visible across samples.
