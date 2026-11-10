# Lab 10A — From Amplicon Reads to an ASV Table

**Applied Bioinformatics · Week 10, Session 1**
**Runs alongside the lecture**

> Twenty V4 samples from a two-group study, plus two extraction blanks and a mock community. You will denoise them, name what you find, and then ask what the numbers on your barplot actually are.

> **Work through Part 1, then submit the denoising job as soon as you have chosen your truncation lengths.** It runs while you answer the questions.

---

## Before you start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname                      # must NOT say login-p0x

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED=$COURSE/shared/wk10
export AMP=$SHARED/amplicon

mkdir -p $MYWORK/week-10/session-01/logs && cd $MYWORK/week-10/session-01
export TMPDIR=$SCRATCH/tmp && mkdir -p $TMPDIR
```

Load QIIME 2:

```bash
module load qiime2
qiime info | head -4
```

The module also provides `biom`, which you will use at the end. Every QIIME 2 output is either a `.qza` (data) or a `.qzv` (a visualization). To open a `.qzv`, download it through OnDemand (**Files**) and drag it onto [view.qiime2.org](https://view.qiime2.org): nothing is installed, and the file never leaves your computer.

```bash
head -4 $AMP/metadata.tsv
```

Note the `sample_type` column: two samples are blanks and one is the mock community. You will come back to all three.

---

## Part 1 — Import, remove primers, and look before you denoise

### 1.1 Import

```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $AMP/manifest.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path demux.qza
```

### 1.2 Remove the primers

Primer sequence is not biology. Left in place, it distorts the error model that DADA2 learns.

```bash
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --p-discard-untrimmed --p-cores 8 \
  --o-trimmed-sequences trimmed.qza
```

These are the Earth Microbiome Project 515F and 806R primers for V4.

### 1.3 Look at the reads

```bash
qiime demux summarize --i-data trimmed.qza --o-visualization trimmed.qzv
qiime tools export --input-path trimmed.qzv --output-path trimmed_summary

# reads per sample, fewest first
awk -F'\t' 'NR>1 {print $2"\t"$1}' trimmed_summary/per-sample-fastq-counts.tsv | sort -n | head -6
```

Now open `trimmed.qzv` and look at the **Interactive Quality Plot**. It is FastQC's per-base quality panel from Week 2, with one difference: the numbers you read off it go straight into the next command.

**Question 1.** Which samples have the fewest reads? Is that a problem?

::: Answer
The two blanks should be at or near the bottom, with far fewer reads than any real sample. That is the expected and reassuring result: little template went in, so little came out.

A blank with as many reads as a real sample would be a finding, and not a good one. It would mean contamination on the scale of a real sample, or a sample swap. A real sample near the blanks' level is a failed library and probably should not be analysed.
:::

**Question 2.** Pick a forward truncation length `F` and a reverse length `R` from the quality plot. Then check them: what do `F` and `R` have to satisfy for the pairs to merge?

::: Answer
The V4 region between these primers is about 253 bp, and DADA2 needs at least 12 bp of overlap to merge a pair. So:

`F + R − 253 ≥ 12`, plus a margin of 20 or so for length variation, so `F + R` should be at least about 285.

Within that limit, truncate where the median quality starts to fall. With 2 × 250 reads (about 230 bases each once primers are removed), something like `F = 220, R = 180` merges comfortably. Truncating both reads hard because the plot "looks bad" is the classic mistake: the reads stop overlapping, and most of them are lost at merging.

`F` and `R` cannot exceed the trimmed read lengths, or DADA2 discards every read.
:::

> **CHECKPOINT 1** — You have `trimmed.qza`, and two truncation lengths written down that you can defend.

---

## Part 2 — Denoise, as a job

```bash
cat > denoise.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=dada2
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
: "${TRUNC_F:?set TRUNC_F before submitting}"
: "${TRUNC_R:?set TRUNC_R before submitting}"

module load qiime2
export COURSE=/cluster/tufts/yzhang85/bioinfo_course/appliedBio
export MYWORK="$COURSE/$USER"
export TMPDIR=$SCRATCH/tmp && mkdir -p $TMPDIR
cd $MYWORK/week-10/session-01

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs trimmed.qza \
    --p-trunc-len-f $TRUNC_F --p-trunc-len-r $TRUNC_R \
    --p-n-threads ${SLURM_CPUS_PER_TASK} \
    --o-table table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats stats.qza
EOF
```

Fill in your two numbers, then submit:

```bash
export TRUNC_F=      # your forward truncation length
export TRUNC_R=      # your reverse truncation length

bash -n denoise.sh
DADA2_JOB=$(sbatch --parsable --export=ALL,TRUNC_F=$TRUNC_F,TRUNC_R=$TRUNC_R denoise.sh)
printf '%s\n' "$DADA2_JOB" > dada2_job_id.txt
squeue -j "$DADA2_JOB"
```

The script refuses to run if either number is empty, so a forgotten value fails immediately instead of after the job has queued.

**Question 3.** Your neighbour read the same quality plot and chose different numbers. Whose reads will merge better, and how could you tell without rerunning?

::: Answer
You can tell from the arithmetic in Question 2. If either pair of numbers leaves less than about 30 bp of overlap, that pair will lose reads at merging. If both leave plenty of overlap, the one that cut more of the low-quality tail will usually keep slightly more reads through denoising.

The general point: this choice can be checked before you run anything, and should be.
:::

**Question 4.** Why must samples from different sequencing runs be denoised separately?

::: Answer
DADA2 learns an error model from the data, and every run has its own error profile. Pooling runs would fit one model to two different error processes, which fits neither well. Denoise each run separately with the same settings, then merge the feature tables.
:::

> **CHECKPOINT 2** — Your denoising job is queued, and its ID is saved.

---

## Part 3 — Where the reads went

```bash
squeue -j "$(cat dada2_job_id.txt)"
ls -lh table.qza rep-seqs.qza stats.qza
```

If the job failed, read `logs/dada2_*.err` and fix it. If time is short, use the staged results: `cp $SHARED/fallback/{table,rep-seqs,stats}.qza .`

```bash
qiime tools export --input-path stats.qza --output-path stats
awk -F'\t' 'NR>2 {printf "%-12s %7d  merged %5.1f%%  kept %5.1f%%%s\n", $1, $2, $7, $9, ($9 < 50 ? "   <-- check" : "")}' stats/stats.tsv
```

**Question 5.** Look at the merged column first. What fraction of input reads survived to the end, and if any sample is flagged, what went wrong?

::: Answer
A healthy run keeps roughly 60–80% of input reads through to non-chimeric. A large drop between the denoised and merged columns means the truncation left too little overlap: change `TRUNC_F` and `TRUNC_R` and resubmit, which is the only fix.

The blanks may be flagged simply because they have very few reads. That is not a truncation problem.

Note what did not happen: nothing errored. A failed merge produces a perfectly good-looking ASV table with fewer reads in it, and this table is the only evidence.
:::

> **CHECKPOINT 3** — Everyone has a feature table and has checked where the reads went.

---

## Part 4 — Name the ASVs

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier $SHARED/silva-v4-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza \
  --m-metadata-file $AMP/metadata.tsv --o-visualization barplot.qzv
```

This uses one job on purpose: each parallel job loads its own copy of the classifier, and several copies can exhaust a node's memory.

How far down did the names go?

```bash
qiime tools export --input-path taxonomy.qza --output-path taxonomy
awk -F'\t' 'NR>1 {n++; if ($2 ~ /g__[^;]+/) g++; if ($2 ~ /s__[^;]+/) sp++}
  END {printf "%d ASVs: %.0f%% named to genus, %.0f%% to species\n", n, 100*g/n, 100*sp/n}' taxonomy/taxonomy.tsv
```

**Question 6.** Why is the species figure what it is, and why is "named to genus" still generous?

::: Answer
The classifier is trained on SILVA's V4 region without species labels, so the species figure is zero by design. Many species share an identical V4 sequence, and SILVA's own standard classifiers now leave species out. Genus is the working limit for 253 bp.

"Named to genus" is generous because some genus labels are placeholders, such as `g__uncultured`, which name a database entry rather than an organism. Unassigned and placeholder names are honest results about the reference database. Draw your figures at genus level.
:::

> **CHECKPOINT 4** — Everyone has `taxonomy.qza` and `barplot.qzv`.

---

## Part 5 — The blanks, the mock, and the y-axis

Open `barplot.qzv` and set the level to 6 (genus). Then collapse the table to genus and export it, so you can check what you see:

```bash
qiime taxa collapse --i-table table.qza --i-taxonomy taxonomy.qza \
  --p-level 6 --o-collapsed-table genus.qza
qiime tools export --input-path genus.qza --output-path genus
biom convert -i genus/feature-table.biom -o genus.tsv --to-tsv
```

### 5.1 What is in the blanks?

```bash
# reads in the blanks, number of real samples containing the genus, genus
awk -F'\t' 'NR==2 {for (i=2;i<=NF;i++) {if ($i ~ /^blank/) b[i]=1; else if ($i ~ /^mock/) m[i]=1}; next}
  NR>2 {bl=0; n=0
        for (i=2;i<=NF;i++) {if (i in b) bl+=$i; else if (!(i in m) && $i>0) n++}
        if (bl>0) {k=split($1, t, ";"); printf "%8d  %3d  %s\n", bl, n, t[k]}}' genus.tsv | sort -k1,1nr | head
```

**Question 7.** Name the genera in your blanks. Several also appear in real samples. What can you conclude about those genera in the real samples?

::: Answer
Expect recognisable reagent genera: *Ralstonia*, *Bradyrhizobium*, *Pseudomonas*, *Sphingomonas* and others from the lecture's list. Where a blank genus also appears in real samples, you cannot yet tell whether it came from the sample or the kit.

This is exactly the situation decontam exists for, and exactly where it cannot be certain. If a blank genus is one of your group differences, it needs independent evidence before it means anything. Contamination is a design problem: the blanks are what make the question answerable at all.
:::

### 5.2 Does the mock look like the mock?

```bash
awk -F'\t' 'NR==FNR {if (FNR>1) e[$1]=$2; next}
  FNR==2 {for (i=2;i<=NF;i++) if ($i ~ /^mock/) c=i; next}
  FNR>2 && c {k=split($1, t, ";"); g=t[k]; sub(/^g__/, "", g); o[g]+=$c; tot+=$c}
  END {for (g in e) printf "%-24s expected %5.1f%%   observed %5.1f%%\n", g, 100*e[g], (tot ? 100*o[g]/tot : 0)}' \
  $AMP/mock_expected.tsv genus.tsv
```

**Question 8.** The mock contains known organisms at known ratios. Where does your observed composition differ from the expected, and what could cause it?

::: Answer
Expect real distortions of several-fold for some members. Three causes, all from the lecture:

- **16S copy number** varies between genomes, so equal numbers of cells give unequal numbers of 16S molecules
- **Primer mismatch** amplifies some taxa more efficiently than others
- **Extraction bias**: some cell walls break more easily than others

None of this is a bug in your pipeline. It is the size of the bias your pipeline has, measured. That is why a mock belongs on every run.
:::

### 5.3 What is on the y-axis?

Take the first real sample, and imagine its most abundant genus quadrupled while nothing else changed at all:

```bash
awk -F'\t' 'NR==2 {for (i=2;i<=NF;i++) if ($i !~ /^(blank|mock)/) {c=i; name=$i; break}; next}
  NR>2 {k=split($1, t, ";"); g[NR]=t[k]; x[NR]=$c; tot+=$c; if ($c>max) {max=$c; top=NR}}
  END {printf "%s: share of reads, before and after its top genus quadruples\n", name
       nt = tot + 3*max
       for (r in x) if (x[r] >= 0.03*tot)
         printf "  %-26s %5.1f%%  ->  %5.1f%%\n", g[r], 100*x[r]/tot, 100*(r==top ? 4*x[r] : x[r])/nt}' genus.tsv
```

**Question 9.** Only one genus changed. What happened to every other genus's share, and what does that mean for a barplot that seems to show a genus "decreasing"?

::: Answer
Every other genus's share fell, though their counts did not change at all. Everything on the barplot sums to one, so one real increase produces an apparent decrease in everything else.

So a genus that looks lower in one group may have decreased, or something else may have increased, and sequencing alone cannot tell you which. That is the compositional problem from the lecture, and it is why the statistics in Session 2 use methods built for it.
:::

> **CHECKPOINT 5** — Everyone can say what their blanks contain, how far the mock is off, and what is on the y-axis.

---

## Takeaways

| Task | Command |
|---|---|
| Import paired reads | `qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'` |
| Remove primers | `qiime cutadapt trim-paired --p-discard-untrimmed` |
| Choose truncation | `qiime demux summarize`, then F + R − 253 ≥ 12 plus a margin |
| Denoise | `qiime dada2 denoise-paired`, as a batch job |
| Check where reads went | export `stats.qza`; read the merged column first |
| Name and plot | `qiime feature-classifier classify-sklearn`, `qiime taxa barplot` |
| Get a plain table | `qiime taxa collapse`, export, `biom convert --to-tsv` |

**Three things to carry forward:**

1. Truncation length is the decision that matters, and it can be checked with arithmetic before you run anything.
2. A blank you sequenced is evidence. A genus in your blanks cannot be a finding until something independent supports it.
3. Everything on a barplot sums to one. One genus increasing makes every other genus look as if it decreased.
