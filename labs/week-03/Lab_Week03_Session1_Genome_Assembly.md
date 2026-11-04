# Lab 3A — Assembling a Genome

**Applied Bioinformatics · Week 3, Session 1**
**Runs alongside the lecture**

> **Read Part 1 first and submit the assembly immediately.** It runs while you do the rest of the lab. Do not read the whole handout before starting.

---

## Before you start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname                      # must NOT say login-p0x

export COURSE=/cluster/tufts/yzhang85/bioinfo_course/appliedBio
export MYWORK="$COURSE/$USER"
export SHARED=/cluster/tufts/bio_appbio/shared/wk3
export FETCH=$MYWORK/week-01/session-02/fetch
export GENOME=12100000        # expected genome size in bases

mkdir -p $MYWORK/week-03/session-01 && cd $MYWORK/week-03/session-01

# the PacBio HiFi run you fetched in Week 1
export HIFI=$(ls $FETCH/results/fastq/*SRR18210286*.fastq.gz 2>/dev/null)
ls -lh "$HIFI"
```

These are the Week 1 paths: your work under `$COURSE`, shared course files in the separate `bio_appbio` allocation. If the allocation announced in class differs, change `COURSE` to match.

If `$HIFI` is empty, your Week 1 fetchngs download did not finish. Use the staged copy:

```bash
export HIFI=$(ls $SHARED/fallback_hifi/*.fastq.gz | head -1)
ls -lh $HIFI
```

---

## Part 1 — Submit the assembly, now

### 1.1 Set up the environment

```bash
module load miniforge
conda create -y -p $MYWORK/envs/wk3 -c conda-forge -c bioconda hifiasm seqkit compleasm rasusa
conda activate $MYWORK/envs/wk3
```

This creates one environment for all of Week 3 under `$MYWORK`, where Week 1 told you environments belong. Lab 3B adds to it.

### 1.2 Subsample first

Your HiFi run is roughly 379× coverage. HiFi assembly is comfortable at 20–30×, and extra
depth only makes hifiasm slower and hungrier for memory.

```bash
rasusa reads --coverage 30 --genome-size $GENOME \
  $HIFI -o hifi_30x.fastq.gz

seqkit stats -a $HIFI hifi_30x.fastq.gz
```

**Question 1a.** Why subsample rather than assemble everything? What would 379× buy you?

::: Answer
Nothing useful. Above about 30×, HiFi assembly quality plateaus — the reads are already
long and accurate, so extra copies of the same information do not resolve anything new.

What it costs is real: runtime and memory scale with input, so 379× turns a quick
assembly into a much longer one for an identical result.

Knowing when more data stops helping is a genuinely useful instinct, and it is the same
reasoning as the coverage arithmetic in Week 2.
:::

### 1.3 Submit it

```bash
mkdir -p logs

cat > assemble.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# a batch job starts a fresh, non-interactive shell: define paths and set up conda inside the script
export COURSE=/cluster/tufts/yzhang85/bioinfo_course/appliedBio
export MYWORK="$COURSE/$USER"
module load miniforge
eval "$(conda shell.bash hook)"
conda activate $MYWORK/envs/wk3
set -euo pipefail

hifiasm -o asm -t ${SLURM_CPUS_PER_TASK} hifi_30x.fastq.gz

# hifiasm writes GFA; convert the primary contigs to FASTA
awk '/^S/{print ">"$2"\n"$3}' asm.bp.p_ctg.gfa | fold > asm.p_ctg.fa
EOF

bash -n assemble.sh
sbatch assemble.sh
squeue -u $USER
```

**Job ID: ______________**

Two lines in that script are worth understanding. `eval "$(conda shell.bash hook)"` makes `conda activate` work in a batch job, which does not run the setup your login shell does; it is the Week 1 lesson that a job does not inherit your session. And because the heredoc delimiter is quoted, the script is written exactly as shown, so awk receives `"\n"` and prints a real newline between each contig name and its sequence.

> Check back at Part 3.

**Question 1b.** Why does the script use `${SLURM_CPUS_PER_TASK}` instead of writing `-t 16`?

::: Answer
So the thread count can never disagree with what SLURM allocated. If you later change `--cpus-per-task` and forget to change the hard-coded number, you either waste cores you paid for or oversubscribe the node and slow everything down.

This is the same principle as `${task.cpus}` in Nextflow, which you will meet properly later.
:::

> **CHECKPOINT 1** — Everyone's assembly is submitted and queued. Now we wait, productively.

---

## Part 2 — Interrogate the reads while it runs

You should know what you are assembling before you look at the assembly.

### 2.1 Basic statistics

```bash
seqkit stats -a $HIFI
```

Record: `num_seqs`, `sum_len`, `avg_len`, `N50`, `Q20(%)`, `Q30(%)`.

### 2.2 Coverage

```bash
seqkit stats -T $HIFI | awk -v g=$GENOME 'NR>1 {s+=$5} END {printf "coverage: %.1fx\n", s/g}'
```

**Question 2.** You have some coverage figure. For HiFi, is it enough? How does the answer differ from what you would need with nanopore or Illumina?

::: Answer
HiFi at 20–30× comfortably assembles a small eukaryote; 15× is workable. That is remarkably low, and it is because HiFi reads are *both* long and accurate — each read is close to trustworthy on its own.

Noisy long reads need far more depth, because consensus across many reads is what removes the error. Illumina needs more still for assembly and cannot span repeats at any depth.

Depth and read length buy different things. This is the Week 2 point again.
:::

### 2.3 The read-length distribution

```bash
seqkit fx2tab -nl $HIFI | awk '{print $NF}' | sort -n | \
  awk '{a[NR]=$1; s+=$1} END {
        printf "min    %d\n", a[1];
        printf "median %d\n", a[int((NR+1)/2)];
        for (i=NR; i>0; i--) { c+=a[i]; if (c>=s/2) { printf "N50    %d\n", a[i]; break } }
        printf "max    %d\n", a[NR] }'
```

**Question 3.** HiFi reads are typically 15–25 kb. Which repeats in a fungal genome will these span, and which will they not?

::: Answer
They will span most transposable elements, most gene-family tandem duplications, and — usually — rRNA operons, which are around 8–9 kb.

They will not span long tandem arrays: the rDNA array, which runs to hundreds of kilobases or more, and centromeric and subtelomeric satellite. Those are where your assembly will break, and the breaks in your contig set are a map of them.

**Every gap in your assembly is a repeat your reads could not span.** That is the whole of assembly theory in one sentence.
:::

### 2.4 Check on the job

```bash
squeue -u $USER
tail -5 logs/hifiasm_*.out
```

> **CHECKPOINT 2** — Read statistics recorded. Assembly should be running or finished.

---

## Part 3 — Look at what came out

### 3.1 Contiguity

```bash
ls -lh asm*
seqkit stats -a asm.p_ctg.fa
```

Fill in:

| Metric | Your value | Expected |
|---|---|---|
| Number of contigs | | tens, not thousands |
| Total length | | close to `$GENOME` |
| Largest contig | | megabases |
| N50 | | megabases |

**Question 4.** Compare total assembly length against the expected genome size. What does it mean if yours is substantially *larger*?

::: Answer
Usually unresolved heterozygosity. hifiasm is a phasing assembler — for a diploid it tries to keep haplotypes apart rather than collapsing them. If the two haplotypes are divergent enough, both end up in the primary assembly and total length inflates toward 2×.

`asm.bp.hap1.p_ctg.gfa` and `hap2` are the separated haplotypes; `asm.bp.p_ctg.gfa` is the primary set. Which you want depends on your question — and for a haploid organism this should not arise at all.

Substantially *smaller* than expected means you lost sequence: too little coverage, or repeats collapsed.
:::

### 3.2 Look at the graph, not just the FASTA

```bash
grep -c '^S' asm.bp.p_ctg.gfa      # segments
grep -c '^L' asm.bp.p_ctg.gfa      # links between them
```

A FASTA is the path the assembler chose. The GFA is everything it considered, including the ambiguities it had to break rather than resolve.

If you have Bandage available, render it:

```bash
Bandage image asm.bp.p_ctg.gfa asm_graph.png --height 1200 2>/dev/null \
  || echo "Bandage not installed — view the GFA in Bandage on your laptop instead"
```

**Question 5.** Many links relative to segments means what?

::: Answer
Unresolved structure. Each link is a junction the assembler could not resolve into a single path, so it stopped and left the decision visible in the graph.

A clean assembly is mostly long segments with few links. A tangle is many short segments densely connected — which is what a short-read assembly of the same genome would look like.
:::

### 3.3 Completeness

```bash
compleasm run -a asm.p_ctg.fa -o compleasm -l fungi_odb10 -t 8
cat compleasm/summary.txt
```

compleasm downloads the lineage on first use. If that fails because the compute node has no internet, add `-L $SHARED/mb_downloads` to use the staged copy.

Record `S` (single), `D` (duplicated), `F` (fragmented), `M` (missing).

**Question 6.** Your `D` is elevated — say 15% rather than under 2%. What does that suggest, and is it a problem?

::: Answer
Unresolved haplotypes. Each locus appears twice, so each single-copy marker is found twice and gets scored as duplicated.

Whether it is a problem depends on your goal. For a haplotype-resolved assembly it is the intended result. For a reference you plan to annotate and map against, it is a problem — you would purge duplicates or use one haplotype set.

This is the same signal as an inflated total length in §3.1, seen through a different instrument. Two metrics agreeing is how you become confident in a diagnosis.
:::

### 3.4 The three questions

You have now measured two of the three dimensions from the lecture:

| Question | Metric | Your answer |
|---|---|---|
| **Contiguity** — how few pieces? | contigs, N50 | |
| **Completeness** — is the content there? | compleasm S% | |
| **Correctness** — are the bases right? | Merqury QV | *(homework)* |

**Question 7.** Why can you not measure correctness right now, and what would you need?

::: Answer
You need reads from a second, more accurate technology — or at least reads more accurate than the assembly — to build a k-mer set to compare against. That is what Merqury does: k-mers present in the assembly but absent from the reads have no evidential support and are counted as errors.

Your BGISEQ reads from Week 1 (SRR17374239) are exactly that second technology. Not the Illumina run: it is RNA-seq, whose k-mers cover only transcribed sequence, so every correct intergenic base would look like an error. Homework 3 does this properly.

This is also why the hybrid design keeps appearing, and why Week 10 will tell you correctness is the question you largely *cannot* answer for a metagenome-assembled genome.
:::

> **CHECKPOINT 3** — Everyone has an assembly with contiguity and completeness numbers.

---

## Part 4 — Fetch the evidence for next session

Session 2 annotates this genome. Structural annotation of a eukaryote is far better with RNA evidence, so fetch some now.

### 4.1 Check what you are asking for

```bash
export RNA_PROJECT=PRJNA000000        # ← REPLACE THIS with the accession on the board

ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS="run_accession,instrument_platform,library_strategy,library_layout,read_count,base_count"

curl -sS --fail --max-time 120 \
  "${ENA_API}?accession=${RNA_PROJECT}&result=read_run&fields=${FIELDS}&format=tsv" \
  -o rna_runs.tsv

# Header plus at least one run?
[ "$(wc -l < rna_runs.tsv)" -ge 2 ] \
  && echo "OK: $(($(wc -l < rna_runs.tsv) - 1)) runs" \
  || echo "NO DATA — see troubleshooting below"

# Aligned table, without needing `column`
awk -F'\t' '{ for (i=1; i<=NF; i++) printf "%-22s", ($i=="" ? "-" : $i); print "" }' rna_runs.tsv
```

::: If that produced no output or an error
Work through these in order:

1. **`curl: command not found`** — use `wget` with the same URL instead.
2. **`curl: (22) ... error: 400`** — the accession is wrong. Check it against the board; the placeholder is not a real accession.
3. **Empty file, no error** — compute nodes may have no outbound internet. Test with
   a bare `curl` to the ENA API endpoint, discarding the output, to see if it responds.
   If unreachable, run this one command from a **login** node (it is tiny), or use the staged copy:
   `cp $SHARED/runs.tsv rna_runs.tsv`
4. **Times out** — ENA is busy. Re-run; it is safe to repeat.
:::

> **Why `-sS --fail` rather than `-s`?** `curl -s` is silent about everything, including a 404 — you
> would get an empty file and no explanation. `-sS --fail` stays quiet on success and complains
> loudly on failure. And `column` is not installed everywhere, and the bash tab-escape syntax is a bash-ism that
> silently produces the wrong separator under `sh`; `awk` with an explicit field separator avoids both problems.

`library_strategy` should say `RNA-Seq`.

**Question 8.** Why does RNA-seq help annotate a genome, and what can it not tell you?

::: Answer
RNA-seq reads come from transcripts, so where they align tells you directly where exons are and — from reads spanning junctions — where introns are. That is evidence, not prediction, and it is the single largest improvement available to a eukaryotic annotation.

What it cannot do: reveal genes that were not expressed in the tissue or condition sequenced. A gene silent in your sample is invisible to the evidence, and only the *ab initio* model will find it. This is why modern pipelines combine both.
:::

### 4.2 Submit it

```bash
mkdir -p $MYWORK/week-03/session-01/fetch_rna/logs
cd $MYWORK/week-03/session-01/fetch_rna
cp $SHARED/ids_rna.csv ids.csv
cat ids.csv

cat > fetch_rna.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=fetch_rna
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=18G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load nextflow singularity

export COURSE=/cluster/tufts/yzhang85/bioinfo_course/appliedBio
export MYWORK="$COURSE/$USER"
cd "$MYWORK/week-03/session-01/fetch_rna"

export NXF_HOME="$SCRATCH/.nextflow"
export NXF_SINGULARITY_CACHEDIR="$SCRATCH/.singularity_cache"
export SINGULARITY_CACHEDIR="$MYWORK/.singularity/oci"
export NXF_OPTS='-Xms512m -Xmx4g'
mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR"

nextflow run nf-core/fetchngs \
  -r 1.13.0 \
  -profile singularity \
  -c "$MYWORK/week-01/session-02/lab.config" \
  --input ids.csv \
  --outdir results \
  --download_method sratools \
  -resume
EOF

bash -n fetch_rna.sh
RNA_JOB_ID=$(sbatch --parsable fetch_rna.sh)
printf '%s\n' "$RNA_JOB_ID" > fetch_rna_job_id.txt
squeue -j "$RNA_JOB_ID"
```

This is your Week 1 download job with three changes: the job name, the directory, and the list of runs. It reuses the `lab.config` you wrote in Week 1, so the same resource limits apply.

**Job ID: ______________** (also saved in `fetch_rna_job_id.txt`)

### 4.3 Before Thursday

```bash
cd $MYWORK/week-03/session-01/fetch_rna
seff "$(cat fetch_rna_job_id.txt)"
ls -lh results/fastq/
```

Also make sure your assembly survived:

```bash
ls -lh $MYWORK/week-03/session-01/asm.p_ctg.fa
```

**Bring both to Session 2.** Fallbacks are at `$SHARED/fallback_rna/` and `$SHARED/fallback_assembly/`.

> **CHECKPOINT 4** — Assembly done, RNA-seq job submitted, both job IDs written down.

---

## Takeaways

Assemble HiFi reads, then turn the graph into a FASTA:

```bash
hifiasm -o asm -t 16 reads.fastq.gz
awk '/^S/{print ">"$2"\n"$3}' asm.bp.p_ctg.gfa | fold > asm.p_ctg.fa
```

Measure the three dimensions:

```bash
seqkit stats -a asm.p_ctg.fa                 # contiguity
compleasm run -a asm.p_ctg.fa -l LINEAGE_odb10   # completeness
# correctness needs Merqury and a second, more accurate read set
```

Coverage is always base count divided by genome size. Inspect the graph with segment and
link counts, or with Bandage.

Three things to carry forward:

1. Every gap in your assembly is a repeat your reads could not span. Read length, not
   depth, decides contiguity.
2. Assembly length above expectation and elevated BUSCO duplication are the same finding,
   unresolved haplotypes, seen two ways.
3. Contiguity, completeness and correctness are three independent questions. You measured
   two today; the third needs a second technology.
