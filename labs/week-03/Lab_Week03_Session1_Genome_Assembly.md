# Lab 3A — Assembling a Genome

**Applied Bioinformatics · Week 3, Session 1**
**Time: 30 minutes · Runs alongside the lecture**

> **Read Part 1 first and submit the assembly immediately.** It runs while you do the rest of the lab. Do not read the whole handout before starting.

---

## Before you start

```bash
srun --pty -p interactive -n 8 --mem=32g --time=0-03:00:00 bash
hostname                      # must NOT say login-p0x

export COURSE=/cluster/tufts/bio_appbio/shared/wk3
export HIFI=~/appbio/wk3/fetch_hifi/results/fastq
export GENOME=15000000        # ← genome size in bases, from the board

mkdir -p ~/appbio/wk3/assembly && cd ~/appbio/wk3/assembly
ls -lh $HIFI
```

If `$HIFI` is empty, your fetchngs job from Lab 2B did not finish:

```bash
export HIFI=$COURSE/fallback_hifi
ls -lh $HIFI
```

---

## Part 1 — Submit the assembly, now (5 minutes)

### 1.1 Set up the environment

```bash
module load anaconda
source activate appbio
conda install -y -c conda-forge -c bioconda hifiasm seqkit compleasm
```

### 1.2 Subsample first

Your HiFi run is roughly 379× coverage. HiFi assembly is comfortable at 20–30×, and extra
depth only makes hifiasm slower and hungrier for memory.

```bash
conda install -y -c conda-forge -c bioconda rasusa

rasusa reads --coverage 30 --genome-size 12.1mb \
  $HIFI/*.fastq.gz -o hifi_30x.fastq.gz

seqkit stats -a $HIFI/*.fastq.gz hifi_30x.fastq.gz
```

**Question 1a.** Why subsample rather than assemble everything? What would 379× buy you?

::: Answer
Nothing useful. Above about 30×, HiFi assembly quality plateaus — the reads are already
long and accurate, so extra copies of the same information do not resolve anything new.

What it costs is real: runtime and memory scale with input, so 379× turns a three-minute
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

module load anaconda
source activate appbio

hifiasm -o asm -t ${SLURM_CPUS_PER_TASK} hifi_30x.fastq.gz

# hifiasm writes GFA; convert the primary contigs to FASTA
awk '/^S/{print ">"$2"\\n"$3}' asm.bp.p_ctg.gfa | fold > asm.p_ctg.fa
EOF

sbatch assemble.sh
squeue -u $USER
```

**Job ID: ______________**

> At 30× on a 12 Mb genome, expect **2–5 minutes**. Check back at Part 3.

**Question 1b.** Why does the script use `${SLURM_CPUS_PER_TASK}` instead of writing `-t 16`?

::: Answer
So the thread count can never disagree with what SLURM allocated. If you later change `--cpus-per-task` and forget to change the hard-coded number, you either waste cores you paid for or oversubscribe the node and slow everything down.

This is the same principle as `${task.cpus}` in Nextflow, which you will meet properly later.
:::

> **CHECKPOINT 1** — Everyone's assembly is submitted and queued. Now we wait, productively.

---

## Part 2 — Interrogate the reads while it runs (9 minutes)

You should know what you are assembling before you look at the assembly.

### 2.1 Basic statistics

```bash
seqkit stats -a $HIFI/*.fastq.gz
```

Record: `num_seqs`, `sum_len`, `avg_len`, `N50`, `Q20(%)`, `Q30(%)`.

### 2.2 Coverage

```bash
seqkit stats -T $HIFI/*.fastq.gz | awk -v g=$GENOME 'NR>1 {s+=$5} END {printf "coverage: %.1fx\n", s/g}'
```

**Question 2.** You have some coverage figure. For HiFi, is it enough? How does the answer differ from what you would need with nanopore or Illumina?

::: Answer
HiFi at 20–30× comfortably assembles a small eukaryote; 15× is workable. That is remarkably low, and it is because HiFi reads are *both* long and accurate — each read is close to trustworthy on its own.

Noisy long reads need far more depth, because consensus across many reads is what removes the error. Illumina needs more still for assembly and cannot span repeats at any depth.

Depth and read length buy different things. This is the Week 2 point again.
:::

### 2.3 The read-length distribution

```bash
seqkit fx2tab -nl $HIFI/*.fastq.gz | awk '{print $NF}' | sort -n | \
  awk '{a[NR]=$1} END {
        printf "min    %d\n", a[1];
        printf "median %d\n", a[int(NR*0.5)];
        printf "N50-ish%d\n", a[int(NR*0.5)];
        printf "max    %d\n", a[NR] }'
```

**Question 3.** HiFi reads are typically 15–25 kb. Which repeats in a fungal genome will these span, and which will they not?

::: Answer
They will span most transposable elements, most gene-family tandem duplications, and — usually — rRNA operons, which are around 8–9 kb.

They will not span long tandem arrays: rDNA repeat clusters that run to hundreds of kilobases, centromeric and subtelomeric satellite. Those are where your assembly will break, and the breaks in your contig set are a map of them.

**Every gap in your assembly is a repeat your reads could not span.** That is the whole of assembly theory in one sentence.
:::

### 2.4 Check on the job

```bash
squeue -u $USER
tail -5 logs/hifiasm_*.out
```

> **CHECKPOINT 2** — Read statistics recorded. Assembly should be running or finished.

---

## Part 3 — Look at what came out (11 minutes)

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

Your Illumina reads from Week 1 are exactly that second technology. Homework 3 does this properly.

This is also why the hybrid design keeps appearing, and why Week 10 will tell you correctness is the question you largely *cannot* answer for a metagenome-assembled genome.
:::

> **CHECKPOINT 3** — Everyone has an assembly with contiguity and completeness numbers.

---

## Part 4 — Fetch the evidence for next session (5 minutes)

Session 2 annotates this genome. Structural annotation of a eukaryote is far better with RNA evidence, so fetch some now.

### 4.1 Check what you are asking for

```bash
export RNA_PROJECT=PRJNA000000        # ← REPLACE THIS with the accession on the board

ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS="run_accession,instrument_platform,library_strategy,library_layout,read_count,base_count"

curl -sS --fail --max-time 120 \
  "${ENA_API}?accession=${RNA_PROJECT}&result=read_run&fields=${FIELDS}&format=tsv" \
  -o runs.tsv

# Header plus at least one run?
[ "$(wc -l < runs.tsv)" -ge 2 ] \
  && echo "OK: $(($(wc -l < runs.tsv) - 1)) runs" \
  || echo "NO DATA — see troubleshooting below"

# Aligned table, without needing `column`
awk -F'\t' '{ for (i=1; i<=NF; i++) printf "%-22s", ($i=="" ? "-" : $i); print "" }' runs.tsv
```

::: If that produced no output or an error
Work through these in order:

1. **`curl: command not found`** — use `wget` with the same URL instead.
2. **`curl: (22) ... error: 400`** — the accession is wrong. Check it against the board; the placeholder is not a real accession.
3. **Empty file, no error** — compute nodes may have no outbound internet. Test with
   a bare `curl` to the ENA API endpoint, discarding the output, to see if it responds.
   If unreachable, run this one command from a **login** node (it is tiny), or use the staged copy:
   `cp $COURSE/runs.tsv .`
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
mkdir -p ~/appbio/wk3/fetch_rna/logs && cd ~/appbio/wk3/fetch_rna
cp $COURSE/ids_rna.csv ids.csv
cat ids.csv

cat > fetch_rna.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=fetch_rna
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load nextflow
module load singularity
export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity_cache

nextflow run nf-core/fetchngs \
    -r 1.12.0 \
    -profile singularity \
    --input ids.csv \
    --outdir results \
    --download_method ftp \
    -resume
EOF

sbatch fetch_rna.sh
squeue -u $USER
```

**Job ID: ______________**

### 4.3 Before Thursday

```bash
cd ~/appbio/wk3/fetch_rna && seff <jobid> && ls -lh results/fastq/
```

Also make sure your assembly survived:

```bash
ls -lh ~/appbio/wk3/assembly/asm.p_ctg.fa
```

**Bring both to Session 2.** Fallbacks are at `$COURSE/fallback_rna/` and `$COURSE/fallback_assembly/`.

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
