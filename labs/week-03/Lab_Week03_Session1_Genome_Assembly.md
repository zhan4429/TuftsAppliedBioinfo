# Week 3 Session 1: Genome Assembly

Applied Bioinformatics - Tufts University Department of Biology

> **Submit both batch jobs in class:** the direct assembly in Part 1, then the pipeline in Part 5. They can run after you leave. Work on the read statistics and RNA-seq preparation while they run; complete the assembly checks and comparison after both jobs finish.

---

## Before You Start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname                      # must NOT say login-p0x

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED=$COURSE/shared/wk3
export FETCH=$MYWORK/week-01/session-02/fetch
export GENOME=12100000        # expected genome size in bases

mkdir -p $MYWORK/week-03/session-01 && cd $MYWORK/week-03/session-01

# the PacBio HiFi reads you cleaned with fastplong in Week 2
export HIFI=$(ls $MYWORK/week-02/session-02/trimmed/*SRR18210286*.filt.fq.gz 2>/dev/null)
ls -lh "$HIFI"
```

These are the same paths as Weeks 1 and 2: your work under `$MYWORK`, shared course files under `$COURSE/shared`. If the allocation announced in class differs, change `COURSE` to match.

nf-core/genomeassembler runs the same two steps in the same order: `fastplong` on the reads, then `hifiasm`.

If `$HIFI` is empty, you did not finish Week 2 Part 2. Use the reads you downloaded in Week 1 instead. In Week 2, `fastplong` removed almost nothing from the HiFi reads, so the assembly will be essentially the same:

```bash
export HIFI=$(ls $FETCH/results/fastq/*SRR18210286*.fastq.gz 2>/dev/null)
ls -lh "$HIFI"
```

If that is empty too, your Week 1 download did not finish. Use the staged copy:

```bash
export HIFI=$(ls $SHARED/fallback_hifi/*.fastq.gz | head -1)
ls -lh $HIFI
```

---

## Part 1: Submit The Assembly, Now

### 1.1 Load The Tools

```bash
module load hifiasm seqkit compleasm
module list
```

Every tool in this course is an environment module. `module avail <tool>` shows the versions installed; write down the ones you used.

### 1.2 Use All HiFi Reads

Use all the cleaned HiFi reads selected in `$HIFI`. This run contains about 4.58
billion bases for a 12.1 Mb genome, or roughly 379× coverage. Subsampling can reduce
runtime and memory use, but it is not required for this lab: the batch job can run
after class.

Create a link with a short, stable name. This points to the original file without
copying the sequencing data:

```bash
ln -sf "$HIFI" hifi.fastq.gz
ls -lhL hifi.fastq.gz
```

Both assembly routes will use this same input. More coverage does not guarantee a
better assembly or resolve every repeat; you will evaluate the results rather than
assuming that extra reads make them correct.

### 1.3 Submit It

```bash
mkdir -p logs

cat > assemble.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# load what the job needs inside the script, so it never depends on the shell that submitted it
module load hifiasm
set -euo pipefail

hifiasm -i -o asm -t ${SLURM_CPUS_PER_TASK} hifi.fastq.gz

# hifiasm writes GFA; convert the primary contigs to FASTA
awk '/^S/{print ">"$2"\n"$3}' asm.bp.p_ctg.gfa | fold > asm.p_ctg.fa
EOF

bash -n assemble.sh
sbatch assemble.sh
squeue --me
```

**Job ID: ______________**

Two lines in that script are worth understanding. `module load hifiasm` sits inside the script even though you loaded it interactively: `sbatch` passes on the environment of the shell you submit from, so without it the job would work from one shell and fail from a fresh login. And because the heredoc delimiter is quoted, the script is written exactly as shown, so awk receives `"\n"` and prints a real newline between each contig name and its sequence.

`-i` tells hifiasm to ignore earlier `.bin` cache files with the same output prefix.
This ensures that an earlier trial with fewer reads is not reused for the full-read
assembly. It recomputes the assembly without requiring you to delete those files.

> Check back at Part 3.

**Question 1.** Why does the script use `${SLURM_CPUS_PER_TASK}` instead of writing `-t 16`?

<details>
<summary>Answer</summary>

So the thread count can never disagree with what SLURM allocated. If you later change `--cpus-per-task` and forget to change the hard-coded number, you either waste cores you paid for or oversubscribe the node and slow everything down.

This is the same principle as `${task.cpus}` in Nextflow, which you will meet properly later.

</details>

**Checkpoint 1.** Everyone's assembly is submitted and queued. Now we wait, productively.

---

## Part 2: Interrogate The Reads While It Runs

You should know what you are assembling before you look at the assembly.

### 2.1 Basic Statistics

```bash
seqkit stats -a $HIFI
```

Record: `num_seqs`, `sum_len`, `avg_len`, `N50`, `Q20(%)`, `Q30(%)`.

### 2.2 Coverage

```bash
seqkit stats -T $HIFI | awk -v g=$GENOME 'NR>1 {s+=$5} END {printf "coverage: %.1fx\n", s/g}'
```

**Question 2.** You have some coverage figure. For HiFi, is it enough? How does the answer differ from what you would need with nanopore or Illumina?

<details>
<summary>Answer</summary>

The roughly 379× average depth is ample for attempting this small yeast assembly.
It does not establish that every region has the same coverage or that the result
will be complete and correct.

Coverage requirements depend on read accuracy and length, repeats, heterozygosity,
and the assembler. HiFi reads combine long sequences with high per-base accuracy.
Noisier reads need more consensus support, while short reads may be unable to bridge
long repeats even at high depth. Coverage and read length solve different problems;
there is no universal 30× cutoff above which additional data can never help.

</details>

### 2.3 The Read-length Distribution

```bash
seqkit fx2tab -nl $HIFI | awk '{print $NF}' | sort -n | \
  awk '{a[NR]=$1; s+=$1} END {
        printf "min    %d\n", a[1];
        printf "median %d\n", a[int((NR+1)/2)];
        for (i=NR; i>0; i--) { c+=a[i]; if (c>=s/2) { printf "N50    %d\n", a[i]; break } }
        printf "max    %d\n", a[NR] }'
```

**Question 3.** The HiFi reads are mostly tens of kilobases long. Which repeat structures in *S. cerevisiae* could an individual read span, and which are much longer than a typical read?

<details>
<summary>Answer</summary>

A read long enough to include a repeat plus unique sequence on both sides can place that repeat unambiguously. Many HiFi reads can span a roughly 6 kb Ty element or a single roughly 9 kb rDNA repeat unit.

Spanning one rDNA unit is different from spanning the entire tandem rDNA array, which contains many copies. That array and some subtelomeric repeat structures can remain unresolved. Budding yeast has short point centromeres, not the long centromeric satellite arrays found in many other eukaryotes.

Repeat length, similarity among copies and the available flanking sequence all matter. Repeats are a common cause of assembly breaks, but missing coverage, contamination and other assembly problems can also create breaks.

</details>

### 2.4 Check On The Job

```bash
squeue --me
tail -5 logs/hifiasm_*.out
```

**Checkpoint 2.** Read statistics recorded and the assembly job checked. It may still be queued or running.

---

## Part 3: Look At What Came Out

Run this part after the direct assembly job finishes. If it is still queued or
running during class, finish Part 4 and confirm that you submitted Part 5, then
return here after class.
Do not replace the running job's outputs with a staged assembly.

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

<details>
<summary>Answer</summary>

Usually unresolved heterozygosity. hifiasm is a phasing assembler — for a diploid it tries to keep haplotypes apart rather than collapsing them. If the two haplotypes are divergent enough, both end up in the primary assembly and total length inflates toward 2×.

`asm.bp.hap1.p_ctg.gfa` and `hap2` are the separated haplotypes; `asm.bp.p_ctg.gfa` is the primary set. Which you want depends on your question — and for a haploid organism this should not arise at all.

Substantially *smaller* than expected means you lost sequence: too little coverage, or repeats collapsed.

</details>

### 3.2 Look At The Graph, Not Just The FASTA

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

<details>
<summary>Answer</summary>

Unresolved structure. Each link is a junction the assembler could not resolve into a single path, so it stopped and left the decision visible in the graph.

A clean assembly is mostly long segments with few links. A tangle is many short segments densely connected — which is what a short-read assembly of the same genome would look like.

</details>

### 3.3 Completeness

```bash
compleasm run -a asm.p_ctg.fa -o compleasm -l fungi_odb10 -t 8
cat compleasm/summary.txt
```

Record the compleasm version and the **actual lineage ID and marker count** in its summary, not just the requested command. Some releases can translate an `odb10` request to an `odb12` lineage. If you use the staged cache with `-L "$SHARED/mb_downloads"`, check which lineage version it contains. The installed module and cache versions need to be confirmed on Pax before comparing scores.

Record `S` (single), `D` (duplicated), `F` (fragmented), `M` (missing).

**Question 6.** Your `D` is elevated — say 15% rather than under 2%. What does that suggest, and is it a problem?

<details>
<summary>Answer</summary>

Unresolved haplotypes. Each locus appears twice, so each single-copy marker is found twice and gets scored as duplicated.

Whether it is a problem depends on your goal. For a haplotype-resolved assembly it is the intended result. For a reference you plan to annotate and map against, it is a problem — you would purge duplicates or use one haplotype set.

This is the same signal as an inflated total length in §3.1, seen through a different instrument. Two metrics agreeing is how you become confident in a diagnosis.

</details>

### 3.4 The Three Questions

You have now measured two of the three dimensions from the lecture:

| Question | Metric | Your answer |
|---|---|---|
| **Contiguity** — how few pieces? | contigs, N50 | |
| **Completeness** — is the content there? | compleasm S% | |
| **Correctness** — are the bases right? | Merqury QV | *(homework)* |

**Question 7.** Why can you not measure correctness right now, and what would you need?

<details>
<summary>Answer</summary>

Merqury evaluates assembly k-mers against a database of accurate genomic read k-mers. An independent sequencing technology can provide a useful check on errors that the assembly reads and assembler share, provided coverage and sample matching are adequate.

Use BGISEQ genomic run **SRR17374239**, rather than Illumina RNA-seq. RNA-seq does not uniformly cover intergenic DNA, so it cannot provide a whole-genome k-mer reference.

The BGISEQ and HiFi BioSamples both say **S288C**, but they have different identifiers and collection dates. They are not documented as the same DNA preparation. A biological difference between the cultures can also produce unsupported assembly k-mers. Report that limitation when interpreting the Homework 3 QV; it is not a proven assembly-error rate.

</details>

**Checkpoint 3.** Everyone has an assembly with contiguity and completeness numbers.

---

## Part 4: Prepare The Evidence For Next Session

Session 2 uses RNA-seq as evidence for gene annotation. You already downloaded the
Illumina RNA-seq run **SRR27956204** in Week 1 and inspected it in Week 2; reuse those
files rather than submitting another download job.

### 4.1 Check The RNA-seq Files

```bash
export RNA="$FETCH/results/fastq"
export RNA_R1="$RNA/SRX23612877_SRR27956204_1.fastq.gz"
export RNA_R2="$RNA/SRX23612877_SRR27956204_2.fastq.gz"
ls -lh "$RNA_R1" "$RNA_R2"
```

The `_1` and `_2` files are the two mates of the same paired-end RNA-seq run. The
BGISEQ run, **SRR17374239**, is genomic DNA and must not be substituted for these
RNA-seq reads. Both runs are in the same download folder, so select the filenames
explicitly.

[The ENA record for SRR27956204](https://www.ebi.ac.uk/ena/browser/view/SRR27956204)
identifies it as paired-end Illumina RNA-seq from *Saccharomyces cerevisiae*. You do
not need a second fetchngs job or a new configuration file.

**Question 8.** Why does RNA-seq help annotate a genome, and what can it not tell you?

<details>
<summary>Answer</summary>

RNA-seq reads provide evidence for transcribed regions. Reads that span splice
junctions help locate exon boundaries and introns.

A gene that was not expressed in the sampled condition may have no RNA-seq support.
Absence of coverage does not establish absence of a gene, so annotation also uses
sequence-based predictions and protein evidence.

</details>

### 4.2 Before Session 2

Check the RNA-seq paths during class. After the assembly job finishes, also check that the assembly exists; complete that check after class if the job is still running:

```bash
ls -lh "$MYWORK/week-03/session-01/asm.p_ctg.fa"
ls -lh "$RNA_R1" "$RNA_R2"
```

Keep the Week 1 `results/fastq` files. Cleaning a completed Nextflow run's intermediate
`work/` directory, as taught in Week 1, is different from deleting its published
`results/` files; later labs still need these reads.

If the RNA files are missing, check the Week 1 download or ask for the staged copies
in `$SHARED/fallback_rna/`. The staged assembly is in `$SHARED/fallback_assembly/`.

**Checkpoint 4.** Both RNA-seq mate files located and their paths recorded. Verify the completed assembly before Session 2.

---

## Part 5: Submit The Same Assembly As A Pipeline

Submit this second batch job once during class, immediately after Part 1. It can run
alongside the direct assembly because it writes into its own `genomeassembler/`
directory. You do not need the direct assembly to finish before submitting it.

nf-core/genomeassembler runs `fastplong` on the reads, assembles them with `hifiasm`,
measures the result with QUAST and BUSCO, and writes a report. Give it the same full
HiFi input as the direct assembly, then compare the two results after class.

Session 2 uses the assembly you made by hand, so Thursday does not depend on this run.

### 5.1 Write The Samplesheet

```bash
export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
mkdir -p $MYWORK/week-03/session-01/genomeassembler/logs
cd $MYWORK/week-03/session-01/genomeassembler

cat > samplesheet.csv <<EOF
sample,hifireads,hifi_fastplong_args
yeast_hifi,$MYWORK/week-03/session-01/hifi.fastq.gz,--mean_qual 10 --length_required 1000 --unqualified_percent_limit 100
EOF
cat samplesheet.csv
```

The pipeline reads one row per sample. `hifireads` points to the same full HiFi file used by the direct assembly job. `hifi_fastplong_args` gives `fastplong` the settings you chose in Week 2; without that column, the pipeline would run `fastplong` with its defaults.

This heredoc is not quoted, unlike the job scripts. Here you want `$MYWORK` replaced by your real path as the file is written, because the pipeline needs full paths. Check that `cat` shows a path starting with `/cluster`.

### 5.2 Submit The Pipeline

Create a small parameter file in the same `genomeassembler/` directory:

```bash
cat > genomeassembler_params.json <<'EOF'
{
  "busco": true,
  "quast": true,
  "busco_lineage": "fungi_odb10"
}
EOF
cat genomeassembler_params.json
```

BUSCO and QUAST are off by default in pipeline version 2.0.1. These settings turn
them on. JSON preserves the unquoted `true` values as booleans; `"true"` would be
text. With newer Nextflow parameter parsing, `--busco true` and even a bare
`--busco` can arrive as strings and fail the pipeline's validation. Use
`-params-file` for these boolean settings, following the
[Nextflow parameter-file guidance](https://training.nextflow.io/latest/nfcore_use/02_configure_execution/).

```bash
cat > run_genomeassembler.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=genomeassembler
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load nextflow singularity

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
cd "$MYWORK/week-03/session-01/genomeassembler"

export NXF_SINGULARITY_CACHEDIR="$MYWORK/.singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

nextflow run nf-core/genomeassembler \
  -r 2.0.1 \
  -profile tufts \
  -params-file genomeassembler_params.json \
  --input samplesheet.csv \
  --outdir results \
  -resume
EOF

bash -n run_genomeassembler.sh
sbatch run_genomeassembler.sh
squeue --me
```

**Job ID: ______________**

On a first run, Nextflow may warn that `-resume` is ignored because there is no
previous execution to resume. That warning does not mean parameter validation
failed. Keep `-resume` for subsequent submissions after fixing a failed run.

This job is small on purpose. In Week 1, `-profile singularity` ran every step inside one large job. `-profile tufts` loads the Tufts cluster settings instead, so Nextflow submits each step as its own SLURM job, with the CPUs and memory that step asks for. The job you submitted only runs Nextflow, which hands out the work and waits. Once it starts, `squeue --me` also lists the jobs it submits; their names begin with `nf-`.

`"busco_lineage": "fungi_odb10"` in the JSON file requests an OrthoDB 10 lineage. Compare its actual lineage ID and marker count with the compleasm summary from Section 3.3. Matching parameter values do not establish that both tools used the same dataset. If they differ, report that difference and defer a direct score comparison until both can be run against the same lineage release. The Tufts profile deletes `work/` automatically when a run succeeds, so there is no clean-up step this time. If the run fails, `work/` stays: read `logs/genomeassembler_*.err` and `.nextflow.log`, fix the problem, and resubmit the same script, and `-resume` reuses the steps that finished.

Return to the main lab directory before continuing with Parts 2–4:

```bash
cd "$MYWORK/week-03/session-01"
```

**Checkpoint 5.** Both batch jobs are submitted and both job IDs are recorded.
Neither job needs to finish before you leave class.

---

## After Class: Compare The Two Assemblies

Work on a compute node. If this is a new login, repeat the environment setup in
Before You Start and load the tools from Section 1.1.

Wait until both jobs have finished successfully. Check their states with `sacct`,
replacing the example IDs with your direct-assembly and pipeline job IDs:

```bash
sacct -j 12345,67890 --format=JobID,JobName,State,ExitCode
```

If a job failed, inspect its logs and resolve the failure before comparing results.
Complete Part 3 on your direct assembly if you have not already done so, then inspect
the pipeline outputs:

```bash
cd $MYWORK/week-03/session-01/genomeassembler
tail -n 5 logs/genomeassembler_*.out
ls results results/yeast_hifi

cat results/yeast_hifi/QC/QUAST/yeast_hifi_assembly/report.txt
cat results/yeast_hifi/QC/BUSCO/short_summary.*.txt
```

The last lines of the log should say the pipeline completed successfully. Open `results/report/report.html` through OnDemand Files for the pipeline's summary. `results/pipeline_info/` records what ran, which software versions were used, and how long each step took.

**Question 9.** Compare QUAST's number of contigs, total length and N50 with your §3.1 values, and BUSCO's `S` and `D` with compleasm's from §3.3. Should they be identical?

<details>
<summary>Answer</summary>

Close, but not necessarily identical, and each difference has a reason you can check.

The pipeline runs QUAST in its large-genome mode, which skips contigs shorter than 3,000 bp. `seqkit stats` counted every sequence, so QUAST can report fewer contigs and a slightly smaller total. The pipeline also ran `hifiasm` from its own container, which may be a different version from the cluster module: compare the version in `results/pipeline_info/nf_core_genomeassembler_software_versions.yml` with `hifiasm --version`.

First compare the actual lineage release and number of marker genes in both summaries. Only if those agree are you comparing the same marker set with different software. If they differ, a score difference mixes software and database effects and cannot be attributed to the assembly or method alone.

A large disagreement would be worth chasing. Two routes to the same assembly should not tell different stories.

</details>

**Question 10.** What did the pipeline do that you did not do by hand, and what did you do that the pipeline did not?

<details>
<summary>Answer</summary>

The pipeline ran `fastplong` on the reads again, mapped the reads back to the assembly for QUAST, ran BUSCO, wrote a report, and recorded every command and software version. Each step got a job sized for that step, and an interrupted run can resume where it stopped.

You inspected the input reads, calculated their coverage, and examined the assembly
graph as well as the FASTA. Both routes received the same full input, but your
interpretation of the metrics is still needed to judge either result.

A pipeline makes a method repeatable. It does not make it understood.

</details>

---

## Takeaways

Assemble HiFi reads, then turn the graph into a FASTA:

```bash
hifiasm -i -o asm -t 16 reads.fastq.gz
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

1. Repeats are a common cause of assembly breaks. Read length, coverage and sequence
   quality all affect whether a repeat can be resolved.
2. Excess assembly length and duplicated markers can support an unresolved-haplotype
   explanation, but also check true duplication, contamination and the ploidy assumption.
3. Contiguity, completeness and correctness are three independent questions. You measured
   two today; the third needs a second technology.

## Dataset References

- [BGISEQ BioSample: S288C, collected 2021-01-21](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN24477750)
- [HiFi BioSample: S288C, collected 2021-12-29](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN24477612)
- [compleasm releases and lineage compatibility](https://github.com/huangnengCSU/compleasm/releases)
