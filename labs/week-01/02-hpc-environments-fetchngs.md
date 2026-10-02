# Week 1 Session 2: HPC, Environments, and Fetching Public Data

Applied Bioinformatics - Tufts University Department of Biology

In this lab you will move from the login node to a compute node, submit a small SLURM batch job, create a reproducible software environment, test Nextflow and container support, and start a real public-data download with `nf-core/fetchngs`.

This lab produces the sequencing data used in Week 2 and Week 3. You will fetch one run from each of the four platforms in a single study, which is what makes next week's QC comparison possible. The download will finish after class, so do not skip the job submission and verification steps.

## Before You Start

Open a terminal on Pax:

- OnDemand: <https://ondemand-prod.pax.tufts.edu/> > **Clusters > Tufts HPC Shell Access**
- SSH: `ssh your_utln@login-prod.pax.tufts.edu`

Set paths and make a work directory:

```bash
export COURSE=/cluster/tufts/bio_appbio
export WEEK1_DATA=$COURSE/shared/wk1
export MYWORK=$COURSE/$USER

mkdir -p "$MYWORK/week-01/session-02"
cd "$MYWORK/week-01/session-02"
pwd
```

Three paths, and you will use all of them every week:

- `$COURSE` — the course directory, shared by everyone in the class
- `$WEEK1_DATA` — this week's staged data, read-only
- `$MYWORK` — your own directory inside the course allocation, named after your UTLN

**Do all your coursework under `$MYWORK`, not in your home directory.** Home has a small
quota that sequencing data will exhaust immediately, and the course allocation is where
the instructor can see and help with your work.

## Part 1: Get Off the Login Node

### 1.1 See Where You Are

```bash
hostname
```

If the hostname starts with `login`, you are on a login node shared by everyone on the cluster. Use login nodes for light work only.

### 1.2 Inspect the Cluster

```bash
sinfo -s
```

You may see partitions such as:

```text
PARTITION    AVAIL  TIMELIMIT   NODES(A/I/O/T) NODELIST
batch*          up 2-00:00:00      46/21/23/90 pax[001-002,004-006,012-019,027-033,035-046,054-061,069,072-075,078-085,107-109,120,145-149,154-171,188-197]
gpu             up 2-00:00:00       14/12/1/27 pax[003,007-011,020-026,049-052,063,105-106,142-144,184-187]
preempt         up 2-00:00:00     65/71/33/169 pax[001-002,004-006,012-019,027-048,053-062,064-075,077-104,107-141,145-183,188-197]
```

### 1.3 Request a Compute Node

```bash
srun --pty -p batch -n 2 --mem=8g --time=0-01:30:00 bash
hostname
```

The hostname should now look like a compute node, not a login node.

**Checkpoint 1.** `hostname` should not start with `login`.

## Part 2: Submit a Batch Job

### 2.1 Write the Job Script

```bash
mkdir -p "$MYWORK/week-01/session-02/logs"
cd "$MYWORK/week-01/session-02"

cat > count_job.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=count
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

echo "Running on: $(hostname)"
grep -vc '^#' "$WEEK1_DATA/anno.gff3"
EOF
```

### 2.2 Spot the Bug Before Submitting

**Question 1.** This job will fail as written. Why?

<details>
<summary>Answer</summary>

`WEEK1_DATA` is set in your current shell, but a batch job starts a fresh shell. Unless you define the variable inside the script, the script does not know where the file is.

This is a common reason a command works interactively but fails in a batch job. The same issue can affect `module load`, `conda activate`, and any other setup step.

</details>

Fix the script by defining `WEEK1_DATA` inside it, then submit:

```bash
sed -i "/^echo \"Running on/a export WEEK1_DATA=$WEEK1_DATA" count_job.sh
cat count_job.sh
sbatch count_job.sh
squeue -u "$USER"
```

The `cat` is not decoration. Read the file and confirm the `export` line landed _before_ the `grep` line that needs it.

### 2.3 Read the Result and the Cost

Replace JOBID below with your own job ID.

```bash
cat logs/count_*.out
seff JOBID
```

`seff` reports what the job actually used against what you asked for. Memory efficiency
will be close to zero here, because counting lines in a text file needs almost nothing.

Over-requesting is not free: SLURM has to find a node with that much memory available, so
your job waits longer in the queue for resources it never touches. Run once with a
generous request, read `seff`, then ask for roughly the observed peak plus a margin.

**Checkpoint 2.** You have a completed job, an output log, and a `seff` report.

## Part 3: Two Ways to Get Software

### 3.1 Conda for Your Own Tools

Conda is available on Pax as a module. Load it first:

```bash
module load miniforge
conda --version
```

Then build an environment for this course and record it:

```bash
conda create -y -n appbio-week1 -c conda-forge -c bioconda seqkit
conda activate appbio-week1
seqkit version
conda env export --from-history > "$MYWORK/week-01/session-02/environment.yml"
cat "$MYWORK/week-01/session-02/environment.yml"
```

Channel order matters. The usual bioinformatics order is `conda-forge` first, then
`bioconda`, and reversing it causes dependency conflicts that are painful to debug.

That `environment.yml` is the deliverable, not the environment itself. It is what lets
somebody else rebuild what you had:

```bash
conda env create -f environment.yml
```

If you want to see what is available in the channels, search for a package:

```bash
conda search -c bioconda seqkit
```

We will not use this environment for the rest of the lab, so deactivate it:

```bash
conda deactivate
```

### 3.2 Nextflow and Containers for Pipelines

You are about to test a real pipeline. Nextflow manages the workflow, and containers provide the pipeline software.

```bash
module load nextflow
nextflow -version

module load singularity
singularity --version
```

Conda and containers both provide software, but they are not equally reproducible. A
conda environment is _solved_ at install time against channels that keep moving, so
rebuilding the same specification in two years can give you different versions. A
container is a fixed filesystem image pinned by tag or digest, identical whenever you
pull it.

Use containers where you can and conda where you must. The pipeline you are about to run
provisions all of its own software as containers, which is why you installed none of the
bioinformatics tools yourself.

## Part 4: Fetch Public Data with `nf-core/fetchngs`

### 4.1 Understand What You Are Asking For

Papers often cite a BioProject accession such as `PRJNA...`. Download tools usually need run accessions such as `SRR...`.

```text
BioProject (PRJNA) -> BioSample (SAMN) -> Experiment (SRX) -> Run (SRR)
```

Set up a clean download directory:

```bash
mkdir -p "$MYWORK/week-01/session-02/fetch"
cd "$MYWORK/week-01/session-02/fetch"
```

Set the project accession from class, then query ENA before downloading anything large:

```bash
export PROJECT=PRJNA792930

ENA_API="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS="run_accession,instrument_platform,library_layout,read_count,base_count"

curl -sS --fail --max-time 120 \
  "${ENA_API}?accession=${PROJECT}&result=read_run&fields=${FIELDS}&format=tsv" \
  -o runs.tsv

[ "$(wc -l < runs.tsv)" -ge 2 ] \
  && echo "OK: $(($(wc -l < runs.tsv) - 1)) runs" \
  || echo "NO DATA: ask for help before continuing"

awk -F'\t' '{ for (i=1; i<=NF; i++) printf "%-22s", ($i=="" ? "-" : $i); print "" }' runs.tsv
```

Expected shape:

```text
run_accession   instrument_platform   library_layout   read_count   base_count
SRR27956204     ILLUMINA              PAIRED           ...
SRR18210286     PACBIO_SMRT           SINGLE           ...
SRR17374240     OXFORD_NANOPORE       SINGLE           ...
```

### 4.2 Choose One Run Per Platform

This project sequenced the same yeast genome on four different platforms. That is unusual
and useful: next week you will compare QC output across all four, and any difference you
see is the platform rather than the sample.

Look at your table again. There are five runs but only four platforms, because the study
includes two Illumina runs. You want one run from each platform.

Write the four accessions into the file fetchngs expects:

```bash
cat > ids.csv <<'EOF'
SRR27956204
SRR18210286
SRR17374239
SRR17374240
EOF

cat ids.csv
```

One run accession per line, no header, no commas. That is the whole format.

Which platform each one is:

```text
SRR27956204   Illumina          paired-end short reads
SRR18210286   PacBio HiFi       long and accurate
SRR17374239   BGISEQ            a second short-read technology
SRR17374240   Oxford Nanopore   long, less accurate
```

**Question 2.** Work out the coverage each of these runs gives, for a 12.1 Mb genome, and
the total you are about to download.

<details>
<summary>Answer</summary>

```bash
grep -f ids.csv runs.tsv | awk -F'\t' '{
  printf "%-13s %-16s %6.0fx  %5.1f Gbp\n", $1, $2, $5/12100000, $5/1e9
  t += $5
} END { printf "\nTOTAL %.1f Gbp, roughly %.0f GB once gzipped\n", t/1e9, t*0.35/1e9 }'
```

`grep -f ids.csv` keeps only the rows whose accession appears in your list, so the totals
describe what you are actually fetching.

Every run is enormously oversequenced for a 12 Mb genome — hundreds of times coverage.
That is normal for deposited data and it is why you will subsample before doing anything
with it. Nobody needs 1,100x to run FastQC, and nobody needs more than about 30x to
assemble HiFi reads.

About 11 GB in total. That is why this is a batch job and not something you watch.

</details>

### 4.3 Test the Pipeline First

```bash
nextflow run nf-core/fetchngs -r 1.13.0 -profile test,singularity --outdir test_out
```

This uses a tiny built-in dataset to confirm that Nextflow runs, containers work, and the cluster can reach the internet.

Always run a test profile before using a new pipeline. Check <https://nf-co.re/fetchngs> for the current release, then pin it with `-r`.

The test has done its job. Clean up before the real run:

```bash
ls -d test_out work .nextflow*
du -sh work
rm -rf test_out work .nextflow.log* .nextflow/
```

`test_out` holds the published results. `work` is where Nextflow actually ran every task,
and it is the one that matters: each task gets its own directory there containing the
staged input files, the resolved script, the logs and the outputs. Even this tiny test
leaves dozens of directories behind.

Three reasons to make this a habit.

**It is the directory that fills your quota.** Published results are a curated copy of a
few files. `work` holds everything, including intermediates you will never look at, and on
a real run it is routinely several times larger than the results. A term of pipeline runs
with no cleanup will exhaust a shared allocation.

**It is reproducible by definition.** Everything in `work` can be regenerated by running
the pipeline again. That is the point of a workflow manager. Deleting it costs you compute
time, not information, which makes it the first thing to remove when you need space.

**Leftover state causes confusing failures.** `.nextflow/` holds the run history and cache.
Starting a new run in a directory still carrying the previous run's cache is a reliable way
to get behaviour you cannot explain.

<details>
<summary>One important exception</summary>

Do **not** delete `work` while a real run is unfinished or might need resuming. `-resume`
reuses completed tasks by looking them up in `work` and `.nextflow/cache`. Delete either
one and the next run starts from zero.

The rule: clean up after a run you are done with, never during one you might continue.
Your fetchngs download in §4.4 should keep its `work` directory until you have confirmed
the FASTQ files are complete in §4.6.

</details>

When you do want to clean up a real run selectively, Nextflow has a command for it:

```bash
nextflow log
nextflow clean -f -before RUN_NAME
```

`nextflow log` lists your previous runs by name; `clean -f -before` removes the work
directories for everything older than the run you name, leaving the most recent intact.

### 4.4 Submit the Real Download as a Job

Do not run the real download interactively. Four runs is roughly 11 GB and will take hours, which is longer than your session and longer than your patience.

```bash
mkdir -p logs

cat > fetch_job.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=fetchngs
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

module load nextflow
module load singularity 2>/dev/null || module load apptainer

export COURSE=/cluster/tufts/bio_appbio
export MYWORK="$COURSE/$USER"
export NXF_SINGULARITY_CACHEDIR="$MYWORK/.singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

nextflow run nf-core/fetchngs \
  -r 1.13.0 \
  -profile singularity \
  --input ids.csv \
  --outdir results \
  --download_method sratools \
  -resume
EOF

cat fetch_job.sh
sbatch fetch_job.sh
squeue -u "$USER"
```

Three details in that script are worth understanding rather than copying.

The heredoc delimiter is quoted, as shown on the `cat` line above. That writes the file exactly as you see it, leaving every variable to be resolved when the job runs.

`COURSE` and `MYWORK` are therefore defined _inside_ the script. This is Question 2 again: the job does not reliably inherit what you set in your own shell, so anything it needs must be set where it runs.

`set -euo pipefail` makes the job stop at the first failure. Without it, a failed `mkdir` would not stop anything, and you would get a confusing container error several minutes later instead of a clear permissions error immediately.

### 4.5 Check After Class

Replace JOBID below with your own job ID.

```bash
cd "$MYWORK/week-01/session-02/fetch"
squeue -u "$USER"
seff JOBID
tail -20 logs/fetchngs_*.out
tail -20 logs/fetchngs_*.err
```

When the job completes:

```bash
ls results/
ls -lh results/fastq/
awk -F, '{ for (i=1; i<=6 && i<=NF; i++) printf "%-20s", $i; print "" }' \
  results/samplesheet/samplesheet.csv | head
```

Expected result directories:

```text
fastq/          the reads, one or two files per run
samplesheet/    samplesheet.csv
metadata/       ENA metadata
pipeline_info/  execution reports, timeline, and software versions
```

You should have six FASTQ files from four runs: the two paired-end runs give `_1` and
`_2`, and the two long-read runs give one file each.

### 4.6 Verify the Download

```bash
cd "$MYWORK/week-01/session-02/fetch/results/fastq"

for f in *.fastq.gz; do
  gzip -t "$f" && echo "OK   $f" || echo "FAIL $f"
done
```

Every file must pass. Now count the reads in each, and compare against what ENA claimed:

```bash
for f in *.fastq.gz; do
  n=$(zcat "$f" | awk 'NR % 4 == 1' | wc -l)
  printf "%-30s %12d reads\n" "$f" "$n"
done

grep -f ../../ids.csv ../../runs.tsv | awk -F'\t' '{printf "%-13s %-16s %12d reads (ENA)\n", $1, $2, $4}'
```

For the paired runs the two mates should have identical counts, and each should match the
ENA figure. For the single-end runs the one file should match directly.

Counting header lines rather than dividing the line count by four is deliberate. A
truncated file gives a line count that is not a multiple of four, so dividing returns a
plausible wrong number with no warning.

**Question 3.** Look at the read counts and the file sizes together. Which platform
produced the fewest reads, and why is that not a sign that anything went wrong?

<details>
<summary>Answer</summary>

The long-read platforms, by a wide margin. PacBio produced a few hundred thousand reads
where Illumina produced tens of millions.

Read _count_ and sequence _yield_ are different things. A nanopore read averaging well
over 10 kb carries as much sequence as a hundred Illumina reads. Comparing platforms on
read count alone is meaningless; compare total bases, which is what you did in Question 4.

This is the distinction next week's QC session is built on.

</details>

**Checkpoint 3.** You have a fetchngs job ID, six FASTQ files, and read counts that match what ENA reported.

> **Do not delete these.** Week 2 runs QC on all four platforms and Week 3 assembles the PacBio run. If you need the space back later, the subsampled copies you make next week are the ones to keep.

## If Your Job Fails

In order:

1. Read the error: `tail -50 logs/fetchngs_*.err`
2. Read the standard output: `tail -50 logs/fetchngs_*.out`
3. Re-run the same job with `-resume` still present.
4. Ask for help with the exact error message.

Two failures are common enough to name.

**Out of disk space.** The `sratools` download method writes large temporary files while converting, often several times the size of the final FASTQ. Check the filesystem with `df -h .`, and your own usage with a `du -sh` on your work directory.

**Permission denied on the container cache.** If `NXF_SINGULARITY_CACHEDIR` points somewhere you cannot write, Singularity fails while pulling images. Confirm the directory exists and is yours:

```bash
ls -ld "$MYWORK/.singularity_cache"
```

If it does not exist and you cannot create it, ask the instructor — your directory under the course path has to be set up for you.

A failed download is normal. Not noticing that it failed is the problem.

## Takeaways

Get a compute node, then submit real work to the scheduler:

```bash
srun --pty -p batch -n 2 --mem=8g --time=0-01:30:00 bash
sbatch script.sh
squeue -u "$USER"
seff JOBID
```

Build an environment and record it so it can be rebuilt:

```bash
module load miniforge
conda create -n env -c conda-forge -c bioconda tool=version
conda env export --from-history > environment.yml
```

Query ENA, then keep the first run of each platform:

```bash
curl -sS --fail "${ENA_API}?accession=${PROJECT}&result=read_run&fields=${FIELDS}&format=tsv" -o runs.tsv
cat > ids.csv <<'EOF'      # one accession per line
SRR27956204
SRR18210286
EOF
```

Test a pipeline before trusting it, then submit the real run as a job:

```bash
nextflow run nf-core/fetchngs -r 1.13.0 -profile test,singularity --outdir test_out
sbatch fetch_job.sh
```

Write a script literally by quoting the heredoc delimiter:

```bash
cat > job.sh <<'EOF'
...
EOF
```

Verify what you downloaded:

```bash
gzip -t file.fastq.gz
zcat file.fastq.gz | awk 'NR % 4 == 1' | wc -l
```

Six habits will save you time this semester:

1. Batch jobs do not inherit your interactive shell setup.
2. Create log directories before submitting jobs.
3. Read the script you just generated before you submit it.
4. Run `-profile test` before using a new workflow.
5. Put long jobs in the scheduler, not in an interactive terminal.
6. Keep everything under `$MYWORK`. Home directory quotas are small and sequencing data is not.
