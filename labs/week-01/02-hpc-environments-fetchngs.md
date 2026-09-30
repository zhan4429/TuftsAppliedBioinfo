# Week 1 Session 2: HPC, Environments, and Fetching Public Data

Applied Bioinformatics - Tufts University Department of Biology

In this lab you will move from the login node to a compute node, submit a small SLURM batch job, create a reproducible software environment, test Nextflow and container support, and start a real public-data download with `nf-core/fetchngs`.

This lab produces PacBio HiFi data used later in the course. The final download may finish after class, so do not skip the job submission and verification steps.

## Before You Start

Open a terminal on Pax:

- OnDemand: <https://ondemand-prod.pax.tufts.edu/> > **Clusters > Tufts HPC Shell Access**
- SSH: `ssh your_utln@login-prod.pax.tufts.edu`

Set paths and make a work directory:

```bash
export COURSE=/cluster/tufts/bio_appbio
export WEEK1_DATA=$COURSE/shared/wk1
mkdir -p ~/appbio/week-01/session-02
cd ~/appbio/week-01/session-02
```

If your instructor gives a different course path, change only `COURSE`.

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
PARTITION    AVAIL  TIMELIMIT   NODES(A/I/O/T)
batch*       up     7-00:00:00  ...
gpu          up     7-00:00:00  ...
interactive  up     4:00:00     ...
largemem     up     7-00:00:00  ...
preempt      up     7-00:00:00  ...
```

**Question 1.** Which partition is appropriate for a 3-hour debugging session? Which is more appropriate for a 5-day assembly?

<details>
<summary>Answer</summary>

Use `interactive` for a 3-hour debugging session. Use `batch` or another instructor-approved long-running partition for a 5-day assembly.

`preempt` can provide more access, but the job can be interrupted by higher-priority work. It is useful for short, resumable jobs and risky for long, non-resumable jobs.

</details>

### 1.3 Request a Compute Node

```bash
srun --pty -p interactive -n 2 --mem=8g --time=0-02:00:00 bash
hostname
```

The hostname should now look like a compute node, not a login node.

**Checkpoint 1.** `hostname` should not start with `login`.

## Part 2: Submit a Batch Job

### 2.1 Write the Job Script

```bash
mkdir -p ~/appbio/week-01/session-02/logs
cd ~/appbio/week-01/session-02

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

**Question 2.** This job will fail as written. Why?

<details>
<summary>Answer</summary>

`WEEK1_DATA` is set in your current shell, but a batch job starts a fresh shell. Unless you define the variable inside the script, the script does not know where the file is.

This is a common reason a command works interactively but fails in a batch job. The same issue can affect `module load`, `conda activate`, and any other setup step.

</details>

Fix the script by defining `WEEK1_DATA` inside it, then submit:

```bash
sed -i "/^echo \"Running on/a export WEEK1_DATA=$WEEK1_DATA" count_job.sh
sbatch count_job.sh
squeue -u "$USER"
```

### 2.3 Read the Result and the Cost

Replace `<jobid>` with your job ID.

```bash
cat logs/count_*.out
seff <jobid>
```

**Question 3.** Memory efficiency will be near zero. Why can over-requesting memory still cost you something?

<details>
<summary>Answer</summary>

SLURM must find a node with the requested memory free. If you ask for much more than you need, your job may wait longer and block resources that someone else could use.

Good habit: run once with a reasonable request, check `seff`, then request roughly the observed peak plus a buffer.

</details>

**Checkpoint 2.** You have a completed job, an output log, and a `seff` report.

## Part 3: Two Ways to Get Software

### 3.1 Conda or Mamba for Your Own Tools

```bash
module load anaconda
conda create -y -n appbio-week1 -c conda-forge -c bioconda seqkit=2.14.0
conda activate appbio-week1
seqkit version
conda env export --from-history > ~/appbio/week-01/session-02/environment.yml
```

Channel order matters. The usual bioinformatics order is `conda-forge` first, then `bioconda`.

If `mamba` is available on Pax, the create step can be faster:

```bash
mamba create -y -n appbio-week1 -c conda-forge -c bioconda seqkit=2.14.0
```

### 3.2 Nextflow and Containers for Pipelines

You are about to test a real pipeline. Nextflow manages the workflow, and containers provide the pipeline software.

```bash
module load nextflow
nextflow -version

module load singularity 2>/dev/null || module load apptainer
singularity --version 2>/dev/null || apptainer --version
```

**Question 4.** Conda environments and containers both provide software. Why is a container often more reproducible?

<details>
<summary>Answer</summary>

A conda environment is solved at install time against channels that change. Re-solving the same specification later can produce different versions.

A container is a fixed filesystem image pinned by tag or digest. Use containers where you can and conda where you must.

</details>

## Part 4: Fetch Public Data with `nf-core/fetchngs`

### 4.1 Understand What You Are Asking For

Papers often cite a BioProject accession such as `PRJNA...`. Download tools usually need run accessions such as `SRR...`.

```text
BioProject (PRJNA) -> BioSample (SAMN) -> Experiment (SRX) -> Run (SRR)
```

Set up a clean download directory:

```bash
mkdir -p ~/appbio/week-01/session-02/fetch
cd ~/appbio/week-01/session-02/fetch
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

### 4.2 Extract the Run You Need

This is Session 1's `awk` pattern applied to a real task:

```bash
awk -F'\t' 'NR>1 && $2=="PACBIO_SMRT" {print $1}' runs.tsv > ids.csv
cat ids.csv
```

Expected:

```text
SRR18210286
```

`ids.csv` should contain one run accession per line, with no header and no commas.

**Question 5.** Compute the coverage for this run assuming a 12.1 Mb genome. Is it more than you need?

<details>
<summary>Answer</summary>

```bash
awk -F'\t' '$2=="PACBIO_SMRT" {printf "%.0fx\n", $5/12100000}' runs.tsv
```

The result is hundreds of times coverage. HiFi assembly is usually comfortable at roughly 20-30x because the reads are long and accurate, so this dataset will be subsampled before assembly.

</details>

### 4.3 Test the Pipeline First

```bash
nextflow run nf-core/fetchngs -r 1.13.0 -profile test,singularity --outdir test_out
```

This uses a tiny built-in dataset to confirm that Nextflow runs, containers work, and the cluster can reach the internet.

Always run a test profile before using a new pipeline. Check <https://nf-co.re/fetchngs> for the current release, then pin it with `-r`.

### 4.4 Submit the Real Download as a Job

Do not run the real download interactively. It is large enough to take longer than your session.

```bash
mkdir -p logs

cat > fetch_job.sh <<EOF
#!/bin/bash
#SBATCH --job-name=fetchngs
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load nextflow
module load singularity 2>/dev/null || module load apptainer

export COURSE="$COURSE"
export NXF_SINGULARITY_CACHEDIR="$COURSE/$USER/.singularity_cache"
mkdir -p "\$NXF_SINGULARITY_CACHEDIR"

nextflow run nf-core/fetchngs \
  -r 1.13.0 \
  -profile singularity \
  --input ids.csv \
  --outdir results \
  --download_method sratools \
  -resume
EOF

sbatch fetch_job.sh
squeue -u "$USER"
```

**Question 6.** Why ask for twelve hours of wall time for a download, and why include `-resume`?

<details>
<summary>Answer</summary>

Public archives can be unpredictable, especially if many students start downloads at once. Generous wall time ends early if the job finishes early; wall time that is too short kills a job near the end.

`-resume` lets a failed or interrupted Nextflow run reuse completed work instead of starting over.

</details>

### 4.5 Check After Class

Replace `<jobid>` with your job ID.

```bash
cd ~/appbio/week-01/session-02/fetch
squeue -u "$USER"
seff <jobid>
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
fastq/          HiFi reads used later in the course
samplesheet/    samplesheet.csv
metadata/       ENA metadata
pipeline_info/  execution reports, timeline, and software versions
```

### 4.6 Verify the Download

```bash
cd ~/appbio/week-01/session-02/fetch/results/fastq
for f in *.fastq.gz; do
  gzip -t "$f" && echo "OK $f" || echo "FAIL $f"
done

zcat *.fastq.gz | awk 'NR % 4 == 1' | wc -l
```

Compare the read count with the `read_count` column in `../../runs.tsv`.

**Checkpoint 3.** You have a fetchngs job ID, and you know how to check whether the download succeeded.

## If Your Job Fails

In order:

1. Read the error: `tail -50 logs/fetchngs_*.err`
2. Read the standard output: `tail -50 logs/fetchngs_*.out`
3. Re-run the same job with `-resume` still present.
4. Ask for help with the exact error message.

A failed download is normal. Not noticing that it failed is the problem.

## Takeaways

| Task | Command |
| --- | --- |
| Request a compute node | `srun --pty -p interactive -n 2 --mem=8g --time=0-02:00:00 bash` |
| Submit and watch | `sbatch script.sh`, `squeue -u "$USER"` |
| Check resource use | `seff <jobid>` |
| Create an environment | `conda create -n env -c conda-forge -c bioconda tool=version` |
| Record an environment | `conda env export --from-history > environment.yml` |
| Query ENA runs | ENA portal API with `result=read_run` |
| Filter to one platform | `awk -F'\t' 'NR>1 && $2=="PACBIO_SMRT" {print $1}' runs.tsv` |
| Test a pipeline | `nextflow run nf-core/fetchngs -r 1.13.0 -profile test,singularity` |
| Fetch data properly | Submit `nf-core/fetchngs` as a batch job |
| Verify FASTQ files | `gzip -t`, then count FASTQ records |

Four habits will save you time this semester:

1. Batch jobs do not inherit your interactive shell setup.
2. Create log directories before submitting jobs.
3. Run `-profile test` before using a new workflow.
4. Put long jobs in the scheduler, not in an interactive terminal.
