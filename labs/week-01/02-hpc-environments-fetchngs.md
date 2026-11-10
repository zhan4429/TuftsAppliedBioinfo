# Week 1 Session 2: HPC, Environments, and Fetching Public Data

Applied Bioinformatics - Tufts University Department of Biology

In this lab you will request a compute node, load software modules, and submit **one
fetchngs batch job** to download public sequencing data. You do not need to wait for the
download to finish during class.


## Before You Start

Open a terminal on Pax:

- OnDemand: <https://ondemand-prod.pax.tufts.edu/> > **Clusters > Tufts HPC Shell Access**
- SSH: `ssh your_utln@login-prod.pax.tufts.edu`

Set the course paths and make your work directory:

```bash
export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"

mkdir -p "$MYWORK/week-01/session-02/fetch/logs"
cd "$MYWORK/week-01/session-02/fetch"
pwd
```

`$COURSE` is the shared course allocation. `$MYWORK` is your own directory, named after
your UTLN. Keep the downloaded reads and results there. Ask for help if you cannot create
your directory. If you reconnect later, repeat these path definitions.

## Part 1: Request a Compute Node

```bash
hostname
sinfo -s
```

A hostname beginning with `login` identifies a shared login node. `sinfo -s` shows the
cluster partitions and their availability. We will use the CPU partition `batch`.

Request an interactive shell on a compute node:

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=2 \
  --mem=8G --time=01:30:00 --pty bash

hostname
echo "Interactive job ID: $SLURM_JOB_ID"
```

The request may wait in the queue. When it starts, the hostname should identify a compute
node. `--ntasks=1` requests one task, `--cpus-per-task=2` gives it two CPUs, and `--mem=8G`
requests 8 GB of memory for the allocation.

**Stay in this shell for Parts 2 and 3.** We will use it to check software and prepare and
submit the batch script. The download will run in its own batch allocation with the
resources requested in that script.

**Checkpoint 1.** You are on a compute node. There is no need to request another
interactive session during class.

## Part 2: Load the Software

```bash
module load nextflow singularity
nextflow -version
singularity --version
module list
```

`module load` makes installed software available in your current shell. `module list`
shows the modules you have loaded. Nextflow coordinates the steps of the pipeline;
Singularity provides the containers containing the tools used by those steps.

This lab uses **nf-core/fetchngs 1.13.0**, which requires **Nextflow 25.10.4 or later**.


## Part 3: Submit the Fetchngs Download

### 3.1 Choose the Data

We will download these four runs from BioProject `PRJNA792930`:

| Run | Platform | Library |
| --- | --- | --- |
| SRR27956204 | Illumina | Paired-end RNA-seq |
| SRR18210286 | PacBio HiFi | Single-end genomic sequencing |
| SRR17374239 | BGISEQ | Paired-end genomic sequencing |
| SRR17374240 | Oxford Nanopore | Single-end genomic sequencing |

These data support Week 2 QC; the PacBio HiFi run also supports Week 3 assembly.
**The Illumina run is RNA-seq, while the other three are genomic sequencing.** Their QC
reports can differ because of library type as well as platform.

Create the input file:

```bash
cd "$MYWORK/week-01/session-02/fetch"

cat > ids.csv <<'EOF'
SRR27956204
SRR18210286
SRR17374239
SRR17374240
EOF

cat ids.csv
```

Despite its `.csv` extension, this input file has **one run accession per line, with no
header or commas**. Using explicit run IDs controls which data we download.

ENA advertised about 22 GB of compressed FASTQ for these runs when checked on
2026-11-03. SRA conversion and intermediate files require additional space. 

### 3.2 Write the Job Script

Write the batch script below. The pipeline will use its default configuration, and its
tasks will run inside the batch allocation. We will not run the pipeline directly in
our smaller interactive allocation.

```bash
mkdir -p logs

cat > fetch_job.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=fetchngs
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load nextflow singularity

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
cd "$MYWORK/week-01/session-02/fetch"

export NXF_SINGULARITY_CACHEDIR="$COURSE/$USER/.singularity_cache"
export SINGULARITY_CACHEDIR="$MYWORK/.singularity/oci"
mkdir -p "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR"

nextflow run nf-core/fetchngs \
  -r 1.13.0 \
  -profile singularity \
  --input ids.csv \
  --outdir results \
  --download_method sratools \
  -resume
EOF
```

Read the script from top to bottom:

- `#SBATCH` lines request resources and name the log files. `%x` becomes the job name;
  `%j` becomes the job ID. Create `logs/` before submitting.
- `module load` loads software inside the batch job, and `cd` selects the work directory.
- `NXF_SINGULARITY_CACHEDIR` keeps downloaded container images in your own course
  directory so Nextflow can reuse them. Nextflow otherwise uses its default home and
  Java settings.
- `-r` selects the pipeline version and `-profile singularity` selects containers.
- `--input` names the accession list and `--outdir` names the results directory.
  `--download_method sratools` downloads and converts the SRA data to FASTQ. We use this
  method because ENA lists only one FASTQ URL for the paired BGISEQ run.
- `-resume` allows completed pipeline steps to be reused if the job needs restarting.

The quoted `<<'EOF'` writes the script literally, so variables inside it are evaluated
when the batch job runs. If your course directory differs, change the `COURSE` line in
both the setup and the script.

### 3.3 Submit Once and Copy the Job ID

From your current compute-node shell:

```bash
bash -n fetch_job.sh
cat fetch_job.sh
sbatch fetch_job.sh
squeue --me
```

`bash -n` checks the script's shell syntax; `cat` lets you review it. `sbatch` submits the
job and prints a message such as:

```text
Submitted batch job 123456
```

Your job ID will be different. `squeue --me` lists your pending and running jobs. The
equivalent command is:

```bash
squeue -u "$USER"
```

Find the row whose `NAME` is **`fetchngs`** and copy its `JOBID` into your notes.
You may also see your interactive shell job; that is a different allocation.
`PD` means pending and `R` means running. A pending job has already been submitted;
do not submit another copy just because it has not started yet.

**In-class checkpoint.** You have `ids.csv`, `fetch_job.sh`, and the job ID
for your fetchngs batch job. You are done with the in-class work.

### 3.4 Leave the Interactive Session

```bash
exit
```

This releases the interactive allocation and returns you to the login node. The fetchngs
batch job has its own allocation and can run after you leave the interactive shell or
log out of Pax.

## After Class: Check the Download

### Check the Job and Its Logs

Reconnect to Pax, repeat the path definitions in **Before You Start**, and return to the
download directory. **Replace `123456` below with your fetchngs job ID.** These commands
inspect small logs and can run on the login node:

```bash
cd "$MYWORK/week-01/session-02/fetch"
squeue --me
sacct -j 123456 --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
tail -n 30 logs/fetchngs_123456.out
tail -n 30 logs/fetchngs_123456.err
```

A job disappearing from `squeue` does not prove success. Check `sacct` for `COMPLETED`
and exit code `0:0`; accounting information can take a short time to appear.
Once it finishes successfully:

```bash
seff 123456
ls -lh results/fastq/
head -n 5 results/samplesheet/samplesheet.csv
```

`seff` reports resource use. `results/fastq/` contains the reads, `samplesheet/` connects
files to samples, and `metadata/` and `pipeline_info/` contain metadata and run reports.
Expect six FASTQ files: two mates each for Illumina and BGISEQ, and one file each for
PacBio and Nanopore. Use the samplesheet to identify them; filenames may use experiment
identifiers rather than the run IDs in `ids.csv`.

### Verify the FASTQ Files on a Compute Node

Do this after the download has completed. Request a fresh interactive session for these
large-file checks:

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=6 \
  --mem=18G --time=01:30:00 --pty bash
cd "$MYWORK/week-01/session-02/fetch"
module load seqkit

ls results/fastq/*.fastq.gz | wc -l
gzip -t results/fastq/*.fastq.gz
seqkit stats -T -j 6 results/fastq/*.fastq.gz > fastq_stats.tsv
cat fastq_stats.tsv
```

The file count should be six. `gzip -t` prints nothing when all files pass; stop and ask
for help if it reports an error. SeqKit should also finish without a parsing error.
`num_seqs` counts reads and `sum_len` counts bases. The two mate files of each paired-end
library should have equal read counts. Long-read files can have fewer reads while still
containing many bases because each read is longer.

When you have reviewed the statistics, release the compute node:

```bash
exit
```

Keep the reads and `fastq_stats.tsv` for Week 2. The actual input path is
`$MYWORK/week-01/session-02/fetch/results/fastq/`; use that path when a later handout asks
where your data are.

## If Your Job Fails

- **Pending:** read the reason in `squeue --me` and wait or ask for help. Do not submit
  duplicate copies.
- **Failed or timed out:** read the job's `.err` and `.out` logs and `.nextflow.log` in
  the download directory. Give the instructor the job ID and exact error message.
- **Disk or cache errors:** check the allocation quota and your container-cache directory with
  the instructor. Conversion needs more space than the final compressed FASTQ files.

After fixing the problem, resubmit `fetch_job.sh` from the same download directory with
`sbatch fetch_job.sh`, and record the new job ID. Keep `work/` and `.nextflow/` so that
`-resume` can reuse completed steps. Keep these directories until your results are
verified; cleanup is not part of this class session.

## What to Keep

- `ids.csv` and `fetch_job.sh`.
- Your fetchngs job ID and its `.out` and `.err` logs.
- The downloaded results and, after verification, `fastq_stats.tsv`.

Keep scripts and notes in your coursework repository. Keep large reads, caches, and
results on cluster storage.

## References

- [Tufts HPC: running software with SLURM](https://rtguides.it.tufts.edu/hpc/access/60-brief-slurm.html)
- [nf-core/fetchngs usage](https://nf-co.re/fetchngs/usage)
- [nf-core/fetchngs 1.13.0 runtime requirement](https://github.com/nf-core/fetchngs/blob/1.13.0/nextflow.config)
- [SeqKit usage and FASTQ statistics](https://bioinf.shenwei.me/seqkit/usage/)
