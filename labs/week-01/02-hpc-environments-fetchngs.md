# Week 1 Session 2: HPC, Environments, and Fetching Public Data

Applied Bioinformatics - Tufts University Department of Biology

In this lab you will request a compute node, submit a SLURM batch job, record a software environment, and use Nextflow with containers to download public sequencing data. You will leave with a job ID and a way to verify the results after class.

These data support the four-platform QC exercise in Week 2; the PacBio HiFi run also supports genome assembly in Week 3. **The selected Illumina run is RNA-seq, while the other three runs are genomic sequencing.** You can compare file formats and QC reports, but differences cannot be attributed to sequencing platform alone.

## Before You Start

Open a terminal on Pax:

- OnDemand: <https://ondemand-prod.pax.tufts.edu/> > **Clusters > Tufts HPC Shell Access**
- SSH: `ssh your_utln@login-prod.pax.tufts.edu`

Use the course allocation announced in class. The shared annotation path comes from Session 1 and may be in a different allocation from your personal work:

```bash
export COURSE=/cluster/tufts/c27sp1bio0196
export WEEK1_DATA=$COURSE/shared/wk1
export MYWORK="$COURSE/$USER"

mkdir -p "$MYWORK/week-01/session-02/logs"
cd "$MYWORK/week-01/session-02"
pwd
ls -lh "$WEEK1_DATA/anno.gff3"
```

If either directory is unavailable, ask for the correct path before continuing. `$COURSE` is shared; `$MYWORK` is your own directory, named after your UTLN. Keep large data, environments, and pipeline outputs under `$MYWORK`.

If you reconnect later, repeat these path definitions before using commands that contain `$MYWORK` or `$WEEK1_DATA`.

## Part 1: Get Off the Login Node

### 1.1 Inspect the Cluster

```bash
hostname
sinfo -s
```

A hostname beginning with `login` identifies a shared login node. Use it to edit files, submit jobs, and inspect small logs. Install environments, run pipelines, and scan large FASTQ files on compute nodes. `sinfo -s` lists partitions and their current availability. We use the CPU partition `batch`; node counts and time limits can change.

### 1.2 Request One Task with Multiple CPUs

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=2 \
  --mem=8G --time=01:30:00 --pty bash

hostname
printf 'Job: %s; CPUs: %s\n' "$SLURM_JOB_ID" "$SLURM_CPUS_PER_TASK"
```

The request may wait in the queue. Once it starts, the hostname should identify a compute node.

`--ntasks=1` requests one task; `--cpus-per-task=2` gives it two CPUs. Using `-n 2` would request two tasks, which is not what we want for one shell or one multithreaded program. `--mem=8G` requests memory for the node allocation, not 8 GB per CPU.

**Checkpoint 1.** You are on a compute node and can identify your job ID and CPU allocation.

## Part 2: Submit a Batch Job

### 2.1 Write a Script with a Clear Input Path

```bash
cd "$MYWORK/week-01/session-02"
mkdir -p logs

cat > count_job.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=count
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

WEEK1_DATA=/cluster/tufts/c27sp1bio0196/shared/wk1
echo "Running on: $(hostname)"
grep -vc '^#' "$WEEK1_DATA/anno.gff3"
EOF

bash -n count_job.sh
cat count_job.sh
COUNT_JOB_ID=$(sbatch --parsable count_job.sh)
printf 'Count job: %s\n' "$COUNT_JOB_ID"
squeue -j "$COUNT_JOB_ID"
```

The quoted heredoc delimiter, `<<'EOF'`, writes the script literally: `$WEEK1_DATA` and `$(hostname)` are evaluated when the job runs. `bash -n` checks syntax without running commands. `sbatch --parsable` returns the job ID, not proof that the job has finished.

SLURM reads `#SBATCH` directives before the shell runs. `$MYWORK` inside a directive would be literal text, so use relative log paths or pass expanded paths as `sbatch` command-line options. Create `logs/` before submission: SLURM opens those files before executing the script.

### 2.2 Understand What the Job Inherits

**Question 1.** Does a batch job inherit `WEEK1_DATA` from your current shell?

<details>
<summary>Answer</summary>

By default, `sbatch` exports the submitting process's environment, including variables you set with `export`. The `export WEEK1_DATA=...` in the setup block would normally be available in the job. An ordinary, unexported shell variable would not.

A new job shell does not reproduce every interactive setup action or shell function. Load the modules needed by each job inside its script. Here, defining `WEEK1_DATA` inside the script makes the input path easy to see and edit.

</details>

### 2.3 Confirm Completion and Read the Logs

```bash
squeue -j "$COUNT_JOB_ID"
sacct -j "$COUNT_JOB_ID" --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

A job disappearing from `squeue` may have completed, failed, or been cancelled. Check `sacct` for `COMPLETED` and exit code `0:0`; accounting information can take a short time to appear. Once the job finishes:

```bash
cat "logs/count_${COUNT_JOB_ID}.out"
cat "logs/count_${COUNT_JOB_ID}.err"
seff "$COUNT_JOB_ID"
```

The output should contain a compute-node hostname and the number of non-comment lines in `anno.gff3`. `seff` compares measured use with requested resources. Very short jobs may have incomplete memory measurements, so do not size future jobs from this tiny example alone. For real workloads, measure representative inputs and leave a margin.

**Checkpoint 2.** You have a completed job, its exact output and error logs, and a resource-use report.

## Part 3: Two Ways to Get Software

Request a new interactive compute session to check the software modules and test the pipeline:

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=6 \
  --mem=18G --time=01:30:00 --pty bash
cd "$MYWORK/week-01/session-02"
```

### 3.1 Modules for Individual Tools

Pax already provides SeqKit as a module. Load it and check its version:

```bash
module avail seqkit
module load seqkit
seqkit version
module list
```

`module avail` shows available software, `module load` makes a tool available in your current shell, and `module list` shows what you have loaded. Record the module names and versions in your notes so you can use the same software later. You do not need to install SeqKit or activate a Conda environment for this lab.

### 3.2 Nextflow and Containers for Pipelines

```bash
module load nextflow singularity
nextflow -version
java -version
singularity --version

export NXF_HOME="$SCRATCH/.nextflow"
export NXF_SINGULARITY_CACHEDIR="$SCRATCH/.singularity_cache"
export SINGULARITY_CACHEDIR="$MYWORK/.singularity/oci"
export NXF_OPTS='-Xms512m -Xmx4g'
mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR"
```

Pax defines `$SCRATCH` for your scratch directory. These settings put Nextflow's downloaded files and cached container images there.

This lab pins **nf-core/fetchngs 1.13.0**, which requires **Nextflow 25.10.4 or later**. If the loaded module is older, ask for a compatible module before continuing. Keep the same pipeline version for the test and the real download.

Nextflow schedules workflow tasks; Singularity supplies their software. A container tag can be changed by its publisher, so a tag alone does not guarantee an immutable image. Record the pipeline release, runtime versions, parameters, and container identifiers; image digests provide stronger identity.

### 3.3 Keep Tasks Inside the Allocation

We will run Nextflow and its tasks inside one SLURM allocation using Nextflow's `local` executor. Write a configuration shared by the test and the real run:

```bash
cat > lab.config <<'EOF'
process {
    executor = 'local'
    resourceLimits = [cpus: 6, memory: 12.GB, time: 23.h]
}
executor {
    cpus = 6
    memory = 12.GB
    queueSize = 1
}
EOF
```

The tasks share six CPUs and a 12 GB memory budget; the 18 GB SLURM allocation leaves room for Nextflow's JVM and other overhead. `queueSize = 1` runs one workflow task at a time. The outer SLURM time limit still applies, including the shorter limit of the interactive test.

Without these limits, the pipeline's default requests can exceed this allocation. A per-task limit alone also does not bound the total resources of several concurrent tasks.

Use `-profile singularity` for this lab. The separate `tufts` profile uses the SLURM executor to submit workflow tasks as additional jobs; that is a different execution model.

## Part 4: Fetch Public Data with nf-core/fetchngs

### 4.1 Inspect Metadata Before Downloading

The pipeline accepts study, sample, experiment, and run accessions. A project accession can expand to many runs, so we will inspect the project and select explicit runs:

```text
BioProject (PRJNA) -> BioSample (SAMN) -> Experiment (SRX) -> Run (SRR)
```

```bash
mkdir -p "$MYWORK/week-01/session-02/fetch"
cd "$MYWORK/week-01/session-02/fetch"

export PROJECT=PRJNA792930
ENA_API='https://www.ebi.ac.uk/ena/portal/api/filereport'
FIELDS='run_accession,experiment_accession,instrument_platform,library_layout,library_strategy,library_source,read_count,base_count,fastq_bytes,experiment_title'

curl --fail --silent --show-error --location --retry 3 --max-time 120 \
  --get "$ENA_API" \
  --data-urlencode "accession=$PROJECT" \
  --data-urlencode 'result=read_run' \
  --data-urlencode "fields=$FIELDS" \
  --data-urlencode 'format=tsv' \
  -o runs.tsv

awk -F'\t' 'NR == 1 { if ($1 != "run_accession") exit 1 }
  END { if (NR < 2) exit 1 }' runs.tsv
cut -f1,3-6 runs.tsv
```

Stop if the request or validation fails; a header-only file is not a valid result. Save `runs.tsv` with your download date. The number of runs in a public project can change.

### 4.2 Select Four Runs and Check Their Meaning

```bash
cat > ids.csv <<'EOF'
SRR27956204
SRR18210286
SRR17374239
SRR17374240
EOF

awk -F'\t' 'NR == FNR { wanted[$1] = 1; next }
  FNR == 1 || ($1 in wanted)' ids.csv runs.tsv > selected_runs.tsv

awk -F'\t' 'NR > 1 { platforms[$3] = 1 }
  END {
    for (p in platforms) n++
    if (NR != 5 || n != 4) exit 1
    print "OK: four selected runs from four platforms"
  }' selected_runs.tsv
```

`ids.csv` contains one accession per line, with no header or commas. The `awk` join matches the complete accession field, rather than matching arbitrary text in a row. Stop if all four selected runs are not present.

The ENA records checked on 2026-11-03 describe:

| Run         | Platform                                  | Library                                 | Use in this course                            |
| ----------- | ----------------------------------------- | --------------------------------------- | --------------------------------------------- |
| SRR27956204 | Illumina NovaSeq 6000                     | Paired-end RNA-seq, no-stress condition | Short-read QC; not a genomic coverage control |
| SRR18210286 | PacBio Sequel, experiment titled HiFi-Seq | Single-end genomic sequencing           | Long-read QC and HiFi assembly                |
| SRR17374239 | BGISEQ-500                                | Paired-end genomic sequencing           | Short-read genomic QC                         |
| SRR17374240 | Oxford Nanopore PromethION                | Single-end genomic sequencing           | Long-read QC                                  |

**Question 2.** Why does a shared BioProject not make this a controlled comparison of sequencing platforms?

<details>
<summary>Answer</summary>

A project groups related submissions; it does not guarantee the same library preparation or experimental conditions. Here, the Illumina data measure RNA abundance, while the other libraries sample genomic DNA. Read length, quality, duplication, and sequence composition can reflect library type, preparation, and processing as well as platform. Inspect `library_strategy`, `library_source`, and experiment descriptions before interpreting differences. Both Illumina runs currently in this project are RNA-seq.

</details>

The PacBio experiment and sample are explicitly labelled HiFi in the archive, even though its ENA FASTQ filename contains `subreads`. Check the experiment metadata and publication rather than classifying the data from the filename alone.

### 4.3 Estimate Yield and Storage

**Question 3.** How many bases and how much compressed data are advertised for the selected runs?

```bash
awk -F'\t' 'NR > 1 {
  bytes = 0
  n = split($9, sizes, ";")
  for (i = 1; i <= n; i++) bytes += sizes[i]
  printf "%-13s %-16s %6.2f Gbp  %6.2f GB advertised FASTQ\n", $1, $3, $8/1e9, bytes/1e9
  total_bases += $8
  total_bytes += bytes
} END {
  printf "TOTAL: %.2f Gbp; %.2f GB advertised FASTQ\n", total_bases/1e9, total_bytes/1e9
}' selected_runs.tsv

df -h .
du -sh "$MYWORK"
```

At the metadata check above, ENA advertised about **31.40 Gbp and 21.67 GB** of compressed FASTQ, not 11 GB. The byte estimate sums every semicolon-separated file size, including both mates where listed. The SRA conversion used below can produce different files and compression sizes, so this is an estimate, not a quota requirement.

SRA downloads, uncompressed FASTQ, temporary conversion files, published copies, and containers require additional space. `df` reports filesystem free space, not necessarily your personal quota; check the course allocation's quota before starting. Ask about staged or shared copies if the class does not have room for individual full downloads.

For a genomic library, nominal coverage is total sequenced bases divided by the assumed genome size, about 12.1 Mb here. **Do not interpret RNA-seq yield divided by genome size as uniform genomic coverage.** Archive counts can also use different conventions; use parsed FASTQ statistics to check the actual yield after conversion. There is no universal compression ratio or assembly coverage target.

### 4.4 Test the Pipeline in Its Own Directory

Still on the interactive compute node:

```bash
mkdir -p "$MYWORK/week-01/session-02/smoke-test"
cd "$MYWORK/week-01/session-02/smoke-test"

nextflow run nf-core/fetchngs \
  -r 1.13.0 \
  -profile test,singularity \
  -c ../lab.config \
  --download_method sratools \
  --outdir results
```

The maintainer-provided test inputs check pipeline startup, containers, and data access. Require successful completion before submitting the real download. This test uses the same release, container engine, resource configuration, and download method as the real run.

Keep the test in `smoke-test/` and the real run in `fetch/`. Each directory has its own `.nextflow/` cache, logs, and `work/` tasks. Previous run state is normal and enables `-resume`; it is not inherently a source of failures.

After the test succeeds, release the interactive allocation:

```bash
exit
hostname
```

### 4.5 Submit the Real Download from the Login Node

Return to the download directory. The input and metadata files you created earlier should still be there:

```bash
cd "$MYWORK/week-01/session-02/fetch"
ls -l ids.csv selected_runs.tsv ../lab.config
mkdir -p logs

cat > fetch_job.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=fetchngs
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=18G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

module load nextflow singularity

export COURSE="/cluster/tufts/c27sp1bio0196"
export MYWORK="$COURSE/$USER"
cd "$MYWORK/week-01/session-02/fetch"

export NXF_HOME="$SCRATCH/.nextflow"
export NXF_SINGULARITY_CACHEDIR="$SCRATCH/.singularity_cache"
export SINGULARITY_CACHEDIR="$MYWORK/.singularity/oci"
export NXF_OPTS='-Xms512m -Xmx4g'
mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR"

nextflow run nf-core/fetchngs \
  -r 1.13.0 \
  -profile singularity \
  -c ../lab.config \
  --input ids.csv \
  --outdir results \
  --download_method sratools \
  -resume

EOF

bash -n fetch_job.sh
cat fetch_job.sh
FETCH_JOB_ID=$(sbatch --parsable fetch_job.sh)
printf '%s\n' "$FETCH_JOB_ID" > fetch_job_id.txt
squeue -j "$FETCH_JOB_ID"
```

We use `sratools` deliberately. ENA currently lists the paired BGISEQ run with a single FASTQ URL; downloading the SRA data and splitting biological reads avoids assuming that ENA's file listing already contains two mate files. The pipeline can also resolve run IDs back to their experiment and merge runs belonging to that experiment. Inspect the generated mappings and samplesheet to confirm what was actually fetched.

Read the script from top to bottom: request resources, load the software, set the work directory and cache locations, then run Nextflow. If your course directory differs, change the `COURSE` line to match **Before You Start**. The cache settings match the test run. The download is the final command, so its exit status determines whether SLURM reports the job as successful. We will check the downloaded files separately below.

### 4.6 Check After Class

Reconnect to Pax and repeat the path definitions in **Before You Start**. These commands inspect small files and can run on the login node:

```bash
cd "$MYWORK/week-01/session-02/fetch"
FETCH_JOB_ID=$(cat fetch_job_id.txt)

squeue -j "$FETCH_JOB_ID"
sacct -j "$FETCH_JOB_ID" --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
tail -n 30 "logs/fetchngs_${FETCH_JOB_ID}.out"
tail -n 30 "logs/fetchngs_${FETCH_JOB_ID}.err"
```

Once the job finishes, require `COMPLETED` with exit code `0:0`, then inspect:

```bash
seff "$FETCH_JOB_ID"
ls -lh results/fastq/
head -n 5 results/samplesheet/samplesheet.csv
```

Expected result directories:

```text
fastq/          published reads
samplesheet/    samplesheet.csv and identifier mappings
metadata/       archive metadata used by fetchngs
pipeline_info/  execution reports, parameters, and software versions
```

For these four experiments, expect six FASTQ files: two mates each for Illumina and BGISEQ, and one file each for PacBio and Nanopore. File names may use experiment identifiers or other pipeline-generated prefixes; use the samplesheet and mappings to connect them to your selected run IDs.

Now request a compute node to check the large files. Keep this separate from the download script so you can inspect each step:

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

The file count should be six. `gzip -t` prints nothing when all files pass; stop and ask for help if it reports an error. SeqKit should also finish without a parsing error. `num_seqs` counts records and `sum_len` counts sequence bases. The two files from each paired-end library should have equal record counts. Compare those figures with `selected_runs.tsv` and fetchngs metadata, checking whether an archive count means reads or read pairs before assuming equality.

When you have reviewed the statistics, release the compute node:

```bash
exit
```

`gzip -t` detects damaged compressed data; it does not establish that all expected reads arrived. Likewise, counting every fourth line is not FASTQ validation: it can return a plausible count for a truncated record. Review parser success, expected file layout, mate counts, and total yield together. For direct archive downloads, also verify the archive's checksums; an ENA checksum will not necessarily match a FASTQ regenerated and recompressed from SRA.

**Question 4.** Why can a long-read file contain fewer records but a similar number of bases to a short-read dataset?

<details>
<summary>Answer</summary>

Each long read contains more sequence. Compare both `num_seqs` and `sum_len`, not read count alone. For paired-end yield, add the bases in both mate files. Remember that RNA-seq and genomic sequencing sample different molecules, so the Illumina yield has a different biological interpretation.

</details>

**Checkpoint 3.** You have the job ID and successful job state, six gzip-tested and parsed FASTQ files, the samplesheet, and a statistics table with plausible mate counts and sequence yields.

> Keep the FASTQ files for Week 2 and the HiFi assembly exercise in Week 3. Give later labs the actual `results/fastq/` path from this run; do not assume a different handout's example path points here.

### 4.7 Resume First; Clean Up After Verification

If a job fails, diagnose and fix the cause before resubmitting. Return to the same `fetch/` directory and use the same inputs, release, configuration, and `-resume` option:

```bash
cd "$MYWORK/week-01/session-02/fetch"
FETCH_JOB_ID=$(sbatch --parsable fetch_job.sh)
printf '%s\n' "$FETCH_JOB_ID" > fetch_job_id.txt
```

Keep both `work/` and `.nextflow/cache/` while you may need to resume. They contain reusable task outputs and the records that locate them. Published results, metadata, and logs are also worth keeping; intermediate files are not automatically disposable just because a workflow produced them.

When the test run is finished and no longer needed, inspect a cleanup dry run from its own directory:

```bash
cd "$MYWORK/week-01/session-02/smoke-test"
nextflow log
nextflow clean RUN_NAME -n -keep-logs
```

Replace `RUN_NAME` with the exact test run name. Only after reviewing the dry run, `nextflow clean RUN_NAME -f -keep-logs` removes that run's temporary task files while retaining execution records and logs. Do not apply cleanup to an active run or to a real download that has not been verified. Retain its reports for troubleshooting.

## If Your Job Fails

Check the job state and exact logs before retrying:

```bash
cd "$MYWORK/week-01/session-02/fetch"
FETCH_JOB_ID=$(cat fetch_job_id.txt)
sacct -j "$FETCH_JOB_ID" --format=JobID,State,ExitCode,Elapsed,MaxRSS
tail -n 50 "logs/fetchngs_${FETCH_JOB_ID}.err"
tail -n 50 "logs/fetchngs_${FETCH_JOB_ID}.out"
tail -n 50 .nextflow.log
```

- **PENDING:** read the reason in `squeue`; the job has not started. Do not submit duplicate copies.
- **Unsupported Nextflow version:** load a module meeting the pinned pipeline's minimum version. Do not silently change the pipeline release just to bypass the error.
- **OUT_OF_MEMORY or a resource request exceeding the allocation:** compare `lab.config`, the SLURM request, and `seff`. Increasing task limits also requires enough memory in the outer allocation, plus room for Nextflow.
- **TIMEOUT:** the allocation ended before completion. Review progress and resubmit with `-resume` after correcting the wall-time request within partition limits.
- **Out of disk space:** SRA conversion can use much more space than the final compressed FASTQ. Check `df -h .`, `du -sh work results`, and the allocation quota; do not delete an active run's work directory.
- **Container-cache permission error:** check `echo "$SCRATCH"` and `ls -ld "$SCRATCH/.singularity_cache"`. The script creates and uses this exact directory. Fix its location or permissions before retrying.
- **Pipeline finished but FASTQ checks failed:** inspect the samplesheet, mappings, and parser error. Preserve the cache and fix the specific problem; do not treat existing output files as proof of success.

If you need help, provide the job ID, the relevant error, and the failed task's work directory from `.nextflow.log`. For a running job you intentionally need to stop, use `scancel JOBID` with that specific job ID.

## What to Keep

- `count_job.sh`, its job ID, and its output/error logs.
- The loaded module names and software versions in your notes.
- `lab.config`, `fetch/ids.csv`, `runs.tsv`, and `selected_runs.tsv`.
- `fetch_job.sh`, `fetch_job_id.txt`, and job logs.
- `results/`, `fastq_stats.tsv`, and the Nextflow cache until verification is complete.
- A short note explaining why these four runs are not a controlled comparison of platform alone.

Keep small scripts, configurations, and analysis notes in your coursework repository. Keep large reads, environments, caches, and generated results on course storage.

## References

- [Tufts HPC: running software with SLURM](https://rtguides.it.tufts.edu/hpc/access/60-brief-slurm.html)
- [SLURM sbatch reference, including environment export and script directives](https://slurm.schedmd.com/sbatch.html)
- [nf-core/fetchngs 1.13.0 source and runtime requirement](https://github.com/nf-core/fetchngs/blob/1.13.0/nextflow.config)
- [nf-core/fetchngs usage and identifier handling](https://nf-co.re/fetchngs/usage)
- [Tufts nf-core profile and its execution models](https://nf-co.re/configs/tufts/)
- [Nextflow caching and resuming](https://docs.seqera.io/nextflow/cache-and-resume/)
- [NCBI fasterq-dump: temporary storage and output formats](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
- [SeqKit usage and FASTQ statistics](https://bioinf.shenwei.me/seqkit/usage/)
- [ENA Illumina experiment SRX23612877](https://www.ebi.ac.uk/ena/browser/view/SRX23612877) and [PacBio experiment SRX14356846](https://www.ebi.ac.uk/ena/browser/view/SRX14356846)
