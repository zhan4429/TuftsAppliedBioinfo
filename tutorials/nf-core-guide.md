# nf-core fundamentals on the Tufts HPC

Run a ready-made bioinformatics pipeline on the Tufts HPC cluster.

> In the [Nextflow guide](nextflow-guide.md) you connected processes and channels yourself. **nf-core provides ready-made pipelines built with Nextflow.** You choose a pipeline, provide your inputs, and let Nextflow run the analysis. This guide is for Applied Bioinformatics students and for anyone using the Tufts HPC.

---

## 1. Understand the main components

| Component                 | What it does                                                          |
| ------------------------- | --------------------------------------------------------------------- |
| **Nextflow**              | Runs the workflow and tracks completed tasks                          |
| **nf-core**               | Community-maintained, peer-reviewed bioinformatics pipelines          |
| **Slurm**                 | Schedules tasks on the cluster                                        |
| **Singularity**           | Runs each pipeline step inside a container, so you install nothing    |
| **Institutional profile** | The Tufts settings: `-profile tufts`                                  |

Your own pipelines in the Nextflow guide used environment modules. nf-core pipelines bring their own software as containers instead: each step names the exact image it was tested with, and Singularity runs it. You only need the `nextflow` and `singularity` modules.

Choose a pipeline for your question at [nf-co.re/pipelines](https://nf-co.re/pipelines), and read the documentation **for the release you will run**. Pipelines this course uses:

| Pipeline                                     | Release | Used for                                    |
| -------------------------------------------- | ------- | ------------------------------------------- |
| [nf-core/demo](https://nf-co.re/demo)         | 1.2.0   | This guide: a small, fast first run         |
| [nf-core/fetchngs](https://nf-co.re/fetchngs) | 1.13.0  | Downloading public data (Week 1)            |
| [nf-core/rnaseq](https://nf-co.re/rnaseq)     | 3.27.0  | RNA-seq: QC, alignment, quantification      |
| [nf-core/sarek](https://nf-co.re/sarek)       | 3.12.0  | Variant calling (Week 4 homework)           |

All four need **Nextflow 25.10.4 or later**.

### Find your way around a pipeline repository

Open the [nf-core/rnaseq repository](https://github.com/nf-core/rnaseq). These are the places to look when you want to understand or adjust a pipeline:

| File or folder                   | What it is for                                                                     |
| -------------------------------- | ---------------------------------------------------------------------------------- |
| `README.md` and `docs/`          | Overview, input requirements, usage, and a description of every output             |
| `main.nf`                        | The entry point that starts the workflow                                           |
| `nextflow.config`                | Default parameters, profiles, and links to the other configuration files           |
| `nextflow_schema.json`           | Every parameter, with its description and the values it accepts                    |
| `conf/base.config`               | Default CPU, memory and time requests, grouped by process label                    |
| `conf/modules/`                  | Settings for particular steps, such as STAR's arguments                            |
| `workflows/`, `subworkflows/`    | How the steps connect; subworkflows group related steps                            |
| `modules/nf-core/`               | Shared tool modules, such as STAR and FastQC, used across many pipelines           |
| `modules/local/`                 | Steps written for this pipeline only                                               |

On GitHub, use the branch and tag menu to select the release you run with `-r`. Files and settings change between releases.

---

## 2. Prepare your workspace

Choose a working directory:

```bash
# Applied Bioinformatics students:
export COURSE=/cluster/tufts/c27sp1bio0196
export NFCORE="$COURSE/$USER/nf-core-guide"

# Other Tufts HPC users: a directory in your lab's storage instead, for example
# export NFCORE=/cluster/tufts/<your_lab>/$USER/nf-core-guide

mkdir -p "$NFCORE/logs" && cd "$NFCORE"
```

Load the software and check the version:

```bash
module purge
module load nextflow singularity
nextflow -version          # 25.10.4 or later
singularity --version
```

### Where container images are kept

Pipelines download two kinds of files: the pipeline code, which is small, and the container images its steps run in. Images are large, and the same image is reused by many runs and many pipelines, so keep one cache for all of them. Nextflow uses the folder named by `NXF_SINGULARITY_CACHEDIR`. Every run script in this guide sets it inside the script, as the labs do, so a job never depends on the shell it was submitted from:

```bash
# Applied Bioinformatics students: the cache the labs use
export NXF_SINGULARITY_CACHEDIR="/cluster/tufts/c27sp1bio0196/$USER/.singularity_cache"

# Other Tufts HPC users: a folder in your lab's storage instead, for example
# export NXF_SINGULARITY_CACHEDIR=/cluster/tufts/<your_lab>/$USER/.singularity_cache
```

Course students: this is the cache Week 1's Lab 2 and the later labs use, so images downloaded there are reused. Other users: change that line in each script.

| Contents                        | Location                                     |
| ------------------------------- | -------------------------------------------- |
| Pipeline code                   | `~/.nextflow/assets/nf-core/<pipeline>/`     |
| Container images, all pipelines | `$NXF_SINGULARITY_CACHEDIR`                  |
| Task files and final outputs    | `work/` and `results/` in your run directory |

The first run of a pipeline is slow while it downloads; later runs reuse the downloads. If the cache is deleted, images are simply downloaded again.

### Read a launch command

```bash
nextflow run nf-core/demo -r 1.2.0 \
  -profile test,tufts \
  --outdir results
```

| Option              | Meaning                                                              |
| ------------------- | -------------------------------------------------------------------- |
| `nf-core/demo`      | The pipeline; Nextflow fetches it from GitHub                        |
| `-r 1.2.0`          | Run exactly this release. Always pin one.                            |
| `-profile test,tufts` | Use the pipeline's test data, and the Tufts cluster settings       |
| `--outdir results`  | Write final results here                                             |

**One dash** for Nextflow's options (`-r`, `-profile`, `-resume`, `-c`). **Two dashes** for the pipeline's parameters (`--input`, `--outdir`). Profiles combine left to right, separated by commas, with no spaces.

`test` makes the pipeline use its own small example inputs. For your own analysis, drop `test` and supply the required inputs (section 8).

---

## 3. Two ways to run on the Tufts HPC

| | **A. `-profile tufts`** | **B. `-profile singularity`, one allocation** |
| --- | --- | --- |
| Where tasks run | Each task is its own Slurm job | All tasks share one Slurm job |
| Your job script requests | A small job for the Nextflow manager | Enough for the largest task, plus overhead |
| Good for | Real analyses with many or large tasks | Small runs; when the queue is busy; teaching labs |
| `work/` after success | **Deleted automatically** | Kept |

Both use Singularity containers. Week 1's Lab 2 used route B for nf-core/fetchngs.

### A. The Tufts profile

`-profile tufts` loads the Tufts settings from [nf-core/configs](https://github.com/nf-core/configs/blob/master/conf/tufts.config). It sets:

| Setting                         | Value                                                              |
| ------------------------------- | ------------------------------------------------------------------ |
| Executor                        | `slurm`: every task becomes its own Slurm job                      |
| Partition                       | `batch`; change it with `--partition <name>`                       |
| GPU steps                       | Steps labelled `process_gpu` go to the `gpu` partition (section 7) |
| Largest request allowed per task | 64 CPUs, 500 GB, 168 hours; larger requests are reduced to these  |
| Concurrent task jobs            | At most 16, submitted no faster than one per ten seconds           |
| Containers                      | Singularity, enabled automatically                                 |
| Reports                         | Execution report, timeline, trace and pipeline diagram             |
| **Cleanup**                     | **`work/` is deleted when a run succeeds**; keep it with `--disable_cleanup true` |

Save this as `run_demo_tufts.sh`:

```bash
cat > run_demo_tufts.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=nfcore-demo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out

set -euo pipefail
module purge
module load nextflow singularity
export NXF_OPTS='-Xms512m -Xmx2g'

# Container cache: the one the labs use. Other Tufts users: a folder in your lab's storage.
export NXF_SINGULARITY_CACHEDIR="/cluster/tufts/c27sp1bio0196/$USER/.singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

nextflow run nf-core/demo -r 1.2.0 \
  -profile test,tufts \
  --outdir "$PWD/results"
EOF

sbatch run_demo_tufts.sh
squeue -u "$USER"
```

The job runs only the **Nextflow manager**, so it needs little CPU or memory. Its tasks appear in `squeue` as separate jobs, each with its own request. The manager must keep running until the pipeline finishes, so give it a time limit that covers the whole run, including queueing.

The script sets the container cache itself, as the labs do, so the job does not depend on the shell you submit from.

### B. One allocation, with Singularity

Here everything runs inside the job itself, using Nextflow's `local` executor. Write a small configuration that tells the pipeline how much of the job it may use:

```bash
cat > local.config <<'EOF'
process {
    executor       = 'local'
    resourceLimits = [cpus: 4, memory: 12.GB, time: 4.h]
}
executor {
    cpus   = 4
    memory = 12.GB
}
EOF
```

`resourceLimits` caps any single task; `executor.cpus` and `memory` cap all running tasks together. Both are needed: without them, the pipeline's default requests can exceed your job, and a per-task limit alone does not stop several tasks from overfilling it. Keep them below the job's request, which leaves room for Nextflow itself.

```bash
cat > run_demo_local.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=nfcore-demo-local
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out

set -euo pipefail
module purge
module load nextflow singularity
export NXF_OPTS='-Xms512m -Xmx2g'

# Container cache: the one the labs use. Other Tufts users: a folder in your lab's storage.
export NXF_SINGULARITY_CACHEDIR="/cluster/tufts/c27sp1bio0196/$USER/.singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

nextflow run nf-core/demo -r 1.2.0 \
  -profile test,singularity \
  -c local.config \
  --outdir "$PWD/results_local"
EOF

sbatch run_demo_local.sh
```

**For real analyses, size the job to the work.** With route B, the job must hold the most demanding task and anything running beside it. With route A, the manager's size does not matter: requesting more CPUs for the manager does not give STAR more CPUs.

---

## 4. Check results and resume

In your run directory:

| File or folder               | What to check                                                   |
| ---------------------------- | --------------------------------------------------------------- |
| `logs/<job>_<id>.out`        | Progress, and the first error if something failed               |
| `results/multiqc/`           | The MultiQC report: open this first                             |
| `results/pipeline_info/`     | Execution report, timeline and trace: what ran, how long, how much memory |
| `work/`                      | Task folders (deleted after success with `-profile tufts`)      |
| `.nextflow.log`              | Full detail when a run fails                                    |

`results/pipeline_info/` is worth a look on every run. Its execution report shows each task's peak memory and CPU use, which is exactly the evidence you need for section 5.

### Resuming, and the cleanup rule

Add `-resume` to the `nextflow run` line and resubmit. Nextflow reuses every unchanged task that finished, and reruns the rest. Keep the same run directory: resume needs both `.nextflow/` and `work/`.

**With `-profile tufts`, a run that succeeds deletes its `work/` directory.** That saves a lot of space, and it is safe because your results are published. But it means:

- `-resume` rescues a run that **failed**, because a failed run's `work/` is kept.
- After a run **succeeds**, `-resume` has nothing to reuse: change one parameter and every step runs again.
- While you are still adjusting parameters, or need to inspect intermediate files, add `--disable_cleanup true`.

In the Nextflow guide you deleted `work/` by hand with `nextflow clean` and watched `-resume` rerun everything. This is the same thing, done for you.

### When a task fails

The log names the failed process and its work directory. In that directory, `.command.err` and `.command.log` hold the tool's own error, and `.command.sh` the exact command. Two common causes on a cluster:

- **Exit code 137, or the job ran out of memory:** the step needs more memory than its default (section 5).
- **Exit code 140, or the task hit its time limit:** the step needs more time.

nf-core pipelines retry many such failures once with a larger request, so a single failed attempt is not always a problem.

---

## 5. Adjust a step's resources

The Tufts profile handles the cluster; the pipeline supplies each step's CPU, memory and time. If a step runs out of memory or time, find its defaults and override them in your own `custom.config`. Never edit the pipeline's files.

### Example: STAR_ALIGN in nf-core/rnaseq 3.27.0

Open [`modules/nf-core/star/align/main.nf`](https://github.com/nf-core/rnaseq/blob/3.27.0/modules/nf-core/star/align/main.nf). Near the top is the process name and its label:

```groovy
process STAR_ALIGN {
    label 'process_high'
    // ...
}
```

Then open [`conf/base.config`](https://github.com/nf-core/rnaseq/blob/3.27.0/conf/base.config) and find that label:

```groovy
withLabel:process_high {
    cpus   = { 12    * task.attempt }
    memory = { 72.GB * task.attempt }
    time   = { 16.h  * task.attempt }
}
```

The path is always the same: **process definition → its label → that label's settings in `conf/base.config`.** STAR's own arguments live separately, in [`conf/modules/align_star.config`](https://github.com/nf-core/rnaseq/blob/3.27.0/conf/modules/align_star.config).

### How `* task.attempt` changes requests on a retry

`task.attempt` is 1 the first time a task runs and goes up each time Nextflow retries it, so each retry asks for more:

| Attempt              | CPUs | Memory | Time limit |
| -------------------- | ---- | ------ | ---------- |
| 1: first try         | 12   | 72 GB  | 16 hours   |
| 2: first retry       | 24   | 144 GB | 32 hours   |
| 3: second retry      | 36   | 216 GB | 48 hours   |

A retry happens only when the error policy allows it. In rnaseq 3.27.0:

```groovy
errorStrategy = { task.exitStatus in ((130..145) + 104 + (175..177)) ? 'retry' : 'finish' }
maxRetries    = 1
```

Only those exit codes, which include running out of memory or time, trigger a retry, and `maxRetries = 1` allows one, so two attempts in total. On the Tufts cluster every request is also capped at 64 CPUs, 500 GB and 168 hours, whatever the attempt.

`task.attempt` counts retries of one task inside one run. Launching the whole pipeline again does not make it 2.

### Override by label or by process name

Save this as `custom.config` in your run directory, and keep only the block you need:

```bash
cat > custom.config <<'EOF'
process {
    // Every step labelled process_high: bigger requests, still growing on retry.
    withLabel: 'process_high' {
        cpus   = { 16    * task.attempt }
        memory = { 96.GB * task.attempt }
        time   = { 24.h  * task.attempt }
    }

    // One particular step, by its process name: fixed requests.
    withName: 'STAR_ALIGN' {
        cpus   = 24
        memory = 128.GB
        time   = 24.h
    }
}
EOF
```

| Selector    | Effect of this example                                                                     |
| ----------- | ------------------------------------------------------------------------------------------ |
| `withLabel` | Every `process_high` step starts at 16 CPUs, 96 GB and 24 hours, and scales on retry        |
| `withName`  | STAR alignment asks for 24 CPUs, 128 GB and 24 hours, the same on every attempt             |

If you keep both, `withName` wins for STAR. Add `-c custom.config` to the `nextflow run` line to apply it.

**These numbers are examples for rnaseq on a large genome, not recommendations.** Base your requests on the execution report in `results/pipeline_info/`: it shows what each step actually used.

---

## 6. Reference genomes: the Tufts iGenomes copy

Many nf-core pipelines accept `--genome <name>` for a standard reference instead of your own FASTA and annotation files. The Tufts profile points those names at a copy on the cluster, `/cluster/tufts/biocontainers/datasets/igenomes/`, so nothing is downloaded. See what is there:

```bash
ls /cluster/tufts/biocontainers/datasets/igenomes/
ls /cluster/tufts/biocontainers/datasets/igenomes/Saccharomyces_cerevisiae/Ensembl/
```

For example, the yeast genome this course uses is `--genome R64-1-1` in nf-core/rnaseq.

**iGenomes annotations are old.** They are convenient for practice and for checking that a pipeline runs, but for an analysis you will report, download a current FASTA and GTF from the same source and release, such as Ensembl, and pass them with `--fasta` and `--gtf` (section 8).

---

## 7. Steps that need a GPU

A few pipelines have GPU-accelerated steps, such as Parabricks-based alignment. The Tufts profile sends any step labelled `process_gpu` to the `gpu` partition with one GPU, so you run the pipeline as usual:

```bash
nextflow run <pipeline> -profile tufts ...
nextflow run <pipeline> -profile tufts --gpu_partition <name> ...    # a different GPU partition
```

Whether any step uses a GPU depends on the pipeline and its options. Most nf-core pipelines run on CPUs only.

---

## 8. Prepare inputs for a real analysis

The `test` profile supplies its own inputs; for your data you supply them. Here is a paired-end **nf-core/rnaseq 3.27.0** example; the [rnaseq usage page](https://nf-co.re/rnaseq/3.27.0/docs/usage/) has the details.

### The samplesheet

Each sample's R1 and R2 files go on one row, and each biological replicate is its own sample. Save a plain CSV file, `samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,strandedness
control_rep1,/full/path/to/control_rep1_R1.fastq.gz,/full/path/to/control_rep1_R2.fastq.gz,auto
treated_rep1,/full/path/to/treated_rep1_R1.fastq.gz,/full/path/to/treated_rep1_R2.fastq.gz,auto
```

| Column         | What to provide                                                                                    |
| -------------- | -------------------------------------------------------------------------------------------------- |
| `sample`       | A name without spaces; each biological replicate separately                                        |
| `fastq_1`      | The full path to the sample's R1 file                                                              |
| `fastq_2`      | The full path to its R2 file; empty for single-end reads                                           |
| `strandedness` | `forward`, `reverse` or `unstranded` if you know it from the library kit; `auto` lets the pipeline infer it |

Keep the header exactly as shown. The [samplesheet schema](https://github.com/nf-core/rnaseq/blob/3.27.0/assets/schema_input.json) lists every accepted value.

**From public data, let nf-core/fetchngs write it.** Add `--nf_core_pipeline rnaseq` when you download with fetchngs, and it writes `results/samplesheet/samplesheet.csv` in exactly this format, with full paths and `strandedness` set to `auto`.

### The reference

For alignment with STAR and quantification with Salmon, provide a genome FASTA and a matching GTF annotation: the same assembly, ideally the same source and release. The pipeline builds the indexes it needs from them.

### The run script

```bash
cat > run_rnaseq.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=rnaseq
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out

set -euo pipefail
module purge
module load nextflow singularity
export NXF_OPTS='-Xms512m -Xmx2g'

# Container cache: the one the labs use. Other Tufts users: a folder in your lab's storage.
export NXF_SINGULARITY_CACHEDIR="/cluster/tufts/c27sp1bio0196/$USER/.singularity_cache"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

nextflow run nf-core/rnaseq -r 3.27.0 \
  -profile tufts \
  --input samplesheet.csv \
  --fasta /full/path/to/genome.fa.gz \
  --gtf /full/path/to/genes.gtf.gz \
  --aligner star_salmon \
  --outdir "$PWD/results" \
  -resume
EOF

sbatch run_rnaseq.sh
```

`--aligner star_salmon` selects STAR alignment with Salmon quantification. Add `-c custom.config` to apply resource overrides from section 5, and `--disable_cleanup true` while you are still adjusting parameters.

The manager job has two days so that it outlives every task's time in the queue; adjust that to your analysis. Copy final results to backed-up storage when the run finishes.

---

## After the guide

### Ask for help

Join the [nf-core Slack](https://nf-co.re/join) and ask in your pipeline's channel, such as `#rnaseq`, or in `#help`. Include the pipeline and its release, your Nextflow version, the profile, the full launch command, and the error message. Remove sensitive sample information before sharing logs.

For accounts, storage, partitions or the scheduler, contact [Tufts Research Technology](https://it.tufts.edu/high-performance-computing).

### Useful links

- [nf-core pipelines](https://nf-co.re/pipelines)
- [Tufts institutional profile](https://github.com/nf-core/configs/blob/master/conf/tufts.config) and [its documentation](https://github.com/nf-core/configs/blob/master/docs/tufts.md)
- [nf-core usage documentation](https://nf-co.re/docs/usage/getting_started/introduction)
- [Nextflow fundamentals on the Tufts HPC](nextflow-guide.md)
