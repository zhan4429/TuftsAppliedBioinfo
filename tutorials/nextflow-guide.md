# Nextflow fundamentals on the Tufts HPC

From a Bash loop to a reproducible workflow on the Tufts HPC cluster.

> This guide is for Applied Bioinformatics students and for anyone using the Tufts HPC. Course students will meet Nextflow in several labs; work through this guide once and those labs will make sense. When you are done, continue with the [nf-core guide](nf-core-guide.md).

---

## Start with the big picture

You already know how to run a tool. A workflow manager helps you connect tools, track their inputs and outputs, and pick up where you left off.

In this guide you run the same analysis three ways: a **Bash loop**, **Nextflow inside one allocation**, and **Nextflow submitting Slurm jobs**. The tools and the input data stay the same.

| Step        | What it does            | What comes out                             |
| ----------- | ----------------------- | ------------------------------------------ |
| **fastp**   | Trim and filter reads   | Trimmed FASTQ files, JSON and HTML reports |
| **FastQC**  | Check the trimmed reads | ZIP and HTML quality reports               |
| **MultiQC** | Combine the QC reports  | One HTML summary                           |

Nextflow describes this as a **dataflow**: FastQC on a file starts as soon as that file's trimmed reads are ready, and MultiQC waits for every report, including the JSON files from fastp.

### Five words you will meet

| Word         | Meaning in this example                                         |
| ------------ | --------------------------------------------------------------- |
| **Process**  | A tool step, such as `FASTP` or `FASTQC`                        |
| **Task**     | One execution of a process, for one input                       |
| **Channel**  | A stream of files or values connecting steps                    |
| **Workflow** | The block that connects processes and channels                  |
| **Executor** | Where tasks run: inside one allocation, or as separate Slurm jobs |

Nextflow and Snakemake both run HPC workflows with containers or modules. Nextflow connects processes through channels in a Groovy-based language; Snakemake connects rules through input and output files in Python. This course uses Nextflow, because the nf-core pipelines are built on it.

---

## Prepare your workspace

Everything runs on compute nodes through `sbatch`. The login node is only for editing files and submitting jobs.

Choose a working directory:

```bash
# Applied Bioinformatics students:
export COURSE=/cluster/tufts/c27sp1bio0196
export NFDIR="$COURSE/$USER/nextflow-guide"

# Other Tufts HPC users: a directory in your lab's storage instead, for example
# export NFDIR=/cluster/tufts/<your_lab>/$USER/nextflow-guide

mkdir -p "$NFDIR/logs" && cd "$NFDIR"
```

Check that Nextflow is available, and its version:

```bash
module load nextflow
nextflow -version
```

The nf-core pipelines used in this course need **Nextflow 25.10.4 or later**. If the default module is older, `module avail nextflow` lists the versions installed; load a newer one by name.

### Download the example data

`setup.sh` downloads eleven small FASTQ files from the nf-core test datasets, at a fixed revision so everyone gets identical files, and checks each one. A small download like this is fine on the login node.

```bash
cat > setup.sh <<'EOF'
#!/bin/bash
set -euo pipefail

mkdir -p data
cd data

# A fixed revision of the nf-core test dataset, so everyone gets identical files.
base_url=https://raw.githubusercontent.com/nf-core/test-datasets/cad884e7fbe5617dcd55bdd07d9ddf65b54febe4/testdata/GSE110004

for file in SRR635{7070..7076}_1.fastq.gz SRR635{7070,7071,7072,7076}_2.fastq.gz; do
    if [[ ! -s "$file" ]]; then
        curl -fL --retry 3 "$base_url/$file" -o "$file.part"
        gzip -t "$file.part"
        mv "$file.part" "$file"
    fi
    gzip -t "$file"
done
ls -lh
EOF

bash setup.sh
```

The script stops if a download or a compressed-file check fails, and skips files it already has, so running it again is safe.

> **A teaching example, one file at a time.** The download is seven R1 files and four R2 files: 11 files, not 11 biological samples. We process each file independently in single-end mode, so the outputs are **not synchronised paired-end reads**. For real paired-end data, run fastp on both mates together with `-i/-I/-o/-O`.

Your workspace will hold:

| Contents            | Folder     |
| ------------------- | ---------- |
| Downloaded FASTQ    | `data/`    |
| Published reports   | `results/` |
| Nextflow task files | `work/`    |
| Slurm logs          | `logs/`    |

---

## Begin with a Bash loop

The familiar version processes one FASTQ file at a time: trim, check quality, move on, then combine the reports.

```bash
cat > qc_forloop.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=qc_forloop
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out

# One FASTQ at a time, in single-end mode. Resubmitting starts over from the first file.
set -euo pipefail
module purge
module load fastp fastqc multiqc

OUTDIR=results/forloop
mkdir -p "$OUTDIR/fastp" "$OUTDIR/fastqc" "$OUTDIR/multiqc"

for fq in data/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    echo "Processing $sample"
    fastp -w 2 -i "$fq" \
        -o "$OUTDIR/fastp/${sample}.trimmed.fastq.gz" \
        -j "$OUTDIR/fastp/${sample}.fastp.json" \
        -h "$OUTDIR/fastp/${sample}.fastp.html"
    fastqc -t 2 -o "$OUTDIR/fastqc" "$OUTDIR/fastp/${sample}.trimmed.fastq.gz"
done

multiqc --force -o "$OUTDIR/multiqc" "$OUTDIR/fastp" "$OUTDIR/fastqc"
EOF

sbatch qc_forloop.sh
```

This works, and it is how most analyses start. Notice what it cannot do. Everything runs one file after another, even though the files are independent. If the job dies on file nine, resubmitting redoes files one to eight. And the script itself spells out every filename and the order of every step. Nextflow takes over all three.

---

## Connect the steps with Nextflow

A Nextflow pipeline is two files: `main.nf` says *what* to run, and `nextflow.config` says *how* and *where*.

### The pipeline: `main.nf`

```bash
cat > main.nf <<'EOF'
// Read QC with Nextflow: fastp -> FastQC -> MultiQC

params.fastq_dir = 'data'
params.outdir    = 'results'

// Each FASTQ file is processed on its own, in single-end mode.
// The outputs are NOT synchronised paired-end reads.

// Step 1: trim reads with fastp
process FASTP {
    tag "${reads.simpleName}"
    cpus 2
    memory '4 GB'
    time '30m'
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    path reads

    output:
    path '*.trimmed.fastq.gz', emit: trimmed
    path '*.fastp.json',       emit: json
    path '*.fastp.html'

    script:
    """
    fastp -w ${task.cpus} -i "${reads.name}" \\
        -o "${reads.simpleName}.trimmed.fastq.gz" \\
        -j "${reads.simpleName}.fastp.json" \\
        -h "${reads.simpleName}.fastp.html"
    """
}

// Step 2: FastQC on the trimmed reads
process FASTQC {
    tag "${reads.simpleName}"
    cpus 2
    memory '4 GB'
    time '30m'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path reads

    output:
    path '*_fastqc.zip', emit: zip
    path '*_fastqc.html'

    script:
    """
    fastqc -t ${task.cpus} "${reads.name}"
    """
}

// Step 3: combine every report with MultiQC
process MULTIQC {
    cpus 1
    memory '4 GB'
    time '30m'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path reports

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    reads_ch = channel.fromPath("${params.fastq_dir}/*.fastq.gz", checkIfExists: true)

    FASTP(reads_ch)
    FASTQC(FASTP.out.trimmed)
    MULTIQC(FASTP.out.json.mix(FASTQC.out.zip).collect(sort: true))
}
EOF
```

Each process declares its resources, inputs, outputs and command:

- `path` stages an input file into the task's own working directory.
- `emit` names an output channel, so the workflow can refer to it.
- `publishDir` copies the declared outputs into the results folder.
- `${task.cpus}` passes the process's CPU request to the tool, so the two always agree.
- `cpus`, `memory` and `time` are the task's request. On the cluster, each one is what Slurm reserves for that task.
- `tag` labels each task in the log with its file name, so you can tell the eleven FASTP tasks apart.

The `workflow` block is the connection map. `FASTP.out.trimmed` feeds FastQC. fastp's JSON reports and FastQC's ZIP reports feed MultiQC: `mix()` merges the two streams, and `collect(sort: true)` waits for every report and puts them in a stable order. Downstream steps take their inputs from channels, never from the results folder.

> **Try it:** which process runs only once? **MultiQC.** With 11 input files you get 11 FASTP tasks, 11 FASTQC tasks and one MultiQC task: **23 tasks** in total.

### The configuration: `nextflow.config`

```bash
cat > nextflow.config <<'EOF'
// Software: every tool is an environment module on the Tufts HPC.
// Nextflow runs `module load <name>` at the start of each task.
process {
    withName: 'FASTP'   { module = 'fastp' }
    withName: 'FASTQC'  { module = 'fastqc' }
    withName: 'MULTIQC' { module = 'multiqc' }
}

// Where tasks run. Choose one with -profile local or -profile slurm.
profiles {
    // Every task runs inside the single Slurm allocation that runs Nextflow.
    local {
        process.executor = 'local'
        executor {
            cpus   = 4
            memory = 12.GB
        }
    }
    // Every task becomes its own Slurm job in the batch partition.
    slurm {
        process.executor = 'slurm'
        process.queue    = 'batch'
        executor {
            queueSize       = 16
            submitRateLimit = '10 sec'
        }
    }
}
EOF
```

Nextflow reads `nextflow.config` automatically from the directory you launch in.

**Software.** The `module` setting tells Nextflow to run `module load fastp` (and so on) at the start of each task, whether the task runs locally or as its own Slurm job. Each task loads only what it needs. To pin a version, name it: `module = 'fastp/0.23.4'`, using a version listed by `module avail fastp`. Keeping software in the configuration, not in `main.nf`, means the same pipeline can run elsewhere with containers just by changing the configuration.

**Profiles.** A profile is a named block of settings you switch on with `-profile`. This file has two, one per way of running. You will meet the same idea in nf-core, where `-profile test,tufts` switches on a test dataset and the Tufts cluster settings.

| Setting                          | What it does                                                               |
| -------------------------------- | -------------------------------------------------------------------------- |
| `executor = 'local'`             | Tasks run inside the allocation that runs Nextflow                         |
| `executor.cpus`, `memory`        | The budget the local tasks share; keep it inside the allocation            |
| `executor = 'slurm'`             | Each task is submitted as its own Slurm job                                |
| `queue = 'batch'`                | The partition for those jobs                                               |
| `queueSize = 16`                 | At most 16 task jobs queued or running at once                             |
| `submitRateLimit = '10 sec'`     | One submission per ten seconds, so the scheduler is not flooded            |

The `slurm` profile mirrors what the Tufts nf-core profile sets for nf-core pipelines.

---

## Choose where tasks run

The pipeline stays the same. The profile changes where its tasks run.

### Local: everything inside one allocation

The launcher requests four CPUs and 16 GB, and the `local` profile lets tasks share four CPUs and 12 GB of it. The difference leaves room for Nextflow itself, which runs on Java. **If you change `--cpus-per-task`, change `executor.cpus` to match.**

```bash
cat > qc_nextflow_local.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=qc_nextflow_local
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out

set -euo pipefail
module purge
module load nextflow

export NXF_OPTS='-Xms512m -Xmx2g'

nextflow run main.nf -profile local \
    --outdir results/local -w work/local -resume -ansi-log false
EOF
```

With four CPUs and two per task, two tasks run at a time.

### Slurm: a separate job for each task

Here the launcher only runs the Nextflow **manager**, which submits each fastp, FastQC and MultiQC task as its own Slurm job. The manager needs little CPU or memory, but it must outlive the whole pipeline, including time its tasks spend waiting in the queue.

```bash
cat > qc_nextflow_slurm.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=qc_nextflow_slurm
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out

set -euo pipefail
module purge
module load nextflow

export NXF_OPTS='-Xms512m -Xmx2g'

nextflow run main.nf -profile slurm \
    --outdir results/slurm -w work/slurm -ansi-log false
EOF
```

### Requested and allocated resources

| Where it is set            | What it describes                        | Local profile                         | Slurm profile                      |
| -------------------------- | ---------------------------------------- | ------------------------------------- | ---------------------------------- |
| `#SBATCH` in the launcher  | The job that runs Nextflow               | Shared by all tasks                   | Only the manager                   |
| `cpus`, `memory`, `time` in `main.nf` | One task                      | A share of the launcher's allocation  | That task's own Slurm job          |

`NXF_OPTS` caps the memory Nextflow's Java process may use.

### Submit the examples

From your workspace, submit them one at a time:

```bash
cd "$NFDIR"
sbatch qc_nextflow_local.sh
squeue -u "$USER"
```

When it finishes, submit the Slurm version and watch the task jobs appear and disappear:

```bash
sbatch qc_nextflow_slurm.sh
squeue -u "$USER"
```

> **Why the Slurm launcher has no `-resume`.** Nextflow keeps one run history per directory, shared by both profiles. With `-resume`, the Slurm run would find the local run's finished tasks, reuse all 23, and submit nothing, which defeats the demonstration. Separate `-w` work folders do not separate the history.

---

## Publish outputs with publishDir

Tasks write their files under `work/`. `publishDir` exposes the declared outputs in a results folder. FASTP uses:

```groovy
publishDir "${params.outdir}/fastp", mode: 'copy'
```

With `--outdir results/local`, Nextflow creates `results/local/fastp/` when it publishes; no `mkdir` is needed.

| Mode                | What appears in results             | Storage and cleanup                                                         |
| ------------------- | ----------------------------------- | --------------------------------------------------------------------------- |
| `copy`              | Independent copies                  | Uses extra space; survives deleting `work/`                                 |
| `symlink` (default) | Absolute symbolic links             | Little extra space; links break when `work/` is cleaned                     |
| `rellink`           | Relative symbolic links             | Survives moving both folders together                                       |
| `link`              | Hard links                          | Same filesystem only; survives deleting `work/`                             |
| `move`              | Files moved out of `work/`          | Only for final outputs nothing downstream needs                             |
| `copyNoFollow`      | Copies that keep symbolic links     | The links' targets must stay available                                      |

**Why this guide uses `copy`:** the reports are small, and copies survive cleaning `work/`. With the default, `symlink`, cleaning `work/` would leave a results folder full of broken links.

Publishing happens in the background. Pass outputs to later steps through channels, never by reading the results folder.

---

## Find your results

| Route          | Combined report                               |
| -------------- | --------------------------------------------- |
| Bash loop      | `results/forloop/multiqc/multiqc_report.html` |
| Nextflow local | `results/local/multiqc/multiqc_report.html`   |
| Nextflow Slurm | `results/slurm/multiqc/multiqc_report.html`   |

The `fastp/` and `fastqc/` folders beside each `multiqc/` hold the individual reports. Open the HTML report through Open OnDemand (**Files**), or copy it to your computer.

```bash
ls results/*/multiqc/multiqc_report.html
du -sh work results
```

A run that finished without errors still needs its reports checked: make sure they describe the eleven files you expected.

---

## Resume and troubleshoot

### Resume

Resubmit the local launcher. It includes `-resume`, so all 23 tasks come back as `Cached` and nothing reruns. Change one input file, or the command in one process, and only the affected tasks run again.

Resume needs **both** `.nextflow/` (the record of past tasks) and `work/` (their outputs). The published results alone are not enough. Changed inputs, scripts or software versions make tasks run again.

To resume a run other than the most recent one, give its name: `-resume RUN_NAME`, using a name from `nextflow log`.

### Read the run history

From your workspace:

```bash
module load nextflow
nextflow log
nextflow log last -f name,status,exit,workdir
```

`last` is the most recent run launched in this directory, whether it succeeded or failed; after switching between the two launchers, check `nextflow log` to see which that is. `-f` chooses columns: `name` is the task, `status` its state, `exit` its exit code (`0` is success), and `workdir` its task directory.

| Command                     | Shows                                            |
| --------------------------- | ------------------------------------------------ |
| `nextflow log`              | Every run: names, session IDs, status            |
| `nextflow log last`         | Task directories of the latest run               |
| `nextflow log RUN_NAME`     | Task directories of one particular run           |
| `nextflow log -quiet`       | Run names only                                   |
| `nextflow log -l`           | The fields you can pass to `-f`                  |

Only the failed tasks of the latest run:

```bash
nextflow log last -F "status == 'FAILED'" -f name,status,exit,workdir
```

**Lowercase `-f` chooses fields; uppercase `-F` filters tasks.** A run that failed before any task started has no tasks to show; read `.nextflow.log` instead.

### Follow a failure to its task folder

Start with the launcher's log in `logs/` and with `.nextflow.log`. The error names the failed task and its work directory. In that directory:

| File           | What to look for                                  |
| -------------- | ------------------------------------------------- |
| `.command.sh`  | The exact command the task ran                    |
| `.command.err` | The tool's error messages                         |
| `.command.log` | Everything the task printed                       |
| `.exitcode`    | The task's exit status                            |
| `.command.run` | The wrapper: module loading, staging, resources   |

| Symptom                                  | Usual cause and fix                                                                                          |
| ---------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| Exit code `127`, "command not found"     | The module did not load. Check its name with `module avail <tool>` and the `module` line in `nextflow.config` |
| Exit code `137` or `OUT_OF_MEMORY`       | The task exceeded its `memory`. Raise it in `main.nf` and resume                                             |
| Exit code `140`, or the task timed out   | The task exceeded its `time`. Raise it and resume                                                            |
| Slurm refuses the task jobs              | A request exceeds the partition's limits; check the `cpus`, `memory` and `time` directives                  |
| The manager job hits its time limit      | Raise `--time` in the launcher and resubmit with `-resume`                                                   |
| A tool cannot open an input that exists  | Rare with modules; add `process.stageInMode = 'copy'` to `nextflow.config` and resume                        |

To debug a tool by hand, start an interactive job and run the task's own script from its directory:

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=01:00:00 --pty bash
cd work/local/ab/cdef12...        # the workdir from the error or from nextflow log
bash .command.run
```

### One dash or two?

**One dash** for Nextflow's own options: `-resume`, `-profile`, `-c`, `-r`, `-w`. **Two dashes** for the pipeline's parameters: `--fastq_dir`, `--outdir`. `--resume` is not the resume option; Nextflow would treat it as a pipeline parameter and ignore it.

---

## Save your work, then clean up

`work/` grows quickly on real data, and keeping it is what makes `-resume` possible. Clean it only when a run is finished and its results are safe.

First save what matters: the published `results/`, the pipeline files, and the logs. Then preview what would be removed:

```bash
nextflow clean -n $(nextflow log -quiet)
```

`-n` only lists what would go. `$(nextflow log -quiet)` names **every** run in this directory, and that matters: once one run has resumed another, they share task directories, and Nextflow removes a shared directory only when every run that uses it is named. Naming one run alone can list nothing at all.

When the preview looks right:

```bash
nextflow clean -f $(nextflow log -quiet)
du -sh work results
```

Published copies survive, but every cleaned task is gone: a later `-resume` will rerun all of them. The Tufts nf-core profile does this automatically after each successful run, which the [nf-core guide](nf-core-guide.md) comes back to.

---

## Keep learning

- [nf-core fundamentals on the Tufts HPC](nf-core-guide.md): run community pipelines with the Tufts profile
- [Nextflow training](https://training.nextflow.io/): hands-on lessons, from the basics up
- [Caching and resuming](https://docs.seqera.io/nextflow/cache-and-resume): what makes a task reusable
- [Process directives](https://docs.seqera.io/nextflow/reference/process): every setting a process accepts
- Tufts HPC questions about accounts, storage or the scheduler: [Tufts Research Technology](https://it.tufts.edu/high-performance-computing)
