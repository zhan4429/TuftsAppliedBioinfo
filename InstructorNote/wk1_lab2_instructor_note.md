# Instructor note — Week 1 Session 2 lab

Planning guidance, checked against ENA and fetchngs 1.13.0 on 2026-11-03.

**30–40 minute class.** Students request one interactive allocation (two CPUs, 8 GB),
load modules, prepare the four-run input and batch script, and submit
one fetchngs batch job from that compute node. The counting-job exercise has been
removed. There is no second interactive allocation or student smoke-test download
during class. Budget 5–10 minutes for connection/allocation, 5 for modules, and 20–25
for job preparation and submission; queue waits may extend this. The batch job runs in
its own allocation. Students record its job ID manually with `squeue --me`, then exit
the interactive shell. Completion checks and FASTQ verification are after-class work.

**Library types matter.** SRR27956204 is RNA-seq from S288C under no stress;
SRR27956203 is RNA-seq under ethanol stress. Neither is a genomic Illumina control.
The selected BGISEQ, PacBio, and Nanopore libraries are WGS. Keep all four for a QC
teaching exercise only if students explicitly discuss this confounding. For a
controlled platform comparison, identify and verify a suitable Illumina WGS run
before changing the accession list. The Week 2 introduction also needs this caveat.

**Scale.** ENA currently advertises 31.40 Gbp and **21.67 GB** of compressed FASTQ
for the four selected runs, or about **650 GB for thirty individual downloads**.
Compute bytes from `fastq_bytes`, summing semicolon-separated entries; do not infer
compressed size from a fixed bases-to-bytes ratio. These are ENA file sizes, not
measured sizes after the lab's SRA conversion. The BGISEQ paired library has one
ENA FASTQ URL, so the lab deliberately uses `sratools` to obtain split reads.
SRA archives, uncompressed FASTQ, conversion scratch, copied results, and containers
increase peak disk use substantially. Check allocation quotas as well as filesystem
free space; published FASTQ size alone is not enough to budget the run.

If it is too much, two ways to cut it without losing the lesson:

*Drop to three platforms.* BGISEQ is the one to lose — it is a second short-read
technology. Keep the RNA-seq caveat for the retained Illumina library. Replace the
list in Section 3.1 with:

```bash
printf 'SRR27956204\nSRR18210286\nSRR17374240\n' > ids.csv
```

ENA currently advertises about 27.75 Gbp and 18.56 GB for these three runs. If you
change the list, also update the selection check, expected file count, and Week 2
inputs; the handout supplies four selected runs and asks students to confirm six
FASTQ files after downloading.

*Split the class.* Assign each student one platform, have them write the result into
an instructor-managed shared directory, and everyone uses the shared copies next
week. One complete shared set instead of thirty personal sets avoids most repeated
downloads. Assign permissions and ownership deliberately; do not have all students
write into the same Nextflow work directory. Stage a fallback either way.

**Week 2 will subsample regardless.** FastQC on a 44-million-read file takes far longer
than a lab slot. Next week opens by cutting each platform down to a workable size, which
is itself worth teaching. RNA-seq yield divided by genome size is not uniform
genomic coverage; interpret the short-read libraries according to their source.

**Stage fallbacks** at `$COURSE/shared/wk2/fallback_fastq/`, subsampled, one per platform.
Some downloads will fail and Week 2 cannot depend on Week 1 having worked.

**Runtime and resources.** The student handout pins fetchngs 1.13.0, whose manifest
requires Nextflow >=25.10.4. Confirm a compatible Nextflow/Java module and Singularity
on Pax before class. The instructor should run the fetchngs test before teaching,
using a separate test directory, `-profile test,singularity`, and
`--download_method sratools`. The test profile has its own resource settings; request
six CPUs and 24 GB to test within the same allocation size as the real job. The students'
two-CPU, 8 GB interactive allocation is only for checking software and submitting the
batch job. It is too small for the pipeline's default conversion step. Students should
not run the pipeline directly in that interactive shell.
The real job requests one task, six CPUs, 24 GB, and 24 hours, as specified for this
course. It uses the Singularity profile. Upstream fetchngs 1.13.0 advertises six CPUs
and 36 GB for the `process_medium` label used by `SRATOOLS_FASTERQDUMP`; the course's
24 GB allocation should be checked with the actual Pax setup and course inputs.
Validate the complete four-run workload with Pax's Nextflow settings before class. This is
the local-executor-inside-one-allocation model, not `-profile tufts` with additional
SLURM submissions. The 24-hour limit is an initial request, not a runtime guarantee.

Resource definitions: [fetchngs 1.13.0 defaults](https://github.com/nf-core/fetchngs/blob/1.13.0/conf/base.config)
and [FASTQ conversion process](https://github.com/nf-core/fetchngs/blob/1.13.0/modules/nf-core/sratools/fasterqdump/main.nf).

**Paths.** Confirm the work allocation
`/cluster/tufts/c27sp1bio0196` before teaching. The counting exercise and its shared
annotation path are no longer needed in Lab 2. Check the paths used by later handouts
as well. Give students the actual verified
`$MYWORK/week-01/session-02/fetch/results/fastq/` location rather than assuming those
examples point to the same files.

**Beginner workflow.** Keep `fetch_job.sh` focused on module loads, plain path
assignments, cache settings, and the Nextflow command. Avoid positional-argument
validation, grouped diagnostic redirection, shell arrays, and `shopt` in this lab.
SeqKit is available as a module on Pax; no Conda installation or activation is needed.

**Verification and retention.** Students run gzip testing and SeqKit statistics
as separate steps in an interactive compute session after the download finishes.
This lets them inspect the file count, errors, and statistics directly. Load SeqKit
with `module load seqkit`. Review mate counts and base totals;
archive `read_count` conventions and SRA-generated filenames need not match a
naive six-file lookup by run accession. Keep `work/` and `.nextflow/cache/` until
the real results are verified. For an instructor test that is no longer needed,
review `nextflow clean RUN_NAME -n -keep-logs` before using the corresponding `-f`
command. Cleanup is not part of the student class session. Deleting history is not
a routine requirement for a new run.
