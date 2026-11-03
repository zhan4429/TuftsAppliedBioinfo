# Instructor note — Week 1 Session 2 lab

Planning guidance, checked against ENA and fetchngs 1.13.0 on 2026-11-03.

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
list in 4.2 with:

```bash
printf 'SRR27956204\nSRR18210286\nSRR17374240\n' > ids.csv
```

ENA currently advertises about 27.75 Gbp and 18.56 GB for these three runs. If you
change the list, also update the selection check, expected file count, and Week 2
inputs; the current student script deliberately requires four platforms and six files.

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
on Pax before class. Test the actual course modules with the provided smoke test.
The real job requests one task, six CPUs, 18 GB, and 24 hours; `lab.config` limits
local workflow tasks to six CPUs and 12 GB in total, leaving JVM overhead. This is
the local-executor-inside-one-allocation model, not `-profile tufts` with additional
SLURM submissions. The 24-hour limit is an initial request, not a runtime guarantee.

**Paths.** Confirm the work allocation
`/cluster/tufts/yzhang85/bioinfo_course/appliedBio` and Session 1 annotation location
`/cluster/tufts/bio_appbio/shared/wk1` before teaching. They are separate settings in
the handout. Week 2 still uses the older course root; Week 3's example `HIFI` path
also differs from this lab's output. Give students the actual verified
`$MYWORK/week-01/session-02/fetch/results/fastq/` location rather than assuming those
examples point to the same files.

**Verification and retention.** Gzip testing and SeqKit parsing now run in the
download batch job, not on the login node. Review mate counts and base totals;
archive `read_count` conventions and SRA-generated filenames need not match a
naive six-file lookup by run accession. Keep `work/` and `.nextflow/cache/` until
the real results are verified. Cleanup examples now use a dry run and a named
completed test run; deleting history is not a routine requirement for a new run.
