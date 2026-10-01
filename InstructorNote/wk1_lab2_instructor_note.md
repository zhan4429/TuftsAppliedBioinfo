# Instructor note — Week 1 Session 2 lab

Removed from the student handout; kept here so the scale guidance is not lost.

**Scale.** Four platforms is about 31 Gbp, roughly 11 GB gzipped per student, so about
330 GB pulled from ENA for a class of thirty. That is the cost of making next week's
cross-platform comparison real. Check the allocation has room before Thursday.

If it is too much, two ways to cut it without losing the lesson:

*Drop to three platforms.* BGISEQ is the one to lose — it is a second short-read
technology and Illumina already covers that ground. Replace the selector in 4.2 with an
explicit list:

```bash
printf 'SRR27956204\nSRR18210286\nSRR17374240\n' > ids.csv
```

That is 27.8 Gbp, about 9.7 GB each.

*Split the class.* Assign each student one platform, have them write the result into
`$COURSE/shared/wk2/<platform>/`, and everyone uses the shared copies next week. One
quarter of the bandwidth, and it models how real collaborations actually divide work.
The cost is that one student's failure blocks others, so stage a fallback either way.

**Week 2 will subsample regardless.** FastQC on a 44-million-read file takes far longer
than a lab slot. Next week opens by cutting each platform down to a workable size, which
is itself worth teaching — nobody needs 1,100x coverage to find out that adapters are
present.

**Stage fallbacks** at `$COURSE/shared/wk2/fallback_fastq/`, subsampled, one per platform.
Some downloads will fail and Week 2 cannot depend on Week 1 having worked.

**Wall time is now 24 hours**, raised from 12. At 11 GB with thirty students competing for
ENA bandwidth, 12 was optimistic.
