# Lab 10B — Shotgun Profiles, Community Statistics, and MAG Quality

**Applied Bioinformatics · Week 10, Session 2**
**Runs alongside the lecture**

> You will profile one shotgun sample yourself, twice, and see what a single parameter does to the species list. Then you will take the whole cohort through Bracken, diversity statistics and differential abundance, and finish by judging the genomes an assembly pipeline recovered from it.

---

## Before you start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED=$COURSE/shared/wk10
export D=$SHARED/shotgun
export DB=$SHARED/k2_standard_08gb

mkdir -p $MYWORK/week-10/session-02/bracken && cd $MYWORK/week-10/session-02

module load kraken2 bracken r
module list
# the R scripts below need four packages; each should print TRUE
Rscript -e 'for (p in c("vegan", "ALDEx2", "ANCOMBC", "TreeSummarizedExperiment")) cat(p, requireNamespace(p, quietly = TRUE), "\n")'
ls $D
```

You will find one sample's reads (`s01_R1`, `s01_R2`, already host-removed), Kraken2 reports for the whole cohort (made with `--confidence 0.1`), the cohort's `metadata.tsv`, and a CheckM2 report from an nf-core/mag run.

---

## Part 1 — Profile one sample, twice

Kraken2 with its default confidence, then with the confidence the lecture recommended:

```bash
for c in 0 0.1; do
  kraken2 --db $DB --paired --threads 8 --memory-mapping \
    --confidence $c --report s01_c$c.kreport --output /dev/null \
    $D/s01_R1.fastq.gz $D/s01_R2.fastq.gz 2> s01_c$c.log
done
```

`--output /dev/null` discards the per-read assignments, which are enormous and rarely needed. The report is what you keep.

```bash
for c in 0 0.1; do
  awk -F'\t' -v c=$c '$4=="U" {u=$2} $4=="R" {r=$2}
    $4=="S" {n++; if ($2>=10) m++}
    END {printf "confidence %-4s classified %5.1f%%   species %4d   with 10+ reads %4d\n", c, 100*r/(u+r), n, m}' s01_c$c.kreport
done
```

The report's columns are: percentage of reads, reads in the clade, reads assigned to this exact node, rank code (`U` unclassified, `R` root, `S` species), NCBI taxonomy ID, and the indented name.

Now look at the tail of the default run:

```bash
awk -F'\t' '$4=="S" && $2<10 {sub(/^ +/, "", $6); print $2"\t"$6}' s01_c0.kreport | sort -n | tail -15
```

**Question 1.** How many species does the default run report, and how many have ten or more reads? What happened to the rest at confidence 0.1?

::: Answer
Most of the default species list sits at single-digit read counts: species whose k-mers a few reads happened to match. Confidence 0.1 requires part of each read's k-mers to agree with the assignment, and the tail shrinks sharply.

Read the tail you printed. If it contains a famous pathogen with a handful of reads, you have not found it: that is k-mer sharing between relatives, and it is how clinical metagenomics has embarrassed itself before. Neither threshold is "correct"; state the one you used, and apply an abundance floor as well.
:::

**Question 2.** The classified percentage also fell. Is that a loss?

::: Answer
Only on paper. Raising the confidence moves reads with weak, ambiguous evidence to unclassified instead of forcing them onto a guess. The reads you lost were the least trustworthy ones.

Report the classified fraction with the database and the threshold beside it. A number like "62% classified" means nothing without both.
:::

> **CHECKPOINT 1** — Everyone has compared the two runs and can explain the long tail.

---

## Part 2 — Bracken, and four numbers per sample

Re-estimate genus abundances for every sample in the cohort:

```bash
for f in $D/kreports/*.kreport; do
  s=$(basename $f .kreport)
  bracken -d $DB -i $f -o bracken/$s.bracken -w bracken/$s.breport -r 150 -l G -t 10 > /dev/null
done
combine_bracken_outputs.py --files bracken/*.bracken -o genus_table.tsv
head -3 genus_table.tsv | cut -f1-5
```

`-r 150` must match the read length the Bracken files were built for, and `-l G` asks for genus: one level per run.

Now the numbers to record for every sample before any statistics:

```bash
printf "%-8s %10s %11s %9s %9s\n" sample reads classified species genera10
for f in $D/kreports/*.kreport; do
  s=$(basename $f .kreport)
  awk -F'\t' -v s=$s -v b=bracken/$s.bracken '
    BEGIN {while ((getline l < b) > 0) {split(l, a, "\t"); if (a[6]+0 >= 10) g++}}
    $4=="U" {u=$2} $4=="R" {r=$2} $4=="S" {sp+=$2}
    END {t=u+r; printf "%-8s %10d %10.1f%% %8.1f%% %9d\n", s, t, 100*r/t, 100*sp/t, g}' $f
done
```

**Question 3.** Sort that table in your head by the classified column. Does any sample stand out, and what could explain it?

::: Answer
One sample classifies at a much lower fraction than the rest. Possible causes: host reads that survived host removal, contamination from another source, or a sample from a different environment whose organisms the database does not hold.

This sample is the one to investigate before it goes into any statistics, and you would not have noticed it in its own report: the Week 2 habit of looking across samples. If it is also in one group only, it can create a difference between groups on its own.
:::

> **CHECKPOINT 2** — Everyone has `genus_table.tsv` and a per-sample table, with any outlier identified.

---

## Part 3 — Diversity, and the confound that invalidates PERMANOVA

```bash
cat > diversity.R <<'EOF'
suppressPackageStartupMessages(library(vegan))
tab    <- read.delim("genus_table.tsv", check.names = FALSE)
num    <- grep("_num$", names(tab))
counts <- t(as.matrix(tab[, num]))
colnames(counts) <- tab$name
rownames(counts) <- sub("(\\.bracken)?_num$", "", rownames(counts))
meta   <- read.delim(file.path(Sys.getenv("D"), "metadata.tsv"), row.names = 1)
meta   <- meta[rownames(counts), , drop = FALSE]
stopifnot(!anyNA(meta$group))

depth <- rowSums(counts)
rich  <- specnumber(counts)
shan  <- diversity(counts, "shannon")
print(data.frame(group = meta$group, depth, richness = rich, shannon = round(shan, 2)))
cat("\nDoes depth differ by group? Wilcoxon p =", signif(wilcox.test(depth ~ meta$group)$p.value, 2), "\n")
cat("Does richness track depth?   Spearman rho =", round(cor(depth, rich, method = "spearman"), 2), "\n\n")

d <- vegdist(counts / depth, "bray")
set.seed(1)
print(adonis2(d ~ group + batch, data = meta, permutations = 999, by = "terms"))
print(permutest(betadisper(d, meta$group), permutations = 999))
EOF
Rscript diversity.R
```

**Question 4.** Before reading the diversity numbers: does depth differ between groups, and does richness rise with depth? Why check those first?

::: Answer
Sequence a sample more deeply and you find more taxa, always, and with Kraken2 richness also counts how many database entries the reads touched. If depth differs between groups, a richness difference between groups is partly or wholly a depth difference. One test and one correlation tell you whether the comparison is fair before you make it.

If depth does track condition, that is a finding about the design before it is a finding about diversity. Compare Shannon, which leans on common taxa, with richness, which leans on rare ones: if only richness differs, suspect depth.
:::

**Question 5.** Read the PERMANOVA (`adonis2`) and the dispersion test (`permutest`) together. Which of the lecture's three outcomes do you have?

::: Answer
- **PERMANOVA significant, dispersion not:** a genuine difference in community composition. The clean result.
- **Both significant:** ambiguous. The centres moved, the spread changed, or both. Report both and say so.
- **Dispersion significant only:** one group is more variable than the other. Often the real biology; instability is a phenotype.

Also read the R² column, not only the p-value: it is the fraction of variation the term explains. Batch is in the model on purpose, and it comes before the conclusion about group, because a batch that lines up with group is a design problem no test can fix.
:::

> **CHECKPOINT 3** — Everyone has alpha diversity, a depth check, and a PERMANOVA read alongside its dispersion test.

---

## Part 4 — Differential abundance, two ways

The lecture's advice was to pick a method for stated reasons and check that your answer survives a second one. Here are two compositionally aware methods, run on the same filtered table:

```bash
cat > da.R <<'EOF'
suppressPackageStartupMessages({library(ALDEx2); library(ANCOMBC); library(TreeSummarizedExperiment)})
tab    <- read.delim("genus_table.tsv", check.names = FALSE)
num    <- grep("_num$", names(tab))
counts <- t(as.matrix(tab[, num]))
colnames(counts) <- tab$name
rownames(counts) <- sub("(\\.bracken)?_num$", "", rownames(counts))
meta   <- read.delim(file.path(Sys.getenv("D"), "metadata.tsv"), row.names = 1)
meta   <- meta[rownames(counts), , drop = FALSE]

keep <- colMeans(counts > 0) >= 0.2
cat(sum(keep), "of", length(keep), "genera are present in at least 20% of samples\n")
x <- t(round(counts[, keep]))                        # genera x samples

set.seed(1)
ald <- aldex(x, conditions = meta$group, mc.samples = 128, test = "t", effect = TRUE, verbose = FALSE)
aldex_hits <- rownames(ald)[which(ald$we.eBH < 0.05)]

tse <- TreeSummarizedExperiment(assays = list(counts = x), colData = meta)
out <- ancombc2(data = tse, assay_name = "counts", fix_formula = "group + batch",
                group = "group", p_adj_method = "holm", prv_cut = 0, lib_cut = 0,
                n_cl = 4, verbose = FALSE)
res   <- out$res
dcol  <- grep("^diff_group", names(res), value = TRUE)[1]
ancom_hits <- res$taxon[which(res[[dcol]])]

cat("\nALDEx2 hits (BH q < 0.05):", length(aldex_hits), "\n"); print(aldex_hits)
cat("\nANCOM-BC2 hits:", length(ancom_hits), "\n");          print(ancom_hits)
cat("\nFound by both:", paste(intersect(aldex_hits, ancom_hits), collapse = ", "), "\n\n")
print(round(ald[order(-abs(ald$effect)), c("effect", "we.eBH")][1:5, ], 3))
EOF
Rscript da.R 2> logs_da.txt
```

The prevalence filter comes first: a genus present in two samples out of twenty cannot support a test.

**Question 6.** Which genera do both methods call, and which only one? How would you report each?

::: Answer
Genera found by both methods, with a sizeable ALDEx2 effect, are your defensible findings. Report the effect size beside the q-value, because with enough samples a trivial difference becomes significant.

A genus significant under one method only is a lead, not a result: say that it did not survive the second method. Disagreement is expected; benchmarks find surprisingly little overlap between methods, partly because they make different choices about compositionality and zeros.

One more check before you believe any of them: is the genus in the outlier sample from Part 2 driving it? And whatever survives is still a change in proportion. To claim the amount changed you would need a spike-in, qPCR or flow cytometry.
:::

> **CHECKPOINT 4** — Everyone has two lists of differentially abundant genera and knows which they would report.

---

## Part 5 — How good are the genomes?

An nf-core/mag run on this cohort binned its assemblies into metagenome-assembled genomes, and CheckM2 scored each bin. Place them in the MIMAG tiers:

```bash
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
  {c=$h["Completeness"]; x=$h["Contamination"]
   if (x >= 10)            t["not a genome (contamination 10%+)"]++
   else if (c > 90 && x < 5) t["high-quality candidate"]++
   else if (c >= 50)       t["medium-quality draft"]++
   else                    t["low-quality draft"]++}
  END {for (k in t) printf "%-36s %d\n", k, t[k]}' $D/checkm2_quality_report.tsv
```

And the best bins, with the numbers that answer Week 3's questions:

```bash
printf "bin\tcompleteness\tcontamination\tN50\tsize\n"
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
  $h["Contamination"] < 5 {print $h["Name"]"\t"$h["Completeness"]"\t"$h["Contamination"]"\t"$h["Contig_N50"]"\t"$h["Genome_Size"]}' \
  $D/checkm2_quality_report.tsv | sort -t "$(printf '\t')" -k2,2nr | head -8
```

**Question 7.** Why are the best bins only high-quality *candidates*, and which of Week 3's three questions has no column in this table at all?

::: Answer
MIMAG's high-quality tier also requires the 23S, 16S and 5S rRNA genes and at least 18 tRNAs. CheckM2 does not check those, and short-read assemblies usually lose the rRNA operons, because they are repetitive and shared between organisms. Long-read MAGs keep them, which is one reason they reach this tier far more often.

The missing question is correctness. You have contiguity (N50, and in the full table the contig count), completeness, and contamination, which is new. There is no base-level quality, because there is no clean set of reads from one organism to measure it against.
:::

**Question 8.** The cohort profile has far more genera than there are good bins. Why, and could you use each bin's read coverage as its abundance?

::: Answer
Recovery favours abundant, simple, strain-poor organisms. Rare members, and species with several co-occurring strains, do not assemble into good bins, so a set of MAGs is not a census.

And no: a MAG's coverage is not its abundance in the community, and forty MAGs do not mean forty organisms. Use the profile from Parts 1 and 2 for abundance and the bins for genomes, and keep the two answers separate.
:::

> **CHECKPOINT 5** — Everyone can place a bin in a MIMAG tier and say what the tier cannot tell you.

---

## Takeaways

| Task | Command |
|---|---|
| Profile reads | `kraken2 --confidence 0.1 --report sample.kreport --output /dev/null` |
| Read a report | rank code in column 4: `U`, `R`, `S`; clade reads in column 2 |
| Re-estimate abundance | `bracken -r <read length> -l G -t 10`, then `combine_bracken_outputs.py` |
| Alpha and beta diversity | vegan: `specnumber`, `diversity`, `vegdist(…, "bray")` |
| Test composition | `adonis2(d ~ group + batch)` with `permutest(betadisper(d, group))` |
| Differential abundance | ALDEx2 and ANCOM-BC2 after a prevalence filter; report both |
| MAG tiers | CheckM2 completeness and contamination against MIMAG |

**Three things to carry forward:**

1. A profile is a statement about a database and a threshold. Report both, with the classified fraction.
2. Check depth and dispersion before believing a diversity result, and check a differential-abundance result with a second method.
3. A MAG is an assembly: contiguity, completeness and contamination can be measured, correctness largely cannot, and coverage is not abundance.
