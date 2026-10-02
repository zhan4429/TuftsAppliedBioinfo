# Week 1 Session 1: Command Line, Text Processing, and Git

Applied Bioinformatics - Tufts University Department of Biology

In this lab you will log in to the Tufts HPC cluster, orient yourself in the filesystem, inspect a real GFF3 annotation file, practice common Unix text-processing patterns, and create your first Git repository for course work.

## Before You Start

You need a terminal on the Tufts HPC cluster, Pax.

Option A: OnDemand

1. Go to <https://ondemand-prod.pax.tufts.edu/>.
2. Log in with your Tufts UTLN and password, then complete 2FA.
3. Open **Clusters > Tufts HPC Shell Access**.

Option B: SSH

```bash
ssh your_utln@login-prod.pax.tufts.edu
```

You should land in your home directory, with a path like `/cluster/home/your_utln`, on a login node with a name like `login-p01`.

## Shared Setup

We can set up a few environment variables to store the course path we will use repeatedly. This is a common practice in bioinformatics pipelines.

```bash
export COURSE=/cluster/tufts/bio_appbio
export WEEK1_DATA=$COURSE/shared/wk1
ls "$WEEK1_DATA"
```

Expected file:

```text
anno.gff3
```

Do not run real software on the login node. The commands in this lab are tiny and safe; Lab2 will introduce how to request a compute node.

## Part 1: Where Am I?

### 1.1 Orient Yourself

```bash
pwd
whoami
hostname
ls -lha
```

**Question 1.** What is the full path to your home directory?

<details>
<summary>Answer</summary>
At Tufts HPC, your home directory is `/cluster/home/your_utln`.
</details>
### 1.2 Make a Working Directory

```bash
mkdir -p $COURSE/week-01/session-01
cd $COURSE/week-01/session-01
pwd
```

### 1.3 Relative and Absolute Paths

```bash
ls "$WEEK1_DATA"      # absolute path through a variable
ls ../                # relative path; depends where you are
cd ..
pwd
cd -                  # return to the previous directory
pwd
```

**Question 2.** If you are in `$COURSE/week-01/session-01`, what does `../..` refer to?

<details>
<summary>Answer</summary>

`..` is `$COURSE/week-01`, so `../..` is `$COURSE`.

</details>

### 1.4 Link the Data Instead of Copying It

```bash
cd $COURSE/week-01/session-01
ln -sf "$WEEK1_DATA/anno.gff3" .
ls -la
```

The arrow in the listing shows that `anno.gff3` is a symbolic link. You can use the file without making a second copy, which matters when the file is 40 GB instead of 40 KB.

**Checkpoint 1.** You should be in `$COURSE/week-01/session-01` with one symbolic link to `anno.gff3`.

## Part 2: Interrogating a Real Annotation File

`anno.gff3` is a genome annotation file: one feature per line, nine tab-separated columns.

### 2.1 Look Before You Leap

```bash
less anno.gff3
wc -l anno.gff3
```

Press `q` to leave `less`.

Lines starting with `#` are not data. Most are the `###` separators Ensembl writes after each gene, not headers. Compare these counts:

```bash
grep -c '' anno.gff3       # count every line
grep -c '^#' anno.gff3     # headers and ### separators — expect thousands, not five
grep -vc '^#' anno.gff3    # count data lines only
```

### 2.2 `cut`: Pull Out Columns

Column 3 is the feature type.

```bash
grep -v '^#' anno.gff3 | cut -f3 | head
```

### 2.3 The Counting Idiom

This four-stage pattern is used constantly in bioinformatics. Build it up one stage at a time and watch what each command adds.

```bash
grep -v '^#' anno.gff3 | cut -f3 | head
grep -v '^#' anno.gff3 | cut -f3 | sort | head
grep -v '^#' anno.gff3 | cut -f3 | sort | uniq -c
grep -v '^#' anno.gff3 | cut -f3 | sort | uniq -c | sort -rn
```

- cut -f3 — pull out the one column you want to count
- sort — bring identical values next to each other
- uniq -c — collapse the runs and count them
- sort -rn — put the biggest counts first

**Question 3.** Remove the `sort` before `uniq -c` and run the command again. What happens, and why?

<details>
<summary>Answer</summary>

`uniq` only collapses lines that are already adjacent. Without `sort`, identical feature types scattered through the file are not brought together, so you get a long list of small counts instead of one count per type.

It produces no error. It just quietly gives the wrong answer, which is the point of the exercise.

</details>

### 2.4 `awk`: Filter on a Column

```bash
awk -F'\t' '$3 == "gene"' anno.gff3 | wc -l
```

The `-F` option sets the field separator to a tab. Use it for GFF3 because the attributes column can contain spaces.

Count genes per chromosome:

```bash
awk -F'\t' '$3 == "gene"' anno.gff3 | cut -f1 | sort | uniq -c | sort -rn
```

### 2.5 `awk`: Do Arithmetic Across Lines

So far `awk` has been filtering lines. It can also accumulate a running total as it reads,
which is how you answer questions about quantity rather than count.

Columns 4 and 5 are the start and end coordinates of each feature. This adds up the length
of every exon in the file:

```bash
awk -F'\t' '$3 == "exon" { s += $5 - $4 + 1 } END { print s }' anno.gff3
```

Four things are happening:

```text
-F'\t'                 split each line on tabs
$3 == "exon"           only act on rows where column 3 is exon
{ s += $5 - $4 + 1 }   add this exon's length to a running total called s
END { print s }        after the last line, print the total once
```

`awk` creates `s` the first time you use it and starts it at zero, so there is nothing to
declare. The block in braces runs once per matching line; the `END` block runs once, at the
very end.

**Question 4.** Why is there a `+ 1`? Run the command without it and see how much the
answer changes.

<details>
<summary>Answer</summary>

GFF3 coordinates are **1-based and inclusive**: the first base of a sequence is base 1, and
both the start and end coordinates are part of the feature. A feature from 100 to 110
therefore covers eleven bases, not ten.

BED format is different. It is **0-based and half-open**, so the same region is written
99 to 110 and its length really is `end - start`, with no `+ 1`.

Drop the `+ 1` and your total is short by exactly one base per exon. With roughly 7,500
exons in this file that is about 7,500 bases missing — a small enough error to look
entirely plausible, which is what makes it dangerous.

The difference between the two answers is the exon count. Try it:

```bash
awk -F'\t' '$3 == "exon" { s += $5 - $4 + 1 } END { print s }' anno.gff3
awk -F'\t' '$3 == "exon" { s += $5 - $4 }     END { print s }' anno.gff3
awk -F'\t' '$3 == "exon"' anno.gff3 | wc -l
```

Mixing coordinate conventions is one of the most common silent errors in genomics, and it
never produces an error message.

</details>

### How much of this genome is transcribed?

A total number of bases is hard to interpret on its own. Divide it by the genome size and
it becomes a proportion you can reason about. The yeast genome is about 12.1 Mb:

```bash
awk -F'\t' '$3 == "exon" { s += $5 - $4 + 1 } END { printf "%.1f%% of the genome\n", 100*s/12100000 }' anno.gff3
```

`printf` works as it does in C: `%.1f` prints a number to one decimal place, and `%%`
prints a literal percent sign.

**Question 5.** You should get a strikingly high number. What would the same calculation
give for the human genome, and why?

<details>
<summary>Answer</summary>

Yeast comes out around 70%. Most of its genome is transcribed, because it has very few
introns, short intergenic regions, and almost no repetitive DNA. It is a compact genome
under selection for fast replication.

Human exons cover only a few percent of the genome. The difference is not that humans have
fewer genes but that human genes are spread
across far more space, with large introns, extensive regulatory regions, and roughly half
the genome made of repetitive elements.

This is why a gene-density figure tells you more about genome architecture than a gene
count does, and it is worth knowing before you interpret any coverage statistic.

</details>

<details>
<summary>One caveat worth knowing</summary>

This sums exon lengths; it does not measure distinct exonic positions. Where a gene has
several annotated transcripts, their shared exons are counted once per transcript, so the
total is inflated.

In yeast that barely matters, because alternative splicing is rare. In human it matters a
great deal.
</details>

### 2.6 `sed`: Make a Targeted Substitution

This file names its chromosomes `I`, `II`, `III`. Many tools, and the UCSC genome browser, expect `chrI`, `chrII`, `chrIII`. Feeding one naming convention to a tool that expects the other produces no error and no overlapping features, which is a slow way to lose an afternoon. Meet the problem now.

Look at the names first:

```bash
grep -v '^#' anno.gff3 | cut -f1 | sort -u
```

Now try the obvious fix, and look carefully at what it did:

```bash
sed 's/^/chr/' anno.gff3 | head -3
```

**Question 6.** Look at the first three lines of that output. What did you just break?

<details>
<summary>Answer</summary>

The header lines. `##gff-version 3` became `chr##gff-version 3`, which is no longer a valid GFF3 header. Tools reading the file will either fail or silently ignore it.

`^` anchors to the start of a line, and comment lines have a start too. The anchor did exactly what you asked; you asked for the wrong thing.

</details>

The fix is to exclude the comment lines:

```bash
sed '/^#/!s/^/chr/' anno.gff3 > anno_chr.gff3
head -3 anno_chr.gff3
grep -v '^#' anno_chr.gff3 | cut -f1 | sort -u
```

The script has two parts, an address and a command:

```text
/^#/        address   lines that start with #
    !       negate    ... that do NOT match
     s/^/chr/         substitute: replace the start of the line with chr
```

`sed` reads one line at a time and runs the command only where the address matches. Here
the `!` inverts it, so header lines and `###` separators pass through untouched while data
lines gain their prefix.

Replacing `^` with text looks odd, because `^` matches a position rather than a character.
That is how you prepend: find the start of the line, put something there. Using `$` instead
would append to the end.

Note the output goes to a **new** file. `sed` streams to standard output and never edits the
input, so redirecting back onto `anno.gff3` would have emptied it before `sed` read a byte of it.

Confirm two things: the headers survived, and every sequence name gained its prefix.


**Question 7.** Your `sort -u` output puts `IX` between `IV` and `Mito`, and `V` after `Mito`. Why, and does it matter here?

<details>
<summary>Answer</summary>

`sort` is comparing text, not Roman numerals. Alphabetically `IX` follows `IV`, and `V` follows `Mito`, so the order is correct as a string sort and meaningless as a chromosome order.

It does not matter for counting unique names, which is all you asked for. It would matter if you were producing a report for a reader, or joining two files on a sorted key. Know what your sort is actually sorting.

</details>

**Checkpoint 2.** You have run `grep`, `cut`, `sort`, `uniq`, `awk`, and `sed` on a real annotation file.

## Part 3: Git

Your homework will use Git repositories. Set one up now.

### 3.1 Configure Git Once Per Account

```bash
git config --global user.name "Your Name"
git config --global user.email "your.name@tufts.edu"
git config --global init.defaultBranch main
git config --list | head
```

### 3.2 Create a Repository

```bash
mkdir -p $COURSE/$USER/hw ## ToDO
cd $COURSE/$USER/hw
git init
git status
```

`git status` will become the command you type most often. It tells you what changed and what is staged for commit.

A file in a Git repository is in one of three places, and two commands move it between them:

```text
working directory  --(git add)-->  staging area  --(git commit)-->  repository
```

Keep that picture in mind; it explains the behaviour of `git diff` later in this lab.

### 3.3 Add `.gitignore` First

```bash
cat > .gitignore <<'EOF'
# Sequencing data: never commit these
*.fastq
*.fastq.gz
*.fq.gz
*.bam
*.bai
*.sam
*.vcf.gz

# Generated analysis directories
data/
results/
work/
.nextflow*

# Editor and operating system noise
.DS_Store
*~
EOF

cat .gitignore
```

Add `.gitignore` before the first commit. Once a file is committed, it is in the repository history and in every future clone. GitHub warns above 50 MB and rejects files above 100 MB; one FASTQ file can be larger than your whole repository should be.

### 3.4 Make Your First Commits

```bash
cat > README.md <<'EOF'
# Applied Bioinformatics Homework

Coursework repository for Applied Bioinformatics.
EOF

git add .gitignore README.md
git status
git commit -m "Add README and gitignore for sequencing data"
```

Make a second commit so your repository has a small history:

```bash
mkdir -p hw1
cat > hw1/notes.md <<'EOF'
# Week 1 Notes

- GFF3 is 1-based and inclusive.
- BED is 0-based and half-open.
- A command can run without being scientifically correct.
EOF

git add hw1/notes.md
git commit -m "Add Week 1 lab notes"
git log --oneline
```

Expected shape:

```text
9f3c1a2 Add Week 1 lab notes
4b7e0d8 Add README and gitignore for sequencing data
```

Your hashes will be different.

### 3.5 Prove `.gitignore` Works

```bash
touch bigfile.fastq.gz
git status
```

`bigfile.fastq.gz` should not appear as untracked. If it does, your `.gitignore` is in the wrong directory.

```bash
rm bigfile.fastq.gz
```

**Question 8.** Write a commit message for "I fixed the awk command that was counting exons wrong." What makes it better than `fixed bug`?

<details>
<summary>Answer</summary>

Something like: `Fix exon length calculation for 1-based inclusive GFF3 coordinates`

It says what changed and why. In six months, `git log` may be the only record of your reasoning, and `fixed bug` tells your future self almost nothing.

</details>

**Checkpoint 3.** `git log --oneline` shows two commits, and a `.fastq.gz` file is ignored.

## Takeaways

### Commands from this lab

Count how many times each value appears in a column:

```bash
grep -v '^#' anno.gff3 | cut -f3 | sort | uniq -c | sort -rn
```

Keep only the rows where a column has a particular value:

```bash
awk -F'\t' '$3 == "gene"' anno.gff3
```

Add up a calculated value across every matching row:

```bash
awk -F'\t' '$3 == "exon" { s += $5 - $4 + 1 } END { print s }' anno.gff3
```

Substitute text at the start of every line, skipping comments:

```bash
sed '/^#/!s/^/chr/' anno.gff3 > anno_chr.gff3
```

Git, in the order you will use it:

```bash
git status              # what has changed, and what is staged
git diff                # exactly what changed, line by line
git add FILE            # stage a change
git commit -m "..."     # record it
git log --oneline       # what has happened so far
```

### Two habits that will save you time all semester

**Check the coordinate convention before doing arithmetic.** GFF3 is 1-based and inclusive, so a feature's length is `end - start + 1`. BED is 0-based and half-open, so the same region's length is `end - start`. Neither format announces which it is, and mixing them gives answers that are wrong by exactly one base per feature.

**Write `.gitignore` before your first commit, not after.** Once a file is committed it is in the history permanently and in every clone. Adding it to `.gitignore` afterwards stops future tracking but does not remove it.

### The idea underneath both

Every command in this lab ran successfully. None of them checked whether the answer made sense. `uniq` without `sort` returns a wrong count with no error. Forgetting `+ 1` returns a total that is off by the number of features. `sed` without the comment filter quietly corrupts your header lines.

Deciding whether an answer is plausible is your job, not the tool's. That is most of what this course is about, and it starts here.
