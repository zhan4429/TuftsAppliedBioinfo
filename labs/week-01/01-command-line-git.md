# Week 1 Session 1: Command Line, Text Processing, and Git

Applied Bioinformatics - Tufts University Department of Biology

In this lab you will log in to the Pax cluster, orient yourself in the filesystem, inspect a real GFF3 annotation file, practice common Unix text-processing patterns, and create your first Git repository for course work.

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

Set the course data path once. If your instructor gives a different course path, change only the first line.

```bash
export COURSE=/cluster/tufts/bio_appbio
export WEEK1_DATA=$COURSE/shared/wk1
ls -hl "$WEEK1_DATA"
```

Expected file:

```text
anno.gff3
```

Do not run real software on the login node. The commands in this lab are tiny and safe; Session 2 covers how to request a compute node.

## Part 1: Where Am I?

### 1.1 Orient Yourself

```bash
pwd
whoami
hostname
ls -lha
```

**Question 1.** What is the full path to your home directory?

### 1.2 Make a Working Directory

```bash
mkdir -p ~/appbio/week-01/session-01
cd ~/appbio/week-01/session-01
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

**Question 2.** If you are in `~/appbio/week-01/session-01`, what does `../..` refer to?

<details>
<summary>Answer</summary>

`..` is `~/appbio/week-01`, so `../..` is `~/appbio`.

</details>

### 1.4 Link the Data Instead of Copying It

```bash
cd ~/appbio/week-01/session-01
ln -sf "$WEEK1_DATA/anno.gff3" .
ls -la
```

The `->` in the listing shows that `anno.gff3` is a symbolic link. You can use the file without making a second copy, which matters when the file is 40 GB instead of 40 KB.

**Checkpoint 1.** You should be in `~/appbio/week-01/session-01` with one symbolic link to `anno.gff3`.

## Part 2: Interrogating a Real Annotation File

`anno.gff3` is a genome annotation file: one feature per line, nine tab-separated columns.

### 2.1 Look Before You Leap

```bash
less anno.gff3
wc -l anno.gff3
```

Press `q` to leave `less`.

Lines starting with `#` are header comments, not data. Compare these counts:

```bash
grep -c '' anno.gff3       # every line
grep -c '^#' anno.gff3     # header lines
grep -vc '^#' anno.gff3    # data lines only
```

**Question 3.** How many header lines does the file have?

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

**Question 4.** Remove the `sort` before `uniq -c` and run the command again. What happens, and why?

<details>
<summary>Answer</summary>

`uniq` only collapses lines that are already adjacent. Without `sort`, identical feature types scattered through the file are not brought together, so you get a long list of small counts instead of one count per type.

It produces no error. It just quietly gives the wrong answer, which is the point of the exercise.

</details>

### 2.4 `awk`: Filter on a Column

```bash
awk -F'\t' '$3 == "gene"' anno.gff3 | wc -l
```

`-F'\t'` sets the field separator to a tab. Use it for GFF3 because the attributes column can contain spaces.

Count genes per sequence:

```bash
awk -F'\t' '$3 == "gene"' anno.gff3 | cut -f1 | sort | uniq -c | sort -rn | head -5
```

### 2.5 `awk`: Do Arithmetic Across Lines

Columns 4 and 5 are start and end.

```bash
awk -F'\t' '$3 == "exon" { s += $5 - $4 + 1 } END { print s }' anno.gff3
```

**Question 5.** Why is there a `+ 1`? Try the command without it and see how much the answer changes.

<details>
<summary>Answer</summary>

GFF3 coordinates are 1-based and inclusive. A feature from 100 to 110 covers eleven bases, not ten.

BED format is different: it is 0-based and half-open. Mixing coordinate conventions is one of the most common silent errors in genomics.

</details>

### 2.6 `sed`: Make a Targeted Substitution

Some tools want `chr1`, while others want `1`. Coordinate naming mismatches can silently break downstream analysis.

```bash
cut -f1 anno.gff3 | grep -v '^#' | sort -u | head
sed 's/^chr//' anno.gff3 > anno_nochr.gff3
cut -f1 anno_nochr.gff3 | grep -v '^#' | sort -u | head
```

**Question 6.** Why does the pattern start with `^`? What would `sed 's/chr//'` do to a line containing the word `chromosome` in its description column?

<details>
<summary>Answer</summary>

`^` anchors the match to the start of the line. Without it, `sed` removes the first `chr` anywhere on the line, so `chromosome` in a description could become `omosome`.

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
mkdir -p ~/appbio/hw
cd ~/appbio/hw
git init
git status
```

`git status` will become the command you type most often. It tells you what changed and what is staged for commit.

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

Add `.gitignore` before the first commit. Once a file is committed, it is in the repository history and every future clone. GitHub warns above 50 MB and rejects files above 100 MB; one FASTQ file can be larger than your whole repository should be.

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
9f3c1a Add Week 1 lab notes
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

**Question 7.** Write a commit message for "I fixed the awk command that was counting exons wrong." What makes it better than `fixed bug`?

<details>
<summary>Answer</summary>

Something like: `Fix exon length calculation for 1-based inclusive GFF3 coordinates`

It says what changed and why. In six months, `git log` may be the only record of your reasoning, and `fixed bug` tells your future self almost nothing.

</details>

**Checkpoint 3.** `git log --oneline` shows two commits, and a `.fastq.gz` file is ignored.

## If You Finish Early

1. Find the longest gene in the annotation:

   ```bash
   awk -F'\t' '$3=="gene" { print $5-$4+1, $1, $4, $5 }' anno.gff3 | sort -rn | head -1
   ```

   Work out what each stage does.

2. Count how many genes are on each strand in column 7.
3. Write a one-liner that produces a tab-separated table of feature type and count.
4. Try `git diff`: edit `README.md`, run `git diff`, then run `git add README.md` and try `git diff --staged`.

## Takeaways

| Skill | Command |
| --- | --- |
| Count occurrences | `... \| sort \| uniq -c \| sort -rn` |
| Filter a column | `awk -F'\t' '$3=="gene"'` |
| Sum across lines | `awk '{s += $5 - $4 + 1} END {print s}'` |
| Substitute at start of line | `sed 's/^chr//'` |
| See what changed | `git status`, `git diff` |
| Record a change | `git add`, `git commit -m "..."` |

Two habits will save you time all semester:

1. Check coordinate conventions before doing arithmetic.
2. Add `.gitignore` before the first commit, not after.

## Homework Connection

Homework 1 extends the command-line and Git habits from this lab. Keep `~/appbio/week-01/session-01` and `~/appbio/hw` unless the instructor tells you to remove them.
