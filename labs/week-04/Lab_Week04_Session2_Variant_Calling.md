# Week 4 Session 2: Variant Calling

Applied Bioinformatics - Tufts University Department of Biology

> S288C is a haploid lab strain: each position in its nuclear genome has one allele. Today you call variants with a haploid model, so every call marks a place where the BGISEQ reads disagree with your HiFi assembly, and you investigate why.

> **Read Part 1 first and submit the caller immediately.** It runs while you do the rest of the lab.

---

## Before You Start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED=$COURSE/shared/wk4

mkdir -p $MYWORK/week-04/session-02/logs && cd $MYWORK/week-04/session-02

# Tuesday's reference and alignment
S1=$MYWORK/week-04/session-01
cp $S1/asm.fa .
ln -sf $S1/sample.bam sample.bam
ln -sf $S1/sample.bam.bai sample.bam.bai
ls -lhL asm.fa sample.bam sample.bam.bai

module load samtools bcftools htslib gatk
module list
```

If Tuesday's alignment is missing, use the staged set. Its files belong together, so take all of them:

```bash
cp $SHARED/fallback/asm.fa .
ln -sf $SHARED/fallback/sample.bam sample.bam
ln -sf $SHARED/fallback/sample.bam.bai sample.bam.bai
```

---

### Dataset Identity Matters

The [BGISEQ BioSample](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN24477750) and
[HiFi BioSample](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN24477612) both identify
**S288C**, but they are separate samples, collected almost a year apart. A difference
between the reads and the assembly could be an assembly error, a real difference
between the two cultures, or a mapping artifact. Treat every call as a candidate to
investigate, not as a correction to make.

## Part 1: Submit The Caller, Now

### 1.1 The Files Every Caller Expects

```bash
samtools faidx asm.fa
gatk CreateSequenceDictionary -R asm.fa -O asm.dict 2> logs/dict.log
ls asm.fa.fai asm.dict
```

`asm.fa.fai` lets tools jump to any position in the assembly; `asm.dict` lists each contig's name and length for GATK.

### 1.2 GATK HaplotypeCaller, Haploid

`--sample-ploidy 1` tells HaplotypeCaller that the genome has one copy of each chromosome. Without it, GATK assumes two copies, as in a human.

```bash
cat > gatk.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=gatk_hc
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
module load gatk
set -euo pipefail
cd $MYWORK/week-04/session-02
mkdir -p "tmp/gatk_$SLURM_JOB_ID"

gatk --java-options "-Xmx12g -Djava.io.tmpdir=$PWD/tmp/gatk_$SLURM_JOB_ID" HaplotypeCaller \
  -R asm.fa -I sample.bam -O calls.vcf.gz \
  --sample-ploidy 1 \
  --native-pair-hmm-threads ${SLURM_CPUS_PER_TASK}
EOF

bash -n gatk.sh
GATK_JOB=$(sbatch --parsable gatk.sh)
echo "$GATK_JOB" > caller_job_id.txt
squeue -j "$GATK_JOB"
```

> Check back at Part 3.

**Question 1.** What does `--sample-ploidy 1` change about how HaplotypeCaller reports a site, and what does it not change?

<details>
<summary>Answer</summary>

With one copy per chromosome, each site has a single allele, so a call is reported as `1` (the alternative allele) instead of choosing between `0/1` and `1/1` as a diploid model would. The setting changes the genotype model only. The reads are the same, and reads that disagree with each other are still there; you will see them in the allele depths in Part 2.

</details>

**Checkpoint 1.** The job is queued and its ID saved.

---

## Part 2: Read A VCF While It Runs

Call 500 kb of the largest contig interactively, with the same haploid setting:

```bash
LARGEST=$(sort -k2,2nr asm.fa.fai | head -1 | cut -f1)
echo "$LARGEST"
gatk HaplotypeCaller -R asm.fa -I sample.bam -L ${LARGEST}:1-500000 \
  --sample-ploidy 1 -O peek.vcf.gz 2> logs/peek.log
```

### 2.1 Header And Records

```bash
bcftools view -h peek.vcf.gz | grep -E '^##(contig|FORMAT)' | cut -c1-80 | head -8
bcftools view -H peek.vcf.gz | head -3 | cut -c1-150
```

The `##contig` lines repeat your assembly's names and lengths: like a BAM's `@SQ`, they bind this VCF to one exact reference.

### 2.2 Genotypes And Read Support

```bash
bcftools query -f '[%GT]\n' peek.vcf.gz | sort | uniq -c
bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t[%AD]\t%INFO/DP\t%INFO/MQ\n' peek.vcf.gz | head
```

Columns: position, reference allele, alternative allele, QUAL, allele depths (reads supporting the reference, then the alternative), total depth, mapping quality.

The **alternative allele fraction** is the share of reads that support the alternative allele: alternative reads / (reference reads + alternative reads). For example, `AD` of `2,38` gives 38 / 40 = 0.95. This prints it for every call:

```bash
bcftools query -f '%POS\t[%AD]\n' peek.vcf.gz \
  | awk -F'[\t,]' '$2+$3 > 0 {printf "%s\t%.2f\n", $1, $3/($2+$3)}'
```

**Question 2.** In a haploid genome, what alternative allele fraction do you expect at a real difference between the reads and the assembly? What could a call with a fraction of only 0.3 to 0.6 mean? Look at its depth and MQ too.

<details>
<summary>Answer</summary>

Close to 1. With one copy of each chromosome, all reads from a real difference should carry the same allele, apart from a few sequencing errors.

A fraction of 0.3 to 0.6 means the reads at that position disagree with each other, which a single haploid copy cannot explain. The usual cause is mismapping: reads from another copy of a repeat were placed here. Such sites often have depth well above the contig's usual level. MQ is lower where the assembly kept several similar copies, but stays at 60 where it collapsed many copies into one, as with the rDNA, because every read then maps uniquely to the one copy that exists. Contamination or a mixed culture are other possibilities. The haploid model still reports these sites as `1`, so the allele depths are where the mixing shows.

</details>

**Checkpoint 2.** Everyone can read a VCF record and calculate an alternative allele fraction.

---

## Part 3: What The Whole Genome Says About Your Assembly

```bash
squeue -j "$(cat caller_job_id.txt)"
ls -lh calls.vcf.gz
```

If the job is still running, use the staged calls, which were made on the staged set: `cp $SHARED/fallback/calls.vcf.gz* .` (and use the staged `asm.fa` too).

### 3.1 Filter

HaplotypeCaller reports every candidate and leaves the FILTER column empty. These hard filters remove calls with weak or suspicious evidence:

| Filter | Removes calls with |
|---|---|
| `QD<2` | low quality for the number of reads covering the site |
| `FS>60`, `SOR>3` | evidence that comes mostly from reads on one strand |
| `MQ<40` | reads that the aligner could not place confidently |

```bash
bcftools filter -s hardfail -e 'QD<2 || FS>60 || MQ<40 || SOR>3' calls.vcf.gz \
  | bcftools view -f PASS -Oz -o calls.pass.vcf.gz
bcftools index -f calls.pass.vcf.gz
```

### 3.2 Count

```bash
for f in calls.vcf.gz calls.pass.vcf.gz; do
  printf "%-18s total %5d   SNPs %5d   indels %5d\n" $f \
    $(bcftools view -H $f | wc -l) \
    $(bcftools view -H -v snps $f | wc -l) \
    $(bcftools view -H -v indels $f | wc -l)
done
```

### 3.3 Clean And Mixed Calls

Split the filtered calls into **clean** calls (alternative allele fraction 0.8 or more) and **mixed** calls (below 0.8), and compare their depth with the genome-wide mean from Tuesday:

```bash
bcftools query -f '%INFO/DP\t[%AD]\n' calls.pass.vcf.gz \
  | awk -F'[\t,]' '$2+$3 > 0 {
      g = ($3/($2+$3) >= 0.8) ? "clean" : "mixed"; n[g]++; d[g] += $1 }
    END {for (g in n) printf "%s calls: %d, mean depth %.0f\n", g, n[g], d[g]/n[g]}'
awk '$1=="total" {print "genome-wide mean depth:", $4}' ../session-01/all.mosdepth.summary.txt
```

And see which contigs hold the calls:

```bash
bcftools query -f '%CHROM\n' calls.pass.vcf.gz | sort | uniq -c | sort -rn | head
```

**Question 3.** How does the mean depth of the mixed calls compare with the clean calls and with the genome-wide mean? Are the calls spread across the genome or concentrated on a few contigs? What does that suggest about where the mixed calls come from?

<details>
<summary>Answer</summary>

Mixed calls usually sit at depth well above the genome-wide mean and cluster on a few contigs or regions: the signature of repeats, where reads from several copies pile onto one. Clean calls usually sit near the genome-wide depth. Check the pattern in your own numbers rather than assuming it; a mixed call at normal depth would need another explanation, such as contamination or a mixed culture.

</details>

**Question 4.** A clean call has normal depth, good MQ and nearly all reads supporting the alternative allele. Give two explanations, and say what evidence would tell them apart.

<details>
<summary>Answer</summary>

It may be an error in the HiFi assembly, or a real difference between the two separately collected S288C cultures. Look at the original HiFi reads at that position: if they also support the alternative allele, the assembly is wrong at that base; if they agree with the assembly, the two cultures probably differ. Do not change the assembly on the BGISEQ evidence alone.

</details>

**Checkpoint 3.** Everyone has counts before and after filtering, and an explanation for the mixed calls.

---

## Going Further (Optional)

For those who finish early:

- **A second caller.** DeepVariant also supports haploid calling through `--haploid_contigs`, which takes the names of this assembly's contigs. Comparing two callers requires normalising both VCFs first (`bcftools norm -f asm.fa -m -any`) and then intersecting them (`bcftools isec`).
- **Which calls touch a gene?** Annotate the filtered calls against your Week 3 BRAKER3 gene models. VEP needs to know which transcripts code for protein, and BRAKER3 does not say, so add that attribute while preparing the file:

```bash
cmp -s asm.fa $MYWORK/week-03/session-02/asm.fa && echo "same assembly" || echo "DIFFERENT: use the staged set"
grep -v '^#' $MYWORK/week-03/session-02/braker/braker.gff3 \
  | awk -F'\t' 'BEGIN {OFS="\t"} $3=="mRNA" || $3=="transcript" {sub(/;$/, "", $9); $9 = $9 ";biotype=protein_coding"} {print}' \
  | sort -k1,1 -k4,4n -k5,5n | bgzip -c > braker.gff3.gz
tabix -p gff braker.gff3.gz
module load vep
vep \
  -i calls.pass.vcf.gz --gff braker.gff3.gz --fasta asm.fa \
  --tab --fields "Location,Allele,Feature,Consequence" \
  -o calls.vep.tsv --force_overwrite --no_stats 2> logs/vep.log
grep -v '^#' calls.vep.tsv | cut -f4 | tr ',' '\n' | sort | uniq -c | sort -rn
```

Staged-set users: replace the BRAKER path with `$SHARED/fallback/braker.gff3`. A label such as `frameshift_variant` describes a consequence if the difference is real; it is not proof that the assembly is wrong.

---

## Takeaways

| Task | Command |
|---|---|
| Reference files for callers | `samtools faidx`, `gatk CreateSequenceDictionary` |
| Call variants, haploid | `gatk HaplotypeCaller --sample-ploidy 1` |
| Read genotypes and allele depths | `bcftools query -f '[%GT]\t[%AD]\n'` |
| Hard filters | `bcftools filter -s hardfail -e 'QD<2 \|\| FS>60 \|\| MQ<40 \|\| SOR>3'` |
| Count SNPs and indels | `bcftools view -H -v snps`, `-v indels` |

**Three things to carry forward:**

1. A genotype is a model's interpretation of the reads. Set the ploidy to match the organism.
2. In a haploid genome, a real difference has an alternative allele fraction close to 1; mixed read support points to repeats or mixed samples.
3. A clean call is a candidate: it can be an assembly error or a real difference between samples, and only more evidence, such as the original HiFi reads, tells them apart.

## Software Reference

- [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller)
- [DeepVariant haploid support](https://github.com/google/deepvariant/blob/r1.6/docs/deepvariant-haploid-support.md)
