# Lab 4B — Calling Variants Against Your Own Genome

**Applied Bioinformatics · Week 4, Session 2**
**Runs alongside the lecture**

> The reads and the reference come from the same haploid yeast. A perfect assembly and a perfect caller would find nothing. Everything you find today is information about one of them.

> **Read Part 1 first and submit both callers immediately.** They run while you do the rest of the lab.

---

## Before you start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname

export COURSE=/cluster/tufts/yzhang85/bioinfo_course/appliedBio
export MYWORK="$COURSE/$USER"
export SHARED=/cluster/tufts/bio_appbio/shared/wk4
export GENOME=12100000

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

## Part 1 — Submit both callers, now

### 1.1 The files every caller expects

```bash
samtools faidx asm.fa
gatk CreateSequenceDictionary -R asm.fa -O asm.dict 2> logs/dict.log
ls asm.fa.fai asm.dict
```

### 1.2 GATK HaplotypeCaller

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

export COURSE=/cluster/tufts/yzhang85/bioinfo_course/appliedBio
export MYWORK="$COURSE/$USER"
module load gatk
set -euo pipefail
cd $MYWORK/week-04/session-02

gatk --java-options "-Xmx12g" HaplotypeCaller \
  -R asm.fa -I sample.bam -O gatk.vcf.gz \
  --native-pair-hmm-threads ${SLURM_CPUS_PER_TASK}
EOF
```

### 1.3 DeepVariant

```bash
cat > dv.sh <<'EOF'
#!/bin/bash
#SBATCH --job-name=deepvariant
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
module load deepvariant
export COURSE=/cluster/tufts/yzhang85/bioinfo_course/appliedBio
export MYWORK="$COURSE/$USER"
cd $MYWORK/week-04/session-02

run_deepvariant \
    --model_type=WGS \
    --ref=asm.fa \
    --reads=sample.bam \
    --output_vcf=dv.vcf.gz \
    --num_shards=${SLURM_CPUS_PER_TASK} \
    --intermediate_results_dir=$SCRATCH/dv_$SLURM_JOB_ID
EOF
```

### 1.4 Submit both

```bash
bash -n gatk.sh && bash -n dv.sh
GATK_JOB=$(sbatch --parsable gatk.sh)
DV_JOB=$(sbatch --parsable dv.sh)
printf '%s\n%s\n' "$GATK_JOB" "$DV_JOB" > caller_job_ids.txt
squeue -j "$GATK_JOB,$DV_JOB"
```

> Check back at Part 3.

**Question 1.** Both callers will run with their default assumption that the sample is diploid. S288C is haploid. Is that a mistake?

::: Answer
It is a deliberate choice, and worth understanding. A diploid model can report a pileup as mixed (0/1). In a haploid genome no site is truly mixed, so every 0/1 call becomes a flag: reads from different places were stacked on one position.

A haploid model (`--sample-ploidy 1` in GATK) must force every site to 0 or 1. That is the correct biology, but it hides the mixed pileups by assigning them to one allele. Part 2 shows the difference. DeepVariant offers only a diploid germline model for this genome, which is one more reason the model is part of the result.
:::

> **CHECKPOINT 1** — Both jobs are queued and their IDs saved.

---

## Part 2 — Read a VCF while they run

Call 500 kb of the largest contig interactively, once with each ploidy:

```bash
LARGEST=$(sort -k2,2nr asm.fa.fai | head -1 | cut -f1)
echo "$LARGEST"

gatk HaplotypeCaller -R asm.fa -I sample.bam -L ${LARGEST}:1-500000 \
  -O peek_p2.vcf.gz 2> logs/peek_p2.log
gatk HaplotypeCaller -R asm.fa -I sample.bam -L ${LARGEST}:1-500000 \
  --sample-ploidy 1 -O peek_p1.vcf.gz 2> logs/peek_p1.log
```

### 2.1 Header and records

```bash
bcftools view -h peek_p2.vcf.gz | grep -E '^##(contig|FORMAT)' | cut -c1-80 | head -8
bcftools view -H peek_p2.vcf.gz | head -3 | cut -c1-150
```

The `##contig` lines repeat your assembly's names and lengths: like a BAM's `@SQ`, they bind this VCF to one exact reference.

### 2.2 Genotypes, under each model

```bash
echo "diploid model:";  bcftools query -f '[%GT]\n' peek_p2.vcf.gz | sort | uniq -c
echo "haploid model:";  bcftools query -f '[%GT]\n' peek_p1.vcf.gz | sort | uniq -c
```

### 2.3 Look at the mixed sites

```bash
bcftools query -i 'GT="het"' -f '%POS\t%REF\t%ALT\t%QUAL\t[%AD]\t%INFO/DP\t%INFO/MQ\n' peek_p2.vcf.gz | head
```

Columns: position, reference allele, alternative allele, QUAL, reads supporting each allele, depth, mapping quality.

**Question 2.** Look at the 0/1 calls. What do their allele depths, depth and MQ have in common, and what happens to the same positions under the haploid model?

::: Answer
Typical patterns: allele depths far from a clean split (for example 30 reference reads against 8 alternative) and depth well above the contig's usual level. Those are the signature of mismapping: reads from another copy of a repeat were placed here, and their differences look like a second allele.

MQ depends on the kind of repeat. Where the assembly kept several partly different copies, MQ drops below 60. Where it collapsed many copies into one, as with the rDNA, every read maps uniquely to the one copy that exists, MQ stays at 60, and the excess depth is the giveaway. (Reads with MAPQ 0 are ignored by the caller, so the fully resolved Ty copies from Tuesday produce silence, not 0/1 calls.)

Under the haploid model the same sites are called 1 (alternative) when the alternative reads win, or are not reported at all when the reference reads win. The evidence of mixing disappears from the output. Same reads, different model, a different and cleaner-looking answer: the lecture's first takeaway.
:::

**Question 3.** Now look at any 1/1 call in the diploid file. Every read disagrees with the reference at that position. In this genome, what does that mean?

::: Answer
The assembly is wrong at that base. The reads and the assembly come from the same DNA, so when every read says something else, the reads are right and the assembly carries an error, most often a one-base insertion or deletion in a homopolymer.

A 1/1 call against your own genome is not a variant. It is a correction.
:::

> **CHECKPOINT 2** — Everyone can read a VCF record and say what a 0/1 and a 1/1 call mean in this genome.

---

## Part 3 — What the whole genome says about your assembly

```bash
squeue -j "$(paste -sd, caller_job_ids.txt)"
ls -lh gatk.vcf.gz dv.vcf.gz
```

If either job is still running, use the staged calls, which were made on the staged set: `cp $SHARED/fallback/{gatk,dv}.vcf.gz* .` (and use the staged `asm.fa` too).

### 3.1 Normalise, then filter

```bash
for v in gatk dv; do
  bcftools norm -f asm.fa -m -any $v.vcf.gz -Oz -o $v.norm.vcf.gz 2> logs/$v.norm.log
done

# GATK leaves FILTER empty: apply the lecture's hard filters. DeepVariant already filters.
bcftools filter -s hardfail -e 'QD<2 || FS>60 || MQ<40 || SOR>3' gatk.norm.vcf.gz \
  | bcftools view -f PASS -Oz -o gatk.pass.vcf.gz
bcftools view -f PASS dv.norm.vcf.gz -Oz -o dv.pass.vcf.gz
for v in gatk dv; do bcftools index -f $v.pass.vcf.gz; done
```

### 3.2 Count by meaning

```bash
for v in gatk dv; do
  printf "%-5s  1/1 SNVs %5d   1/1 indels %5d   0/1 %5d\n" $v \
    $(bcftools view -H -i 'GT="AA"' -v snps   $v.pass.vcf.gz | wc -l) \
    $(bcftools view -H -i 'GT="AA"' -v indels $v.pass.vcf.gz | wc -l) \
    $(bcftools view -H -i 'GT="het"'          $v.pass.vcf.gz | wc -l)
done
```

### 3.3 Turn the errors into a quality score

Use the 1/1 calls both callers agree on:

```bash
bcftools isec -n=2 -w1 -i 'GT="AA"' -i 'GT="AA"' gatk.pass.vcf.gz dv.pass.vcf.gz -Oz -o errors.vcf.gz
N=$(bcftools view -H errors.vcf.gz | wc -l)
awk -v e=$N -v g=$GENOME 'BEGIN {printf "%d errors in %d bases: QV %.1f\n", e, g, -10*log(e/g)/log(10)}'
```

**Question 4.** Compare that QV with the Merqury QV from Homework 3. Why is agreement between them meaningful, and why might they still differ?

::: Answer
They are independent instruments. Merqury compares k-mers in the assembly with k-mers in the BGISEQ reads; this compares aligned reads with the assembly, base by base. Agreement from two different methods is far stronger evidence than either number alone, the same logic as Lab 3A's two signals.

They differ because each sees different things. Variant calling is blind in the MAPQ 0 islands from Tuesday, so errors inside repeats are missed and this QV is optimistic. Merqury counts every assembly k-mer the reads lack, including some in regions the reads cover poorly, so it can be pessimistic. If there are no shared errors at all, the QV is undefined: say that the error count is below what this many reads can detect.
:::

### 3.4 Where the 0/1 calls are

```bash
bcftools query -i 'GT="het"' -f '%CHROM\n' gatk.pass.vcf.gz | sort | uniq -c | sort -rn | head
bcftools query -i 'GT="het"' -f '%INFO/DP\n' gatk.pass.vcf.gz | sort -n | \
  awk '{a[NR]=$1} END {print "median depth at 0/1 calls:", a[int((NR+1)/2)]}'
awk '$1=="total" {print "genome-wide mean depth:", $4}' ../session-01/all.mosdepth.summary.txt
```

**Question 5.** Compare depth at the 0/1 calls with the genome-wide depth, and look at which contigs hold most of them. What is the evidence that they are mismapping rather than biology?

::: Answer
Depth at 0/1 calls is usually well above the genome-wide mean, and the calls cluster on a few contigs, above all the one carrying the collapsed rDNA, where about 150 copies' worth of reads land on the few copies that were assembled. Extra depth means reads from several copies landed on one; clustering means the problem is local to repeats, not spread across the genome as real variation would be.

Note what does *not* produce them: the Ty elements that made Tuesday's MAPQ 0 islands. The caller ignores MAPQ 0 reads, so those windows are silent. Silence is not the same as no errors: an assembly error inside a Ty copy cannot be called at all.

And the decisive argument: S288C is haploid. There is no second haplotype for a real 0/1 call to come from.
:::

### 3.5 One sanity check from the lecture

```bash
bcftools stats gatk.pass.vcf.gz | awk '$1=="TSTV" {print "Ti/Tv:", $5}'
```

**Question 6.** That is the transition/transversion ratio of these calls. Why is it not a measure of mutation here, and what would a value near 0.5 tell you?

::: Answer
These calls are not mutations: they are assembly errors and mismapping artefacts, so their Ti/Tv describes how errors arise, not how genomes evolve.

A value near 0.5 is what random substitutions give, because there are twice as many possible transversions as transitions. It would mean the substitutions in this set behave like noise, which is exactly what you would expect from errors.
:::

> **CHECKPOINT 3** — Everyone has counts, a QV estimate and an explanation for where the 0/1 calls sit.

---

## Part 4 — Which errors broke a gene?

The Week 3 annotation was predicted on this assembly, errors included. First check that it really is the same sequence:

```bash
cmp -s asm.fa $MYWORK/week-03/session-02/asm.fa && echo "same assembly" || echo "DIFFERENT: use the staged set"
```

Annotate the agreed errors against your BRAKER3 gene models. VEP needs to know which transcripts code for protein, and BRAKER3 does not say, so add that attribute while preparing the file:

```bash
grep -v '^#' $MYWORK/week-03/session-02/braker/braker.gff3 \
  | awk -F'\t' 'BEGIN {OFS="\t"} $3=="mRNA" || $3=="transcript" {sub(/;$/, "", $9); $9 = $9 ";biotype=protein_coding"} {print}' \
  | sort -k1,1 -k4,4n -k5,5n | bgzip -c > braker.gff3.gz
tabix -p gff braker.gff3.gz

module load vep
vep \
  -i errors.vcf.gz --gff braker.gff3.gz --fasta asm.fa \
  --tab --fields "Location,Allele,Feature,Consequence" \
  -o errors.vep.tsv --force_overwrite --no_stats 2> logs/vep.log

grep -v '^#' errors.vep.tsv | cut -f4 | tr ',' '\n' | sort | uniq -c | sort -rn
```

Staged-set users: replace the BRAKER path with `$SHARED/fallback/braker.gff3`.

**Question 7.** Suppose some errors are reported as `frameshift_variant`. The alternative allele is the correct sequence, so what has actually gone wrong, and in which file?

::: Answer
The assembly carries an indel inside a gene, and BRAKER3 predicted its gene model around that error. "Frameshift" here means that correcting the assembly would shift the reading frame of the model, so the model is the thing that is wrong: probably truncated, or joined to its neighbour.

That is Week 3's warning made concrete: functional annotation inherits structural errors, and structural annotation inherits assembly errors. One base in the assembly, a wrong protein in every file downstream. The fix belongs at the top: polish the assembly with these calls, then annotate again.
:::

> **CHECKPOINT 4** — Everyone can trace one assembly error into a gene model.

---

## Takeaways

| Task | Command |
|---|---|
| Reference files for callers | `samtools faidx`, `gatk CreateSequenceDictionary` |
| Call variants | `gatk HaplotypeCaller`, DeepVariant in a container |
| Change the ploidy model | `--sample-ploidy 1` |
| Normalise before comparing | `bcftools norm -f asm.fa -m -any` |
| Hard filters | `bcftools filter -s hardfail -e 'QD<2 \|\| FS>60 \|\| MQ<40 \|\| SOR>3'` |
| Calls two callers share | `bcftools isec -n=2 -w1` |
| Consequences on your gene models | `vep --gff genes.gff3.gz --fasta asm.fa` |

**Three things to carry forward:**

1. A variant call is a model's reading of a pileup. A diploid model on a haploid genome turns mismapped repeats into 0/1 calls; a haploid model hides them.
2. Against your own genome, a 1/1 call is a correction and a 0/1 call is a mapping problem. Counted, they measure the assembly as well as Merqury does.
3. Errors travel downstream. An indel in the assembly becomes a broken gene model, then a wrong protein.
