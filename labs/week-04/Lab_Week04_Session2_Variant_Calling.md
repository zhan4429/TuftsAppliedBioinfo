# Week 4 Session 2: Variant Calling

Applied Bioinformatics - Tufts University Department of Biology

> Today you use a haploid model for the S288C nuclear genome and investigate places where the BGISEQ reads disagree with your HiFi assembly. The model is an analysis assumption; read support and sample identity still matter.

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

module load samtools bcftools htslib gatk mosdepth
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

A haploid genotype has one allele number: `0` means reference, `1` the first alternative, and `2` a second alternative if the record has one. In this lab we select sites with exactly one alternative allele, so a non-reference genotype is `1`. Ploidy changes the genotype model and can change call confidence and which variants are reported. It does not change the reads: disagreement among them remains visible in the allele depths.

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

bcftools view -m 2 -M 2 -v snps peek.vcf.gz -Oz -o peek.snps.vcf.gz
```

We focus on **biallelic SNPs**: one reference base and one alternative base.
`-m 2 -M 2` keeps records with exactly two alleles in total; `-v snps` keeps SNPs.
This makes the two-number allele-depth calculation below unambiguous.

If this interval has no SNPs, use the worked allele-depth examples below and
return to the whole-genome calls in Part 3. An empty interval alone does not mean
the caller failed. An instructor-prepared peek file from the same reference and
haploid model can also be used to save class time.

### 2.1 Header And Records

```bash
bcftools view -h peek.snps.vcf.gz | grep -E '^##(contig|FORMAT)' | cut -c1-80 | head -8
bcftools view -H peek.snps.vcf.gz | head -3 | cut -c1-150
```

The `##contig` lines list the assembly's contig names and lengths, like a BAM's `@SQ`. These help check compatibility, but matching names and lengths do not prove that two FASTA files contain identical sequences. Keep the actual reference used for calling.

### 2.2 Genotypes And Read Support

```bash
bcftools query -f '[%GT]\n' peek.snps.vcf.gz | sort | uniq -c
bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t[%AD]\t%INFO/DP\t%INFO/MQ\n' peek.snps.vcf.gz | head
```

Columns: position, reference allele, alternative allele, QUAL, allele depths (reads supporting the reference, then the alternative), total depth, mapping quality.

The **alternative allele fraction** is the share of reads that support the alternative allele: alternative reads / (reference reads + alternative reads). For example, `AD` of `2,38` gives 38 / 40 = 0.95. This prints it for every SNP in the biallelic subset:

```bash
bcftools query -f '%POS\t[%AD]\n' peek.snps.vcf.gz \
  | awk -F'[\t,]' '$2+$3 > 0 {printf "%s\t%.2f\n", $1, $3/($2+$3)}'
```

**Question 2.** In a haploid genome, what alternative allele fraction do you expect at a real difference between the reads and the assembly? What could a call with a fraction of only 0.3 to 0.6 mean? Look at its depth and MQ too.

<details>
<summary>Answer</summary>

Close to 1 at a single-copy nuclear locus in a homogeneous haploid sample: almost all informative reads should support the same allele, apart from sequencing errors.

A fraction of 0.3 to 0.6 shows mixed read support. One possible explanation is
mismapping: reads from different repeat copies were placed at the same locus.
Elevated depth can support that explanation. MQ may be low where several similar
copies remain in the assembly, but can stay high when the assembly has collapsed
them into one apparent location; it is not guaranteed to be 60.

Contamination, a mixed culture, or an incorrect ploidy assumption are other
possibilities. Organelle contigs also need separate interpretation. Haploid calling
may omit mixed sites or report an alternative allele; the genotype alone does not
describe the mixture.

</details>

**Checkpoint 2.** Everyone can read a VCF record and calculate an alternative allele fraction.

---

## Part 3: What The Whole Genome Says About Your Assembly

```bash
sacct -j "$(cat caller_job_id.txt)" --format=JobID,State,ExitCode
ls -lh calls.vcf.gz calls.vcf.gz.tbi
```

Use your own output only after the job reports `COMPLETED` with exit code `0:0`.
A file can exist while it is still being written.

If the job is still queued, running, or failed, use the matched staged example in a
**separate directory**. The staged VCF must have been generated with the same
haploid setting. These commands preserve your running job's reference and outputs:

```bash
mkdir -p "$MYWORK/week-04/session-02/example"
cd "$MYWORK/week-04/session-02/example"
cp "$SHARED/fallback/asm.fa" .
ln -sf "$SHARED/fallback/sample.bam" sample.bam
ln -sf "$SHARED/fallback/sample.bam.bai" sample.bam.bai
cp "$SHARED/fallback/calls.vcf.gz" .
```

Continue the remaining commands in whichever directory you selected. The BAM,
FASTA, and VCF must come from the same analysis. Record whether you used your own
results or the staged example.

### 3.1 Filter

HaplotypeCaller already uses a calling-confidence threshold. Its `FILTER` value
`.` means that site filters have not been applied, not that every candidate passed.

First select biallelic SNPs from the full callset:

```bash
bcftools view -m 2 -M 2 -v snps calls.vcf.gz -Oz -o calls.snps.vcf.gz
```

The following are **illustrative SNP hard filters** for this exercise, not validated
yeast-specific cutoffs. Indels need separate assessment and are outside the core
exercise. A passing SNP still needs interpretation.

| Filter | Removes calls with |
|---|---|
| `QD<2` | low quality for the number of reads covering the site |
| `FS>60`, `SOR>3` | evidence that comes mostly from reads on one strand |
| `MQ<40` | reads that the aligner could not place confidently |

```bash
bcftools filter -s hardfail -e 'QD<2 || FS>60 || MQ<40 || SOR>3' calls.snps.vcf.gz \
  | bcftools view -f PASS -Oz -o calls.pass.snps.vcf.gz
bcftools index -f calls.pass.snps.vcf.gz
```

### 3.2 Count

```bash
for f in calls.snps.vcf.gz calls.pass.snps.vcf.gz; do
  echo "$f"
  bcftools view -H "$f" | wc -l
done
```

### 3.3 Clean And Mixed Calls

For this summary, label SNPs with alternative allele fraction **0.8 or more**
"clean" and the rest "mixed". These are teaching labels for read support; "clean"
does not prove that a call is true or that the assembly is wrong.

Save the group for each SNP so the depth and contig summaries use the same groups:

```bash
bcftools query -f '%CHROM\t%POS\t%INFO/DP\t[%AD]\n' calls.pass.snps.vcf.gz \
  | awk -F'[\t,]' 'BEGIN {OFS="\t"; print "contig", "position", "depth", "group"}
    $4+$5 > 0 {
      g = ($5/($4+$5) >= 0.8) ? "clean" : "mixed"
      print $1, $2, $3, g
    }' > call_groups.tsv
head call_groups.tsv

awk -F'\t' 'BEGIN {n["clean"]=0; n["mixed"]=0}
  NR>1 {n[$4]++; d[$4]+=$3}
  END {for (g in n) {
    if (n[g]>0) printf "%s calls: %d, mean depth %.0f\n", g, n[g], d[g]/n[g]
    else printf "%s calls: 0, mean depth not available\n", g
  }}' call_groups.tsv
```

Calculate a coverage baseline from the **matching BAM in your current directory**.
This also works when you switched to the staged example:

```bash
mosdepth -t 4 --fast-mode call_depth sample.bam
awk '$1=="total" {print "genome-wide mean depth:", $4}' call_depth.mosdepth.summary.txt
```

VCF depth and mosdepth use different read-selection rules, so treat this as context,
not an exact numerical equality. Also consider the typical depth of the relevant
nuclear contig rather than comparing every contig with one genome-wide mean.

Now count **only mixed SNPs** by contig:

```bash
awk -F'\t' 'NR>1 && $4=="mixed" {print $1}' call_groups.tsv \
  | sort | uniq -c | sort -rn | head
```

If there are no mixed SNPs, the last command prints nothing. Record zero; do not
infer a depth or a clustering pattern for a group that is absent.

**Question 3.** If mixed SNPs are present, how does their mean depth compare with the clean SNPs and the coverage baseline? Are the mixed SNPs concentrated on a few contigs? What explanations fit those observations? If there are none, what can you conclude?

<details>
<summary>Answer</summary>

High depth and concentration in repeat-rich regions would support reads from
multiple copies piling onto one assembled locus. High-copy organelle or plasmid
contigs need their own baseline. Mixed support at ordinary nuclear depth could
instead motivate checks of sample mixture, ploidy, or alignment problems.

Report the pattern you actually observe. If the mixed group is absent, say that no
SNPs in this filtered subset met that definition. Haploid calling and filtering can
remove sites with mixed evidence, so zero does not prove that the genome has no
repeats, mixtures, or mapping problems.

</details>

**Question 4.** A clean call has normal depth, good MQ and nearly all reads supporting the alternative allele. Give two explanations, and say what evidence would tell them apart.

<details>
<summary>Answer</summary>

It may be an error in the HiFi assembly, or a real difference between the two separately collected S288C cultures. Consistent, well-mapped original HiFi support for the alternative allele would favor an assembly error; HiFi support for the reference would favor a difference between cultures. Check for mapping artifacts and mixed evidence before treating either pattern as proof. Do not change the assembly on the BGISEQ evidence alone.

</details>

**Checkpoint 3.** Counts before and after SNP filtering are recorded, with an interpretation of the observed support groups, including an absent group if applicable.

---

## Going Further (Optional)

For those who finish early:

- **A second caller.** DeepVariant also supports haploid calling through `--haploid_contigs`, which takes the names of this assembly's contigs. Comparing two callers requires normalising both VCFs first (`bcftools norm -f asm.fa -m -any`) and then intersecting them (`bcftools isec`).
- **Which calls touch a gene?** Annotate the filtered calls against your Week 3 BRAKER3 gene models. VEP needs to know which transcripts code for protein, and BRAKER3 does not say, so add that attribute while preparing the file:

First select the annotation and the exact FASTA it describes. For your own Week 3
annotation:

```bash
ANNOTATION_ASM="$MYWORK/week-03/session-02/asm.fa"
GFF="$MYWORK/week-03/session-02/braker/braker.gff3"
```

If you used the staged callset, use these assignments instead:

```bash
ANNOTATION_ASM="$SHARED/fallback/asm.fa"
GFF="$SHARED/fallback/braker.gff3"
```

The following block runs annotation only if the annotation FASTA matches the
calling reference and the GFF exists:

```bash
if [ -s "$GFF" ] && cmp -s asm.fa "$ANNOTATION_ASM"; then
  samtools faidx asm.fa
  grep -v '^#' "$GFF" \
    | awk -F'\t' 'BEGIN {OFS="\t"} $3=="mRNA" || $3=="transcript" {sub(/;$/, "", $9); $9 = $9 ";biotype=protein_coding"} {print}' \
    | sort -k1,1 -k4,4n -k5,5n | bgzip -c > braker.gff3.gz
  tabix -f -p gff braker.gff3.gz
  module load vep
  vep \
    -i calls.pass.snps.vcf.gz --gff braker.gff3.gz --fasta asm.fa \
    --tab --fields "Location,Allele,Feature,Consequence" \
    -o calls.vep.tsv --force_overwrite --no_stats 2> vep.log
  grep -v '^#' calls.vep.tsv | cut -f4 | tr ',' '\n' | sort | uniq -c | sort -rn
else
  echo "Stop: the annotation is missing or its FASTA differs from the calling reference."
fi
```

The core subset contains SNPs, so expect consequences such as synonymous or
missense changes. Exploring indel consequences such as frameshifts requires a
separately assessed indel callset. A consequence label is not proof that the
assembly is wrong.

---

## Takeaways

| Task | Command |
|---|---|
| Reference files for callers | `samtools faidx`, `gatk CreateSequenceDictionary` |
| Call variants, haploid | `gatk HaplotypeCaller --sample-ploidy 1` |
| Read genotypes and allele depths | `bcftools query -f '[%GT]\t[%AD]\n'` |
| Apply SNP filters | `bcftools filter`, then `bcftools view -f PASS` |
| Select biallelic SNPs | `bcftools view -m 2 -M 2 -v snps` |

**Three things to carry forward:**

1. A genotype is a model's interpretation of the reads. Set the ploidy to match the organism.
2. A real difference at a single-copy haploid nuclear locus should have alternative support close to 1. Mixed support needs investigation; its absence from a filtered callset does not rule out mapping problems.
3. A clean call is a candidate: it can be an assembly error or a real difference between samples, and additional evidence, such as well-mapped original HiFi reads, helps distinguish those explanations.

## Software Reference

- [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller)
- [GATK hard-filtering guidance](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)
- [DeepVariant haploid support](https://github.com/google/deepvariant/blob/r1.6/docs/deepvariant-haploid-support.md)
