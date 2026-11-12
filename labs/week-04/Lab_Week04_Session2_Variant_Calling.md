# Week 4 Session 2: Variant Calling

Applied Bioinformatics - Tufts University Department of Biology

> The reads and the assembly are labeled S288C but come from separate BioSamples. Calls can reflect assembly errors, biological differences between cultures, or mapping and calling artifacts. Today you investigate the evidence for those explanations.

> **Read Part 1 first and submit both callers immediately.** They run while you do the rest of the lab.

---

## Before You Start

```bash
srun --partition=batch --nodes=1 --ntasks=1 --cpus-per-task=8 \
  --mem=32G --time=03:00:00 --pty bash
hostname

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
export SHARED=$COURSE/shared/wk4
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

### Dataset Identity Matters

The [BGISEQ BioSample](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN24477750) and
[HiFi BioSample](https://www.ebi.ac.uk/ena/browser/api/xml/SAMN24477612) both identify
**S288C**, with collection dates 2021-01-21 and 2021-12-29 respectively. These are
separate BioSamples; the records do not establish one shared DNA preparation or
prove identical genomes. We use a haploid nuclear-genome model as a teaching
assumption and interpret disagreements as candidates to investigate, rather than
automatically correcting the assembly.

## Part 1: Submit Both Callers, Now

### 1.1 The Files Every Caller Expects

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

export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
module load gatk
set -euo pipefail
cd $MYWORK/week-04/session-02
mkdir -p "tmp/gatk_$SLURM_JOB_ID"

gatk --java-options "-Xmx12g -Djava.io.tmpdir=$PWD/tmp/gatk_$SLURM_JOB_ID" HaplotypeCaller \
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

module load deepvariant
set -euo pipefail
export COURSE=/cluster/tufts/c27sp1bio0196
export MYWORK="$COURSE/$USER"
cd $MYWORK/week-04/session-02
mkdir -p "tmp/dv_$SLURM_JOB_ID"

run_deepvariant \
    --model_type=WGS \
    --ref=asm.fa \
    --reads=sample.bam \
    --output_vcf=dv.vcf.gz \
    --num_shards=${SLURM_CPUS_PER_TASK} \
    --intermediate_results_dir="$PWD/tmp/dv_$SLURM_JOB_ID"
EOF
```

### 1.4 Submit Both

```bash
bash -n gatk.sh && bash -n dv.sh
GATK_JOB=$(sbatch --parsable gatk.sh)
DV_JOB=$(sbatch --parsable dv.sh)
printf '%s\n%s\n' "$GATK_JOB" "$DV_JOB" > caller_job_ids.txt
squeue -j "$GATK_JOB,$DV_JOB"
```

> Check back at Part 3.

**Question 1.** Both scripts use diploid genotype output to let us inspect mixed read evidence. Why might that be useful for an assumed haploid nuclear genome, and what does it not prove?

<details>
<summary>Answer</summary>

A diploid model can report mixed evidence as 0/1. In a haploid nuclear genome, that is a reason to investigate mismapped repeats, contamination, mixed cells or an incorrect ploidy assumption; it is not proof of any single cause.

GATK can use `--sample-ploidy 1`. DeepVariant 1.6 and later also support `--haploid_contigs` for named contigs; this exercise deliberately leaves that option unset. Those names must match this assembly, not generic human chromosome names. DeepVariant adjusts genotype probabilities for listed haploid contigs during postprocessing, so it is incorrect to describe it as supporting only diploid output.

A haploid setting changes the genotype model. It does not remove contradictory reads or establish which sequence is biologically correct.

</details>

**Checkpoint 1.** Both jobs are queued and their IDs saved.

---

## Part 2: Read A VCF While They Run

Call 500 kb of the largest contig interactively, once with each ploidy:

```bash
LARGEST=$(sort -k2,2nr asm.fa.fai | head -1 | cut -f1)
echo "$LARGEST"

gatk HaplotypeCaller -R asm.fa -I sample.bam -L ${LARGEST}:1-500000 \
  -O peek_p2.vcf.gz 2> logs/peek_p2.log
gatk HaplotypeCaller -R asm.fa -I sample.bam -L ${LARGEST}:1-500000 \
  --sample-ploidy 1 -O peek_p1.vcf.gz 2> logs/peek_p1.log
```

### 2.1 Header And Records

```bash
bcftools view -h peek_p2.vcf.gz | grep -E '^##(contig|FORMAT)' | cut -c1-80 | head -8
bcftools view -H peek_p2.vcf.gz | head -3 | cut -c1-150
```

The `##contig` lines repeat your assembly's names and lengths: like a BAM's `@SQ`, they bind this VCF to one exact reference.

### 2.2 Genotypes, Under Each Model

```bash
echo "diploid model:";  bcftools query -f '[%GT]\n' peek_p2.vcf.gz | sort | uniq -c
echo "haploid model:";  bcftools query -f '[%GT]\n' peek_p1.vcf.gz | sort | uniq -c
```

### 2.3 Look At The Mixed Sites

```bash
bcftools query -i 'GT="het"' -f '%POS\t%REF\t%ALT\t%QUAL\t[%AD]\t%INFO/DP\t%INFO/MQ\n' peek_p2.vcf.gz | head
```

Columns: position, reference allele, alternative allele, QUAL, reads supporting each allele, depth, mapping quality.

**Question 2.** Look at the 0/1 calls. What do their allele depths, depth and MQ have in common, and what happens to the same positions under the haploid model?

<details>
<summary>Answer</summary>

Typical patterns: allele depths far from a clean split (for example 30 reference reads against 8 alternative) and depth well above the contig's usual level. Those are the signature of mismapping: reads from another copy of a repeat were placed here, and their differences look like a second allele.

MQ depends on the kind of repeat. Where the assembly kept several partly different copies, MQ drops below 60. Where it collapsed many copies into one, as with the rDNA, every read maps uniquely to the one copy that exists, MQ stays at 60, and the excess depth is the giveaway. (Reads with MAPQ 0 are ignored by the caller, so the fully resolved Ty copies from Tuesday produce silence, not 0/1 calls.)

Under the haploid model the same sites are called 1 (alternative) when the alternative reads win, or are not reported at all when the reference reads win. The evidence of mixing disappears from the output. Same reads, different model, a different and cleaner-looking answer: the lecture's first takeaway.

</details>

**Question 3.** Inspect a 1/1 call and its allele depths. What are two possible explanations when most well-mapped BGISEQ reads support an allele different from the HiFi assembly?

<details>
<summary>Answer</summary>

One possibility is an assembly error; another is a real difference between the separately collected S288C cultures. Mapping and calling artifacts also remain possible. A 1/1 genotype is a model's conclusion, not a statement that literally every read supports the alternative allele.

Check the BGISEQ and HiFi alignments, depth, strand and mapping evidence, and local repeats. Concordant support in the original HiFi reads strengthens an assembly-error explanation. Do not automatically polish a HiFi assembly to the sequence of a different BioSample.

</details>

**Checkpoint 2.** Everyone can read a VCF record and say what a 0/1 and a 1/1 call mean in this genome.

---

## Part 3: What The Whole Genome Says About Your Assembly

```bash
squeue -j "$(paste -sd, caller_job_ids.txt)"
ls -lh gatk.vcf.gz dv.vcf.gz
```

If either job is still running, use the staged calls, which were made on the staged set: `cp $SHARED/fallback/{gatk,dv}.vcf.gz* .` (and use the staged `asm.fa` too).

### 3.1 Normalise, Then Filter

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

### 3.2 Count By Meaning

```bash
for v in gatk dv; do
  printf "%-5s  1/1 SNVs %5d   1/1 indels %5d   0/1 %5d\n" $v \
    $(bcftools view -H -i 'GT="AA"' -v snps   $v.pass.vcf.gz | wc -l) \
    $(bcftools view -H -i 'GT="AA"' -v indels $v.pass.vcf.gz | wc -l) \
    $(bcftools view -H -i 'GT="het"'          $v.pass.vcf.gz | wc -l)
done
```

### 3.3 Summarize Shared Differences

Collect homozygous-alternative calls both callers share, then summarize shared SNVs:

```bash
bcftools isec -n=2 -w1 -i 'GT="AA"' -i 'GT="AA"' gatk.pass.vcf.gz dv.pass.vcf.gz -Oz -o candidate_differences.vcf.gz
N=$(bcftools view -H -v snps candidate_differences.vcf.gz | wc -l)
L=$(awk '{s+=$2} END {print s}' asm.fa.fai)
awk -v n="$N" -v g="$L" 'BEGIN {if(n>0) printf "%d shared SNVs in %d assembly bases; disagreement score %.1f\n", n,g,-10*log(n/g)/log(10); else print "No shared SNVs detected; no finite score estimated"}'
```

This is a whole-assembly **SNV disagreement score**, not a validated assembly QV.
It excludes indels, misses differences outside callable regions, and can include
biological differences between the BioSamples. Counting an indel as one error would
also confuse variant events with incorrect bases. The two callers use the same reads,
so agreement between callers does not make the evidence independent.

**Question 4.** Compare this disagreement score with the Merqury QV from Homework 3.
Why must you check sample matching, callable regions and what each method counts
before interpreting their agreement or difference?

<details>
<summary>Answer</summary>

Merqury evaluates assembly k-mer support; this exercise counts shared called SNVs.
They have different detection limits and denominators. Both use the BGISEQ library,
so they also share its coverage biases and any genuine differences from the HiFi
culture. Agreement can be informative but is not independent confirmation of an
assembly-error rate. Neither number proves that every unsupported base is an error.

</details>

### 3.4 Where The 0/1 Calls Are

```bash
bcftools query -i 'GT="het"' -f '%CHROM\n' gatk.pass.vcf.gz | sort | uniq -c | sort -rn | head
bcftools query -i 'GT="het"' -f '%INFO/DP\n' gatk.pass.vcf.gz | sort -n | \
  awk '{a[NR]=$1} END {print "median depth at 0/1 calls:", a[int((NR+1)/2)]}'
awk '$1=="total" {print "genome-wide mean depth:", $4}' ../session-01/all.mosdepth.summary.txt
```

**Question 5.** Compare depth at the 0/1 calls with the genome-wide depth, and look at which contigs hold most of them. Do these observations support a mismapping explanation, and what alternatives remain?

<details>
<summary>Answer</summary>

Enrichment at unusually high depth or in repeated regions supports a mismapping or collapsed-repeat explanation. Compare allele depths and mapping quality as well as the contig distribution; the pattern must be observed, not assumed.

Under the haploid nuclear-genome assumption, 0/1 is unexpected, but that assumption alone does not distinguish repeat mapping from contamination, mixed cultures or different ploidy. MAPQ filters can hide some repeat-associated problems; absence of calls there does not validate the assembly.

</details>

### 3.5 One Sanity Check From The Lecture

```bash
bcftools stats gatk.pass.vcf.gz | awk '$1=="TSTV" {print "Ti/Tv:", $5}'
```

**Question 6.** That is the transition/transversion ratio of these calls. Why is it not a measure of mutation here, and what would a value near 0.5 tell you?

<details>
<summary>Answer</summary>

This is the Ti/Tv of called differences between the reads and assembly. It may mix biological differences with assembly and mapping errors, so it cannot be read directly as a mutation process or an assembly-error spectrum.

For equally likely single-base substitutions, there are twice as many possible transversions as transitions, giving 0.5. A value near 0.5 is compatible with such a simple noise model, but is not proof that the calls are artifacts. Examine call types, sequence context, mapping evidence and sample provenance.

</details>

**Checkpoint 3.** Everyone has counts, a disagreement score and an explanation for where the 0/1 calls sit.

---

## Part 4: Which Differences Could Affect A Gene?

The Week 3 annotation was predicted on this assembly, errors included. First check that it really is the same sequence:

```bash
cmp -s asm.fa $MYWORK/week-03/session-02/asm.fa && echo "same assembly" || echo "DIFFERENT: use the staged set"
```

Annotate the shared candidate differences against your BRAKER3 gene models. VEP needs to know which transcripts code for protein, and BRAKER3 does not say, so add that attribute while preparing the file:

```bash
grep -v '^#' $MYWORK/week-03/session-02/braker/braker.gff3 \
  | awk -F'\t' 'BEGIN {OFS="\t"} $3=="mRNA" || $3=="transcript" {sub(/;$/, "", $9); $9 = $9 ";biotype=protein_coding"} {print}' \
  | sort -k1,1 -k4,4n -k5,5n | bgzip -c > braker.gff3.gz
tabix -p gff braker.gff3.gz

module load vep
vep \
  -i candidate_differences.vcf.gz --gff braker.gff3.gz --fasta asm.fa \
  --tab --fields "Location,Allele,Feature,Consequence" \
  -o errors.vep.tsv --force_overwrite --no_stats 2> logs/vep.log

grep -v '^#' errors.vep.tsv | cut -f4 | tr ',' '\n' | sort | uniq -c | sort -rn
```

Staged-set users: replace the BRAKER path with `$SHARED/fallback/braker.gff3`.

**Question 7.** A candidate difference is labeled `frameshift_variant`. What would it mean if further evidence confirmed an assembly error, and what if it were a real difference between the two cultures?

<details>
<summary>Answer</summary>

If independent evidence confirms an assembly indel error, the gene model may have been predicted around the wrong sequence and could need rebuilding after correction. If it is a true difference between cultures, the consequence describes a difference relative to this assembly's gene model, not proof that the assembly is wrong.

Check the transcript and protein evidence, the original HiFi reads and the local alignment before changing the assembly. A predicted frameshift is a consequence annotation, not an instruction to replace the reference allele.

</details>

**Checkpoint 4.** Everyone can trace one candidate difference to a gene model and state what evidence is needed to classify it.

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

1. A genotype is a model-based interpretation of read evidence. Inspect the reads and check the ploidy assumption.
2. Separate BioSamples can differ biologically. Shared calls and Merqury QV do not by themselves establish a true assembly-error rate.
3. A verified assembly error can affect downstream gene models; a consequence label alone does not verify that error.

## Software Reference

- [DeepVariant haploid support](https://github.com/google/deepvariant/blob/r1.6/docs/deepvariant-haploid-support.md)
