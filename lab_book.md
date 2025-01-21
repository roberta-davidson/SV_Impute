# Structural Variant Imputation - Computational Lab book

Notes as I learn

## Literature

### Imputation

[**Noyvert et al. preprint.** - Imputation of structural variants using a
multi-ancestry long-read sequencing panel
enables identification of disease
associations](https://www.medrxiv.org/content/10.1101/2023.12.20.23300308v1.full.pdf)

[**Chen et al. 2021** - Investigating the Effect of Imputed Structural Variants from Whole-Genome Sequence on Genome-Wide Association and Genomic Prediction in Dairy Cattle](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7922624/#notes-a.s.ftitle)

[**Zhang et al. 2023** - The efficient phasing and imputation pipeline of low‐coverage whole genome  sequencing data using a high‐quality and publicly available reference panel in cattle](https://onlinelibrary.wiley.com/doi/epdf/10.1002/aro2.8)

### Structural variation

[**Gustafson et al. preprint.** - Nanopore sequencing of 1000 Genomes Project samples to build a comprehensive catalog of human genetic variation](https://www.medrxiv.org/content/10.1101/2024.03.05.24303792v1)

### Indigenous Australians

[**Silcocks et al. 2023** - Indigenous Australian genomes show deep structure and rich novel variation](https://www.nature.com/articles/s41586-023-06831-w)

[**Reis et al. 2023** - The landscape of genomic structural variation in Indigenous Australians](https://www.nature.com/articles/s41586-023-06842-7)

## QUILT

I chose to use QUILT because it performed the best in [this](https://paperpile.com/c/4Usx5n/NSdm) benchmarking study.

### Running the sample dataset

```bash
#modules
module load singularity
module list
## Sirngularity command to execute quilt
#QUILT.R=singularity exec -B /g/data/ r-quilt_1.0.5--r43h06b5641_0.sif QUILT.R
wd=/g/data/te53/rd8238/data

# test example from github #
rm -rf ./quilt_output 
mkdir ./quilt_output
cd ../data/
singularity exec -B /g/data/ /g/data/te53/rd8238/containters/r-quilt_1.0.5--r43h06b5641_0.sif QUILT.R \
--outputdir=quilt_output \
--chr=chr20 \
--regionStart=2000001 \
--regionEnd=2100000 \
--buffer=10000 \
--bamlist=$wd/package_2021_01_15A/bamlist.1.0.txt \
--posfile=$wd/package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.posfile.txt \
--phasefile=$wd/package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.phasefile.txt \
--reference_haplotype_file=$wd/package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.noNA12878.hap.gz \
--reference_legend_file=$wd/package_2021_01_15A/ALL.chr20_GRCh38.genotypes.20170504.chr20.2000001.2100000.noNA12878.legend.gz \
--genetic_map_file=$wd/package_2021_01_15A/CEU-chr20-final.b38.txt.gz \
--nGen=100 \
--save_prepared_reference=TRUE
```

Input files in the example are:

- text file listing input bams (bamlist)
- haplegend fileset of `.samples` `.hap` and `.legend` where:

 `.samples` looks like

```text
sample population group sex
HG00096 HG00096 HG00096 2
HG00098 HG00098 HG00098 2
```

  `.hap` looks like (row=SNP, col=reference haplotype)

```text
0 0 1 0
1 0 0 0
```

`.legend` looks like

```text
id position a0 a1
chr20:1990005_G_A 1990005 G A
chr20:1990034_T_A 1990034 T A
```

- genetic map file

```text
position COMBINED_rate.cM.Mb. Genetic_Map.cM.
82590 6.80436893751158 0
82603 6.8056503043227 0.0000884567961876505
83158 6.81470108539 0.00386559271508675
```

- position file

```text
chr20   1990005 G       A
chr20   1990034 T       A
chr20   1990057 A       G
```

- phase file

```text
NA12878HT       NA12878ONT      NA12878
0|0     0|0     0|0
0|0     0|0     0|0
```

### 1000 Genome data testing

I have structural variant calls from 1kGP:
`sniffles2_joint_sv_calls.vcf.gz`

To get haplegend fileset I run:

```bash
bcftools convert sniffles2_joint_sv_calls.vcf.gz --haplegendsample sniffles2_joint_sv_calls
```

Output files:

```bash
sniffles2_joint_sv_calls.samples
sniffles2_joint_sv_calls.legend.gz
sniffles2_joint_sv_calls.hap.gz
```

BUT, `*.hap.gz` looks like:

```bash
0* 0* 0* 0* 0* 0* 0* 0* 0* ? ? ? ? 0* 0* 0* 0* 0* 0* 0* 0* 0* 0* ? ? ? ? ? ? 0* 0* 0* 0* 0* 0* 0* 0* 0* 0* 0* 0* 0* 0* ? ? ? ? 0* 0* 0* 0* 0* 0* ? ? ? ? ? ? ? ? 0* 0*
```

Which is an issue reported when the vcf is not phased, and so haploid genotypes are called. \
Therefore the longread reference panel must be a phased VCF. \
NCIG data is not trios so cannot biologically phase, but can work on statistical phasing.

**NB** NCIG data is being base recalibrated, so wll have to feed this data into the process later.

### Phased 1kgp data

1000G data is phased SNV + SV calls [here](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)

From this I create the haplegend files again, and now the `.hap.gz` file looks like it should:
```bash
1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 1 0 0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
```

### Tri-allelic variants

Then, when trying to run quilt, it errors on multiallelic variants appearing as multiple lines with the same position.

to remove in vcf:

```bash
bcftools norm -d all in.vcf > out.vcf
```

# PIVOT!
I realised that QUILT is the wrong tool for the job as it only imputs SNPs and I misinterpreted the paper I read about that.

After reading literature i realise their is no benchmark of Imputation tools for structural variant imputation, and no tool designed intentionally for that. Although some papers do it, they do so to serve a different purpose.

Decide to benchmark tools that have been used for structural variation in the 1000 genomes dataset to adress ancestry variation, and test various coverages.

# SV Imputation Tool Benchmark

## Tools:

FImpute, Beagle and Minimac3

### Minimac3

- chromosome names must not contain "chr"
- input vcf + index files
- run by chromosome separately
- many other options can be modified

```bash
#chromosomes with 'chr'
refHap=/g/data/te53/rd8238/data/1kG/20220422_3202_phased_SNV_INDEL_SV/1kgp_SNV_INDEL_SV_chr22.vcf.gz

#chromosomes as numbers
targetVcf=/g/data/te53/rd8238/data/1kG/shortread/vcf/1kgp_chr22.vcf.gz

#basic command
singularity exec -B /g/data/ /g/data/te53/rd8238/containers/minimac3_2.0.1.sif Minimac3-omp --refHaps $refHap \
    --haps $targetVcf \
    --prefix 1kgp_chr22 \
    --chr 22 \
    --cpus 12
```

### Beagle

- input vcfs, with phased genotypes and GT field
- input genetic map in plink format

```bash
#inputs
#chromosomes as numbers, phased, with GT field
refHap=/g/data/te53/rd8238/data/1kG/20220422_3202_phased_SNV_INDEL_SV/1kgp_SNV_INDEL_SV_chr22.vcf.gz
#chromosomes as numbers, phased, with GT field
targetVcf=/g/data/te53/rd8238/data/1kG/shortread/vcf/1kgp_chr22.vcf.gz
map=/g/data/te53/rd8238/data/genetic_map/plink.chr22.GRCh38.map

#currently with defaults coded in
#singularity command to execute beagle
singularity exec -B /g/data/ /g/data/te53/rd8238/containers/beagle_5.4_22Jul22.46e--hdfd78af_0.sif beagle -Xmx10g gt=$targetVcf \
    ref=$refHap \
    out="1kgp_chr22_beagle" \
    map=$map \
    chrom=22 \
    impute=true \
    imp-states=1600 \
    cluster=0.005 \
    ap=true \
    gp=true \
    nthreads=1
```

### FImpute

- difficult to find package to prepare input files
- FImpute is depreciated anyway
- Discard from study

## calculating r^2

- want to calulate correlation between imputed and true variants
- can i simply count the proportion that are genoytped simlarly?
- want to be able to subset ths by varant type (insertion, deletion, snp, repeats, inversion)

"multiallelc" repeat SVs are recorded lke this:

```bash
22 11724390 22:11724390:C:CTAT C CTAT
22 11724390 22:11724394:C:CTATTAT C CTATTAT
22 11724390 22:11724391:C:CTATTATTAT C CTATTATTAT
22 11724390 22:11724392:C:T C T
22 11724390 22:11724393:CTAT:C CTAT C
22 11724390 22:11724395:CTATTAT:C CTATTAT C
```

**BUT**

- r^2 is recorded in vcf outputs from beagle and minimac3

# Pipeline plan

- Ensure vcfs split by chromosome
- Convert contig names to numbers only, no “chr”
- Left align all variants. + other data cleaning if necessary. (better with ref fasta but not sure how to find)

```bash
bcftools norm --multiallelics +any -o out.norm.vcf -Oz in.vcf 
bcftools index out.norm.vcf
```

- (Optional) mask known STR regions (imputing STRs is pointless waste of compute)
- (Optional) mask centromeric regions as they contain a lot of STRs.
- Impute SVs using Beagle
- Impute SVs using Minimac3
- Extract metrics from all (r2, allele length, allele frequency)
- Combine across chromosomes
- Categorize metric by variant type (Insertion, deletion etc.)
- Secret option 13: downsample short read (short variant) data to various depths and run comparatively.

### Outstanding questions

- do we care about imputing sex chromosomes?
- how do we annotate SVs by type (INSertion, DELetion, DUPlication, INVersion, BNP(breakpoint))
- Categorise by MAF 0<0.005<0.05<0.5 = rare | low frq | common
- subset by ancestry?
- Do I need to perform leave-one-out imputation to benchmark? How do i set up reference panels for this?
- Do I need to re-calculaye r2 concordance to accurately compare different tools which have their own internal r2 that may be variable?
- Also, I think I read somewhere that the structural variants recorded in the VCFs are not the full variant, but an encoding of them. So a 100bp variant may be encoded as only the beginning. Somewhere there is files to decode this, (unfortunately forgot), but it means that the variant length and type should not be calculated simply from reading the vcf but from looking up the variant in a dictionary.
