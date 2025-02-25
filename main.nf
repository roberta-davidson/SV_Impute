#!/usr/bin/env nextflow

//nextflow.enable.dsl=1
nextflow.enable.dsl=2

params.ref_vcf = 'ref.vcf'
params.target_vcf = 'target.vcf'
params.output_dir = 'output'
params.beagle_map = 'plink.map'
params.fasta = 'ref.fa'

process checkInputs {
    input:
    path ref_vcf
    
    output:
    stdout

    script:
    """
    echo "Hello World!"
    echo "Input ref VCF: ${params.ref_vcf}"
   #echo "Input target VCF: ${params.target_vcf}"
    """
}

process REFremoveChrPrefix {
    // to publish output files to output directory
    publishDir "${params.output_dir}/rename_chr_ref_vcf/", mode: 'copy'

    input:
    each path(ref_vcf)

    output:
    path("no_chr*.vcf.gz")

    script:
    """
	bcftools annotate --rename-chrs /g/data/te53/rd8238/SVimp/dev/chr_rename.txt $ref_vcf --threads 4 -Oz -o no_chr_$ref_vcf
    """
}

process REFleftAlignAndCollapse {
    publishDir "${params.output_dir}/left_align_ref_vcf/", mode: 'copy'
    errorStrategy 'ignore'

    input:
    each path(ref_cleaned_vcf)
    path fasta

    output:
    path("alignvar_*.vcf.gz")

    script:
    """
    bcftools annotate --set-id '%CHROM\\_%POS' "$ref_cleaned_vcf" --threads 4 -Oz -o pos_"$ref_cleaned_vcf"
    bcftools norm --fasta-ref $fasta --multiallelics -any pos_"$ref_cleaned_vcf" --threads 4 -Oz -o alignvar_"$ref_cleaned_vcf"
    bcftools index -t alignvar_"$ref_cleaned_vcf"
    """
}
//

process REFcheckAndSplitVCF {
    // to publish output files to output directory
    publishDir "${params.output_dir}/split_ref/", mode: 'copy'

    input:
    each path(ref_alignvar_vcf)

    output:
    path("ref_split_*.vcf.gz")

    script:
    """
    bcftools index "$ref_alignvar_vcf"
            #make variable list of all chromosomes in vcf and separate only those
    chrs=\$(zcat "$ref_alignvar_vcf" | grep -v "#" | awk '{print \$1}' | uniq)
    for chr in \$chrs; do
        bcftools view -r "\${chr}" --threads 4 -Oz -o ref_split_"\${chr}".vcf.gz "$ref_alignvar_vcf"
    done
    """
}

process TARGETremoveChrPrefix {
    publishDir "${params.output_dir}/rename_chr_target_vcf/", mode: 'copy'

    input:
    each path(target_vcf)

    output:
    path("no_chr_*.vcf.gz")

    script:
    """
	bcftools annotate --rename-chrs /g/data/te53/rd8238/SVimp/dev/chr_rename.txt $target_vcf --threads 4 -Oz -o no_chr_$target_vcf
    """
}

process TARGETleftAlignAndCollapse {
    publishDir "${params.output_dir}/left_align_target_vcf/", mode: 'copy'
    input:
    each path(target_cleaned_vcf)
    path fasta

    output:
    path("alignvar_*.vcf.gz")

    script:
    """
    bcftools annotate --set-id '%CHROM\\_%POS' $target_cleaned_vcf --threads 4 -Oz -o pos_$target_cleaned_vcf
    bcftools norm --fasta-ref $fasta --multiallelics -any pos_$target_cleaned_vcf --threads 4 -Oz -o alignvar_$target_cleaned_vcf 
    bcftools index -t alignvar_$target_cleaned_vcf
    """
}

process TARGETcheckAndSplitVCF {
    // to publish output files to output directory
    publishDir "${params.output_dir}/split_target/", mode: 'copy'

    input:
    each path(target_alignvar_vcf)

    output:
    path("target_split_*.vcf.gz")

    script:
    """
    bcftools index "$target_alignvar_vcf"
    chrs=\$(zcat "$target_alignvar_vcf" | grep -v "#" | awk '{print \$1}' | uniq)
    for chr in \$chrs; do
        bcftools view -r "\${chr}" --threads 4 -Oz -o target_split_"\${chr}".vcf.gz "$target_alignvar_vcf"
    done
    """
}

process IMPUTE_Mimimac3 {
    // to publish output files to output directory
    publishDir "${params.output_dir}/minimac3/", mode: 'copy'

    input:
    tuple val(chr), path(map), path(ref), path(target)
    
    output:
    tuple val(chr), path("1kgp_minimac3_*.vcf.gz")

    script:
    """
    Minimac3-omp --refHaps ${ref} --haps ${target} \
    --prefix 1kgp_minimac3_$chr --chr $chr --cpus 24
    """
}

process IMPUTE_Beagle {
    // to publish output files to output directory
    publishDir "${params.output_dir}/beagle/", mode: 'copy'

    input:
    tuple val(chr), path(map), path(ref), path(target)

    output:
    tuple val(chr), path("1kgp_beagle_*.vcf.gz"), emit: beagle_vcf

    script:
    """
    beagle -Xmx40g gt=${target} ref=${ref} out="1kgp_beagle_${chr}" \
    map=${map} \
    chrom=${chr} impute=true imp-states=1600 cluster=0.005 ap=true gp=true nthreads=12
    """
}

process BeaglecollectMetrics {
    publishDir "${params.output_dir}/metrics/", mode: 'copy'
    
    input:
    tuple val(chr), path(beagle_vcf)

    output:
    path("beagle_metric_*.txt")

    script:
    """
    bcftools query -f '%CHROM %POS %REF %ALT %INFO/DR2 %INFO/AF %INFO/IMP\\n' ${beagle_vcf} -o beagle_metric_${chr}.txt
    """
}

process Minimac3collectMetrics {
    publishDir "${params.output_dir}/metrics/", mode: 'copy'
    
    input:
    tuple val(chr), path(minimac3_vcf)

    output:
    tuple val(chr), path("minimac3_metric_*.txt")

    script:
    """
    bcftools query -f '%CHROM %POS %REF %ALT %INFO/R2 %INFO/AF %INFO/MAF %INFO/ER2\\n' ${minimac3_vcf} -o minimac3_metric_${chr}.txt
    """
}

//main workflow
workflow {
	//files
    ref_vcf = file(params.ref_vcf)
	target_vcf = file(params.target_vcf)
    beagle_map_files = channel.fromPath(params.beagle_map)
    fasta = file(params.fasta)
    //chromosomes
    chr = Channel.of(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
	//prep ref
    ref_cleaned_vcf = REFremoveChrPrefix(ref_vcf)
    ref_alignvar_vcf = REFleftAlignAndCollapse(ref_cleaned_vcf,fasta)
    ref_split_vcf_files = REFcheckAndSplitVCF(ref_alignvar_vcf)
    //prep target
    target_cleaned_vcf = TARGETremoveChrPrefix(target_vcf)
    target_alignvar_vcf = TARGETleftAlignAndCollapse(target_cleaned_vcf,fasta)
    target_split_vcf_files = TARGETcheckAndSplitVCF(target_alignvar_vcf)
    // channel combine for imputation
    //beagle_map_files
    beagle_map_files
        .map {
            file ->
                file = file
                chr = file.getName().replaceAll(/plink.chr(\d+)\.GRCh38.map/, '$1')
                [ chr, file ]
        }
        .set{beagle_by_chr}
    //ref
    ref_split_vcf_files
        .map {
            file ->
                file = file
                chr = file.getName().replaceAll(/ref_split_(\d+)\.vcf.gz$/, '$1')
                [ chr, file ]
        }
        .set{refs_by_chr}
    //target
    target_split_vcf_files
        .map {
            file ->
                file = file
                chr = file.getName().replaceAll(/target_split_(\d+)\.vcf.gz$/, '$1')
                [ chr, file ]
        }
        .set{target_by_chr}

    beagle_by_chr.join(refs_by_chr).join(target_by_chr).set{merged}
    merged.view()

    //Impute
   // minimac3_vcf = IMPUTE_Mimimac3(merged)
   // beagle_vcf = IMPUTE_Beagle(merged)
    IMPUTE_Mimimac3(merged).set{minimac3_vcf}
    IMPUTE_Beagle(merged).set{beagle_vcf}
    System.out.println("beagle_vcf is :")
    System.out.println(beagle_vcf)
    
    //Get metrics
    beagle_metrics = BeaglecollectMetrics(beagle_vcf) 
    minimac3_metrics = Minimac3collectMetrics(minimac3_vcf)
