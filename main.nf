#!/usr/bin/env nextflow

params.reference    = "/home/alexandr/Documents/SNP/data/MT.fna"
params.reads        = "/home/alexandr/Documents/SNP/data/*R{1,2}.fastq"
params.outdir       = "results"
params.help = false

// in your_script.nf
if ( params.help ) {
    help = """SNP.nf: A description of your script and maybe some examples of how
             |                to run the script
             |Required arguments:
             |  --input_file  Location of the input file file.
             |                [default: ${params.reads}]
             |
             |Optional arguments:
             |  --use_thing   Do some optional process.
             |  -w            The NextFlow work directory. Delete the directory once the process
             |                is finished [default: ${workDir}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

log.info """\
    ===================================
     S N P   P I P E L I N E
    ===================================
    reference: ${params.reference}
    reads    : ${params.reads}
    outdir   : ${params.outdir}
    """
    .stripIndent(true)

    
/*
 * Define the `REFINDEX` process that creates an index
 * given the reference genome file
 */
process REFINDEX {
    tag "$reference"
    input:
    path reference

    output:
    path "*"

    script:
    """
    bwa index $reference
    """
}

/*
 * Define the `QCONTROL` process that performs quality trimming and filtering of reads
 */
process QCONTROL{
    tag "${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/fastp"

    input:
    tuple val(sid), path(reads)

    output:
    tuple val(sid), path(fq_1_trimmed), path(fq_2_trimmed), emit: trimmed_reads
            file("${sid}.fastp_stats.html")

    script:
    fq_1_trimmed = sid + '_R1_P.fastq.gz'
    fq_2_trimmed = sid + '_R2_P.fastq.gz'
    """
    fastp -q 20 -l 140 --trim_poly_g --thread ${task.cpus} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]}\
    --out1 $fq_1_trimmed \
    --out2 $fq_2_trimmed \
    --json ${sid}.fastp_stats.json \
    --html ${sid}.fastp_stats.html
    """
}

/*
 * Define the `ALIGN` process that aligns reads to the reference genome
 */
process ALIGN {
    tag "$reference ${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/ALIGN"

    input:
    tuple val(sid), path(reads1), path(reads2)
    path reference
    path idx
    
    output:
    path "${sid}.bam"
    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${reads1} ${reads2} > ${sid}.bam
    """
}

/*
 * Define the `GENIDX` process that creates an index
 * given the reference genome file
 */
process GENIDX {
    tag "$reference"
    publishDir "${params.outdir}/genidx"
    
    input:
    path reference

    output:
    path "*"

    script:
    """
    samtools faidx $reference
    """
}

/*
 * Define the `CRSEQDICT` process that creates a sequence dictionary
 * for the given reference genome file
 */
process CRSEQDICT {
    tag "$reference"
    publishDir "${params.outdir}/picard"

    input:
    path reference
    path bamFile

    output:
    path "*.dict"
    path "*.bam"

    script:
    """
    java -jar /usr/local/bin/mm/share/picard-2.27.2-0/picard.jar CreateSequenceDictionary R=${reference} O=${reference.baseName}.dict
    java -jar /usr/local/bin/mm/share/picard-2.27.2-0/picard.jar AddOrReplaceReadGroups \
      I=$bamFile \
      O=${bamFile.baseName}.bam \
      --RGID 4 \
      --RGLB lib1 \
      --RGPL illumina \
      --RGPU unit1 \
      --RGSM 20
    """
}

/*
 * Define the `HAPCALL` process that performs variant calling
 * using GATK HaplotypeCaller
 */
process HAPCALL {
    tag "$reference $bamFile"
    publishDir "${params.outdir}/haplotype_caller"
	debug true
//    errorStrategy 'ignore'
	
    input:
    path reference
    path bamFile
    path idx
    path dict

    output:
    file '*.vcf'

    script:
    """
    gatk HaplotypeCaller \
        -R ${reference} \
        -I ${bamFile} \
        -O ${bamFile.baseName}.vcf
    """
}

process PREPARE {
    tag "$bamFile"
	
    input:
    path bamFile

    output:
    file '*.sorted.bam'

    script:
    """
	samtools sort $bamFile -o ${bamFile.baseName}.sorted.bam
	samtools index ${bamFile.baseName}.sorted.bam
    """
}

process CALLSNP {
    tag "$reference $bamFile"
    publishDir "${params.outdir}/bcftools"
	
    input:
    path reference
    path bamFile

    output:
    file '*.vcf'

    script:
    """
    bcftools mpileup -Ou -f $reference $bamFile | bcftools call -mv > ${bamFile.baseName}.vcf
    """
}

process ANNOTATE {
    tag "$vcf"
    publishDir "${params.outdir}/vep"
	debug true
	
    input:
    path vcf

    output:
    file '*.vep.vcf'

    script:
    """
    vep --database -i $vcf -o ${vcf.baseName}.vep.vcf 
    """
}

input_fastqs = params.reads ? Channel.fromFilePairs(params.reads, checkIfExists: true) : null

workflow {
        REFINDEX(params.reference)
        QCONTROL(input_fastqs)
        ALIGN(QCONTROL.out[0], params.reference, REFINDEX.out)
//        GENIDX(params.reference)
        PREPARE(ALIGN.out)
//        CRSEQDICT(params.reference, PREPARE.out)
//       HAPCALL(params.reference, CRSEQDICT.out[1], GENIDX.out, CRSEQDICT.out[0])
        CALLSNP(params.reference, PREPARE.out)
        ANNOTATE(CALLSNP.out)
}

workflow.onComplete {
    log.info """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
        
    log.info ( workflow.success ? "\nDone" : "\nOops" )
}




