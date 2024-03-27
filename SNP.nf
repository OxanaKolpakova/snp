#!/usr/bin/env nextflow

params.reference = "/home/alexandr/Documents/SNP/data/MT.fna"
params.reads = "/home/alexandr/Documents/SNP/data/*R{1,2}.fastq"
params.outdir = "results"

log.info """\
    ===================================
    S N P   P I P E L I N E
    ===================================
    reference: ${params.reference}
    reads    : ${params.reads}
    outdir   : ${params.outdir}
    """
    .stripIndent(true)

process BWA_TRY {
    
    output:
    stdout

    script:
    """
    which bwa
    """
}
    
/*
 * Define the `BWAINDEX` process that creates an index
 * given the reference genome file
 */
process BWAINDEX {
    tag "$reference"
    publishDir "${params.outdir}/bwaindex"
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
 * Define the `FASTP` process that performs quality trimming and filtering of reads
 */
process FASTP{
    tag "${sid}"
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sid), path(reads)

    output:
    tuple val(sid), path(fq_1_trimmed), path(fq_2_trimmed), emit: trimmed_reads
            file("${sid}.fastp_stats.json")
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
 * Define the `BWAMEM` process that aligns reads to the reference genome
 */
process BWAMEM {
    tag "$reference ${sid}"
    cpus params.cpus
    memory params.memory
    publishDir "${params.outdir}/bwamem", mode: 'copy'
    
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
    publishDir "${params.outdir}/crseqdict"
    debug true

    input:
    path reference

    output:
    path "*"

    script:
    """
    java -jar /home/alexandr/picard.jar CreateSequenceDictionary R=${reference} O=${reference.baseName}.dict
    """
}

/*
 * Define the `HaplotypeCaller` process that performs variant calling
 * using GATK HaplotypeCaller
 */
process HAPCALL {
    tag "$reference $bamFile"
    publishDir "${params.outdir}/haplotype_caller"
    enabled = false

    input:
    path reference
    path bamFile

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

input_fastqs = params.reads ? Channel.fromFilePairs(params.reads, checkIfExists: true) : null

workflow {
//        BWAINDEX(params.reference)
//        FASTP(input_fastqs)
//      BWAMEM(FASTP.out[0], params.reference, BWAINDEX.out)
//        GENIDX(params.reference)
//        CRSEQDICT(params.reference)
//        HAPCALL(params.reference, BWAMEM.out)
        BWA_TRY().view()
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone" : "\nOops" )
}


