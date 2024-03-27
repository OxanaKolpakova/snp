#!/usr/bin/env nextflow

params.reference = "/home/alexandr/Documents/SNP/data/MT.fna"
params.reads = "/home/alexandr/Documents/SNP/data/*R{1,2}.fastq"
params.outdir = "results"
params.short = false // If the 'short' parameter is passed, perform alignment without fastp

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
 * Define the `BWAMEM` process that aligns reads to the reference genome
 */

process FASTP{
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

process BWAMEM {
    tag "BWAMEM"
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


if (params.reads != false) {
            Channel.fromFilePairs(params.reads, checkIfExists: true )
                .set { input_fastqs }
        }

workflow {
        BWAINDEX(params.reference)
        FASTP(input_fastqs)
        BWAMEM(FASTP.out[0], params.reference, BWAINDEX.out)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone" : "\nOops" )
}





