#!/usr/bin/env nextflow

/*
========================================================================================
   fastqc module
========================================================================================
*/

nextflow.enable.dsl=2

process fastqc {
    tag 'QC of samples in ${params.sequencing_directory}'
    label 'process_lowmem'
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    each sample
    path sequencing_directory

    output:
    path "*", emit: fastQC

    script:
    """
    module load fastqc/0.11.5
    mkdir -p $sample
    fastqc -o $sample "${params.sequencing_directory}"/$sample/*.fastq.gz
    """
}