#!/usr/bin/env nextflow

/*
========================================================================================
   CellRangerVDJ module
========================================================================================
*/

nextflow.enable.dsl=2

process CellRangerVDJ {
    tag {"Run CellRanger VDJ for samples in ${params.sample_sheet}"}
    label 'process_lowmem'

    publishDir "${params.outdir}/CellRangerVDJ", mode: 'copy'

    input:
    each sample
    path sequencing_directory 
    path reference_location

    output:
    path "*", emit: CellRangerVDJ

    script:
    """
    # Load module
    module load cellranger/7.1.0

    # Run CellRanger jobs for each sample
    cellranger vdj --id=$sample --reference=${params.reference_location} --fastqs=${params.sequencing_directory}/$sample/ --sample=$sample
    """
}
