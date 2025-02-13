#!/usr/bin/env nextflow

/*
========================================================================================
   CellRangerGEX module
========================================================================================
*/

nextflow.enable.dsl=2

process CellRangerGEX {
    tag {"Run CellRanger count for samples in ${params.sample_sheet}"}
    label 'process_lowmem'

    publishDir "${params.outdir}/CellRangerGEX", mode: 'copy'

    input:
    each sample
    path sequencing_directory 
    path reference_location
    val expect_cells

    output:
    path "*", emit: GEX_$sample

    script:
    """
    # Load module
    module load cellranger/7.1.0

    # Run CellRanger jobs for each sample
    cellranger count --id=$sample --transcriptome=${params.reference_location} --fastqs=${params.sequencing_directory}/$sample/ --sample=$sample --expect-cells=${params.expect_cells}
    """
}
