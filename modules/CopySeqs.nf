#!/usr/bin/env nextflow

/*
========================================================================================
   CopySeqs module
========================================================================================
*/

nextflow.enable.dsl=2

process CopySeqs {
    tag {"copy ${params.source_directory}"}
    label 'process_lowmem'

    input:
    path source_directory  
    path sequencing_directory

    output:
    path sequencing_directory, emit: copy

    script:
    """
    # Create the sequencing_directory if it doesn't exist
    mkdir -p ${params.sequencing_directory}

    # Copy the source_directory to the destination directory using rsync
    rsync -avz --progress ${params.source_directory}/* ${params.sequencing_directory}
    
    """
}
