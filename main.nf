#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Initiate parameters
params.sample_sheet = "${launchDir}/samplesheet.csv"
params.source_directory = null
params.sequencing_directory = null
params.reference_location = null
params.expect_cells = null
params.analysis_type = null
params.skip_copy_step = false  // Set this to true to skip the copy process
params.outdir = 'results'


// Import modules
include { CopySeqs } from "./modules/CopySeqs.nf"
include { fastqc } from "./modules/fastqc.nf" addParams(OUTPUT: "${params.output}/fastqc")
include { CellRangerGEX } from "./modules/CellRangerGEX.nf" addParams(OUTPUT: "${params.output}/CellRangerGEX")
include { CellRangerVDJ } from "./modules/CellRangerVDJ.nf" addParams(OUTPUT: "${params.output}/CellRangerVDJ")
include { Seurat } from "./modules/Seurat.nf" addParams(OUTPUT: "${params.output}/Seurat")


// Orchestrate the process flow
workflow {
    destination_ch = Channel.fromPath( params.sequencing_directory, checkIfExists: false)
    ref_ch = Channel.fromPath( params.reference_location, checkIfExists: true )
    type_ch = Channel.value ( params.analysis_type )

    // Parse samplesheet to extract sample names
    samples_ch = Channel.fromPath(params.sample_sheet)
            .splitCsv(header: false, sep: ',')
            .flatten()

    // CopySeqs run or not and pass the actual input for the next process
    if (!params.skip_copy_step){
        source_ch = Channel.fromPath( params.source_directory, checkIfExists: true)
        CopySeqs( source_ch, destination_ch )

        actual_input = CopySeqs.out.copy

    } else{
        actual_input = destination_ch
    }

    // Run fastqc
    fastqc( samples_ch, actual_input)

    // Rn CellRanger and Seurat depending on the analysis type
    if (params.analysis_type == "GEX") {
        num_ch = Channel.value ( params.expect_cells )
        CellRangerGEX( samples_ch, actual_input, ref_ch, num_ch)
        Seurat( CellRangerGEX.out.collect() )
    } else if (params.analysis_type == "VDJ") {
        CellRangerVDJ( samples_ch, actual_input, ref_ch)
    } else if (params.analysis_type == "GEX_and_VDJ") {
        num_ch = Channel.value ( params.expect_cells )
        CellRangerGEX( samples_ch, actual_input, ref_ch, num_ch)
        CellRangerVDJ( samples_ch, actual_input, ref_ch)
        Seurat( CellRangerGEX.out.CellRangerGEX.combine )
    } else {
        println("Error: Invalid analysis_type provided. Please specify GEX or VDJ.")
         exit 1
    }
}

// Workflow Event Handler

workflow.onComplete {

   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}
