#!/bin/bash

source /etc/profile

#$ -o scRNA-seq.out
#$ -e viralrecon_training.err
#$ -N viralrecon_training
#$ -cwd

module load miniconda3/20230728
conda deactivate
conda activate scRNA-seq

module load nextflow/22.10.6

nextflow run main.nf \
	--source_directory /path/to/original/sample/reads/
	--sequencing_directory /path/to/desired/sample/reads 
	--reference_location /path/to/reference
	--skip_copy_step <false or true> 
	--expect_cells <a number>
	--analysis_type <GEX or VDJ or GEX_and_VDJ> 
	-c nextflow.config