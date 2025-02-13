#!/usr/bin/env nextflow

/*
========================================================================================
   Seurat module
========================================================================================
*/

nextflow.enable.dsl=2

process Seurat {
	publishDir "${params.outdir}/Seurat", mode: 'copy'
	label "process_highmem"
	
	input:
	path CellRanger_samples

	output:
	path "*", emit: Seurat_results

	script:
	"""
	Rscript - <<EOF
	### Load libraries
	library(Seurat)
	library(dplyr)
	library(tidyr)
	library(patchwork)
	library(harmony)
	library(ggplot2)

	# Extract folder/sample names, not the hidden folders
	setwd('${launchDir}/results/CellRangerGEX/')
	print(getwd())
	all_dirs <- list.dirs(path = getwd(), recursive = FALSE, full.names = FALSE)
	names <- Filter(function(dir) substr(dir, 1, 1) != ".", all_dirs)
	print(getwd())
	print(names)

	# Load samples as objects
	for (name in names) {
		matrixDir = getwd()
        read_10x <- paste0(name, '.data <- Read10X(data.dir = "', name, '/outs/raw_feature_bc_matrix")')
		print(read_10x)
		eval(parse(text = read_10x))
        CreatObject <- paste0(name, "<- CreateSeuratObject (counts= ", name, ".data)")
        eval(parse(text = CreatObject))
    }

	### Merge data and add cell ids
	names_vector <- names[2:length(names)]
	names_string <- paste(names_vector, collapse = ', ')
	names_vector_quote <- paste0(\"'\", names, \"'\")
	names_string_quote <- paste(names_vector_quote, collapse = ', ')
	merged <- paste0('merged <- merge(', names[1], ', y = c(', names_string, '), add.cell.ids = c(', names_string_quote, '))')
	eval(parse(text = merged))

	# Creat Seurat folder
	dir.create("../Seurat")
	# Save merged seurat object
	save(merged, file='../Seurat/merged.RData')
	
	### Filter data to keep only high-quality cells
	filtered <- subset(merged, subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 & nCount_RNA >= 200) ## Change threshold accordingly
	counts<- GetAssayData(object = filtered, slot = 'counts')
	nonzero <- counts > 0
	keep_genes <- Matrix::rowSums(nonzero) >= 10
	filtered_counts <- counts[keep_genes, ]
	filtered <- CreateSeuratObject(filtered_counts, meta.data = filtered@meta.data)
	save(filtered,file="../Seurat/filtered.RData")

	### Normalize data using LogNormalize
	filtered <- NormalizeData(filtered, normalization.method = 'LogNormalize', scale.factor = 10000)
	filtered <- FindVariableFeatures(filtered, selection.method = "vst", nfeatures = 2000)
	filtered <- ScaleData(filtered)
	filtered <- RunPCA(object= filtered)
	filtered <- FindNeighbors(filtered, dims = 1:20)
	filtered <- FindClusters(filtered, resolution = 0.4)
	filtered <- RunUMAP(filtered, dims = 1:20)
	dim_plot <- DimPlot(filtered, reduction = "umap")
	ggsave(filename = "../Seurat/filtered_norm_scaled.png", plot = dim_plot)
	save(filtered,file='../Seurat/filtered_norm_scaled.RData')

	### Normalize data using scTransform
	#filtered <- SCTransform(filtered, vars.to.regress = "percent.mt", verbose = FALSE)
	#filtered <- RunPCA(filtered, features = VariableFeatures(object = filtered))
	#filtered <- FindNeighbors(filtered, dims = 1:20)
	#filtered <- FindClusters(filtered, resolution = 0.4)
	#filtered <- RunUMAP(filtered, dims = 1:20)
	#dim_plot <- DimPlot(filtered, reduction = "umap")
	#ggsave(filename = "../Seurat/filtered_scTrans.png", plot = dim_plot)
	#save(filtered, file="../Seurat/filtered_scTrans.RData")

	### Batch correction using harmony
	harmony <- RunHarmony(filtered, group.by.vars = 'sample', dims.use = 1:20, max.iter.harmony=1) ## Adjust number of iterations accordingly
	harmony <-  FindNeighbors(harmony, dims = 1:20, verbose = F, reduction='harmony')

	### Clustering
	resolution_values <- c(0.2, 0.3, 0.4, 0.5, 0.6)  # Add more values as desired
	for (i in resolution_values) {
	harmony <- FindClusters(harmony, resolution = i)
	}
	harmony <-  RunUMAP(harmony, dims = 1:20, reduction='harmony')
	dim_plot <- DimPlot(filtered, reduction = "umap")
	ggsave(filename = "../Seurat/harmony_1.png", plot = dim_plot)
	save(harmony, file="../Seurat/harmony_1.RData")
	EOF
	"""
}
