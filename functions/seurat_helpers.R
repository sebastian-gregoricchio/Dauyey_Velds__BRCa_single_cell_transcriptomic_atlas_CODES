library(Seurat)

# nornalize and cluster a Seurat object
ProcessXenium <- function(obj, config, verbose=TRUE) {
	obj <- subset(obj, subset=(nCount_Xenium >= 40 & nFeature_Xenium >= 15))
	stopifnot(ncol(obj) > 30)

	obj <- obj[rowSums(GetAssayData(obj, layer="counts")) > 0,]

	obj <- NormalizeData(obj, scale.factor=config$scalefactor, normalization.method=config$normmethod, verbose=verbose)
	obj <- ScaleData(obj, do.center=TRUE, do.scale=config$doscale, verbose=verbose)


	VariableFeatures(obj) <- if (config$nvfeatures > 0)
		VariableFeatures(FindVariableFeatures(dt, nfeatures=config$nvfeatures))
	else
		rownames(obj)

	if(config$nocluster) {
		obj
	} else {
		obj <- RunPCA(obj, npcs = 30, verbose=verbose)
		obj <- RunUMAP(obj, dims = 1:config$pcadims, verbose=verbose)
		obj <- FindNeighbors(obj, reduction = "pca", dims = 1:config$pcadims, verbose=verbose)
		obj <- FindClusters(obj, resolution = config$resolution, verbose=verbose)
		obj
	}
	obj
}

# process a list of Seurat objects into a single normalized/clustered Seurat object
ProcessXeniumSet <- function(obj, config, verbose=TRUE) {
	all <- merge(obj[[1]], obj[-1])
	all <- JoinLayers(all)
	ProcessXenium(all, config, verbose)
}

