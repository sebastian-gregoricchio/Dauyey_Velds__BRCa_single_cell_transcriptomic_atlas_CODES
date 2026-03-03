library(tidyverse)
library(infercnv)
library(RColorBrewer)
library(ComplexHeatmap)
library(readxl)
library(patchwork)

scriptpath="Dauyey_Velds__BRCa_single_cell_transcriptomic_atlas_CODES"
source(file.path(scriptpath, "functions", "functions.R"))
source(file.path(scriptpath, "functions", "seurat_helpers.R"))
source(file.path(scriptpath, "functions", "plots.R"))

config <- list(
	normmethod = "LogNormalize",
	scalefactor=100,
	doscale = TRUE,
	nvfeatures = 0,
	pcadims = 30,
	nneighbors = 20,
	resolution = 1,
	nocluster = FALSE
	)

processInferCnvAssay <- function(obj) {
	DefaultAssay(obj) <- "infercnv"
	VariableFeatures(obj) <- rownames(GetAssayData(obj, assay="infercnv"))
	# copy data to scale data

	obj <- RunPCA(obj, assay="infercnv", reduction.name = "pca_cnv")

	obj <- RunUMAP(obj, dims = 1:10, reduction="pca_cnv", reduction.name="umap_cnv", npcs=30)
	obj <- FindNeighbors(obj, assay="infercnv", reduction = "pca_cnv",  dims = 1:10)
	obj <- FindClusters(obj, graph.name="infercnv_snn", cluster.name="cnv_clusters", resolution = 0.5)

}

integrateCNV <- function(obj, addcelltypes=c("Normal Epithelial", "Cancer Epithelial")) {
	mnr <- unique(obj$M.nr)
	stopifnot(length(mnr) == 1)

	# load the infercnv object
	iobj <- readRDS(file.path(infercnvpath, paste0(mnr, "/run.final.infercnv_obj")))

	# save the seurat_expression cluster under a different name
	obj$expr_clusters <- obj$seurat_clusters

	# add the infercnv data (epithelial cells only = tumor)
	ei <- Reduce(c, iobj@observation_grouped_cell_indices[addcelltypes])

	# some cells might have been from filterded cells or cores
	m <- which(colnames(iobj@expr.data) %in% colnames(obj))
	ei <- intersect(ei, m)
	if(length(ei) == 0) return(obj)

	# create seurat assay (use v4, v5 errors in processing.....)
	ia <- CreateAssayObject(data=iobj@expr.data[,ei])
	ia <- AddMetaData(ia, iobj@gene_order)

	obj[["infercnv"]] <- ia
	# copy this data matrix to scale data fort dimensional reduction
	obj@assays$infercnv@scale.data <- as.matrix(obj@assays$infercnv@data)

	# do dimensional reduction clustering and umap
	obj <- processInferCnvAssay(obj)
	obj
}

CnvComplexHeatmap <- function(obj, clustcols=tableau20) {
	chrs <- obj[["infercnv"]]@meta.features[,"chr"]
	ha <- HeatmapAnnotation(
		chr = anno_block(gp = gpar(fill ="lightgray") , labels = unique(chrs))
	)
	cols <-  circlize::colorRamp2(seq(.93, 1.07, length.out=11), rev(brewer.pal(name="RdBu", n=11)))

	# get some cell label metadata
	mat <- t(GetAssayData(obj, assay="infercnv"))
	metad <- obj@meta.data[rownames(mat),]

	o <- order(metad$cnv_clusters)
	mat <- mat[o,]
	metad <- metad[o,]
	nd <- length(unique(metad$Donor.Block.ID))
	ra <- rowAnnotation(cnv_cluster=metad$cnv_cluster, donor_block_id=metad$`Donor.Block.ID`, celltype=droplevels(metad$spacexrmajor_first_type),
		col=list(
			cnv_cluster=clustcols,
			donor_block_id=setNames(brewer.pal(n=9, "Set1")[1:nd], sort(unique(metad$Donor.Block.ID))),
			celltype=setNames(brewer.pal(n=2, "Dark2")[1:2], c("Normal Epithelial","Cancer Epithelial"))
			))

	Heatmap(mat, cluster_rows = FALSE, cluster_columns = FALSE,
		show_row_names = FALSE, show_column_names = FALSE,
		col=cols, column_title=NULL, row_title=NULL,
		column_split=chrs, row_split=metad$cnv_clusters,
		bottom_annotation=ha, left_annotation=ra) |>
	draw() |>
	grid.grabExpr()
}

PlotAll <- function(obj, cols=NULL, save=NULL) {
	d1 <- DimPlot(obj, reduction="umap", group.by="cnv_clusters", cols=tableau20)
	nd <- length(unique(obj$Donor.Block.ID))
	bcols <- setNames(brewer.pal(name="Set1", n=9)[1:nd], sort(unique(obj$Donor.Block.ID)))
	d2 <- DimPlot(obj, reduction="umap", group.by="Donor.Block.ID", cols=bcols)

	cols <- if(is.null(cols)) tableau20 else cols
	# plot the spatial layout of the CNV clusters
	sp <- ImageDimPlotTMAGuide(obj, group.by="cnv_clusters", cols=cols, nrow=2, size=1, flip_xy=FALSE)

	hm <- CnvComplexHeatmap(obj)

	pl <- (( sp / hm) | (d1 / d2)) + plot_layout(width=c(0.75,.25))

	if (!is.null(save)) ggsave(save, plot=pl, width=4000, height=2000, unit="px", dpi=100)
	invisible(pl)
}

# use dark tableau20 colors first, skip gray
tableau20 <-c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6", "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295", "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")
tableau20b <-tableau20[c(seq(1,12,2), seq(15,20,2), seq(2,13,2), seq(16,20,2))]
names(tableau20b) <- as.character(0:17)

# publication figures for the infercnv data
meta <- read_xlsx("./metadata_tma.xlsx")
rdspath <- "RDSfiles"
infercnvpath <- "infercnv"

meta |> filter(!is.na(M.nr)) |>
	group_by(M.nr) |>
	group_walk(function(rows, gr) {
		obj <- meta |> filter(M.nr == gr$M.nr) |> group_by(`Donor Block ID`) |>
			GroupLoadRDS(rdspath)

		mnr <- gr$M.nr
		# merge, but do not process
		obj <- merge(obj[[1]], obj[-1])
		obj <- JoinLayers(obj)
		print(obj)

		# add the spacexr cell types
		obj <- AddSpaceXRFirstType(obj)

		# remove any cores without eptithelial cells
		hasepi <- sapply(unique(obj$TMA_ID), function(tma)
			sum(obj$spacexrmajor_first_type[obj$TMA_ID == tma & obj$nCount_Xenium > 100] %in% c("Normal Epithelial", "Cancer Epithelial")) > 20)
		print(hasepi)

		if (!any(hasepi)) return()
		obj <- subset(obj, cells=which(obj$TMA_ID %in% names(hasepi[hasepi])))

		#process the set
		obj <- try(ProcessXenium(obj, config))
		if(inherits(obj, "try-error")) return()

		# add the inferCNV assay
		obj <- try(integrateCNV(obj))
		if(inherits(obj, "try-error")) return()

		if (!("infercnv" %in% Assays(obj))) return()

		DefaultAssay(obj) <- "Xenium"
		
		# store the CNV cluster metadata for easy access
		SaveMetadata(obj, file.path(infercnvpath, paste0(gr$M.nr, "-meta.Rds")), "cnv_clusters")


		# TMA core plots:
		pname <- function(mnr, ty, format="pdf") {
			file.path(infercnvpath, "pdfs_2", paste0(mnr, "-", ty, ".", format))
		}

		pl <- DimPlot(obj, reduction="umap", group.by="cnv_clusters", cols=tableau20b, raster=FALSE)
		ggsave(pname(mnr, "cnv_umap"), plot=pl, width=6, height=5.5)

		nd <- length(unique(obj$Donor.Block.ID))
		bcols <- setNames(brewer.pal(name="Set1", n=9)[1:nd], sort(unique(obj$Donor.Block.ID)))
		pl <- DimPlot(obj, reduction="umap", group.by="Donor.Block.ID", cols=bcols, raster=FALSE)
		ggsave(pname(mnr, "donor_umap"), plot=pl, width=6, height=5.5)

		# plot the spatial layout of the CNV clusters
		pl <- ImageDimPlotTMAGuide(obj, group.by="cnv_clusters", cols=tableau20b, nrow=2, size=1, flip_xy=FALSE, border.size=NA, na.legend=FALSE)
		nc <- length(unique(obj$TMA_ID))
		ggsave(pname(mnr, "tma_cnvclust"), plot=pl, width=6*ceiling(nc/2), height=12, dpi=125)

		pl <- CnvComplexHeatmap(obj, clustcols=tableau20b)
		ggsave(pname(mnr, "cnv_heatmap"), plot=pl, width=24, height=8)

		PlotAll(obj, cols=tableau20b, save=pname(mnr, "overview", format="png"))

}, .keep=TRUE)

