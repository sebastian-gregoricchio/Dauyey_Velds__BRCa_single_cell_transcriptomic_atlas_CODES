library(tidyverse)
library(readxl)
library(writexl)

scriptpath="Dauyey_Velds__BRCa_single_cell_transcriptomic_atlas_CODES"
source(file.path(scriptpath, "functions", "functions.R"))

source("scripts/load_and_subset_functions.R")
source("scripts/seurat_helpers.R")

# load sample metadata add number of epithelial cells from additional xls metadata
meta <- read_xlsx("./metadata_tma.xlsx")

rdspath <- "RDSfiles"
outpath <- "diffexp"

config <- list(
	normmethod = "LogNormalize",
	scalefactor=100,
	doscale = TRUE,
	nvfeatures = 0,
	pcadims = 30,
	nneighbors = 20,
	resolution = 1,
	nocluster = TRUE
	)

# this selects the comparison groups
gr <- read_xlsx("metadata_comparisons_all.xlsx") |> select(1,2,22:32)

meta <- meta |> left_join(gr |> select(1:2, 3:ncol(gr)))

comps <- colnames(gr)[-c(1:2)]

# Helper function to load the set of cores used in the comparison
# filter out Patients that do not contribute at leas `filtercc` cell to the comparison
LoadCompSet <- function(tab, comp, rdspath, filterbycc=0) {
	sr <- tab |> filter(!is.na(.data[[comp]]))
	print(sr)

	

	obj <- meta |> filter(`Donor Block ID` %in% unique(sr$`Donor Block ID`)) |>
		group_by(`Donor Block ID`) |>
		GroupLoadRDS(rdspath=rdspath)

	# add the spacexr major first type
	obj <- lapply(obj, function(obj) {
		obj <- try(AddSpaceXRFirstType(obj, prefix="spacexrmajor"))
		if (inherits(obj, "try-error")) {
			# this pair was not labeled due to low cell counts
			# return missing NA celltypes:
			obj@meta.data$spacexrmajor_first_type <- NA
			return(obj)
		}
		# Modify the M-12 celltype labelling to force normal to Cancer epithelial
		# This is confirmed by a pathologist
		if (unique(obj$M.nr) == "M-12")
			obj@meta.data$spacexrmajor_first_type[obj@meta.data$spacexrmajor_first_type == "Normal Epithelial"] <- "Cancer Epithelial"

		obj
	})

	obj <- lapply(obj, subset, subset=(spacexrmajor_first_type == "Cancer Epithelial"))
	
	#merge the set
	obj <- merge(obj[[1]], obj[-1])
	obj <- JoinLayers(obj)

	#add the group metadata
	cols <- c("TMA","TMA_ID", comp)
	om <- obj@meta.data[,c("TMA","TMA_ID")]
	om <- om |> left_join(sr |> select(all_of(cols)))
	obj <- AddMetaData(obj, metadata=om[,comp], col.name="Grouping" )

	# filter M nrs that have at least `filterbycc` for both labels (all cores combined)
	mkeep <- obj@meta.data |> count(M.nr, Grouping) |>
		mutate(n=ifelse(n<50, NA, n)) |>
		pivot_wider(names_from=Grouping, values_from=n) |>
		na.omit() |> pull(M.nr)
	obj <- obj[,obj$M.nr %in% mkeep]

	# process the set
	obj <- ProcessXenium(obj, config)

	Idents(obj) <- "Grouping"

	obj
}

# load all sets
all <- lapply(setNames(comps, comps), LoadCompSet, tab=meta, rdspath=rdspath)

# run findAllMarkers on all groups
markersd  <- lapply(all, FindAllMarkers, min.diff.pct=0.3, logfc.threshold=0.6, only.pos=TRUE)
markersfc  <- lapply(all, FindAllMarkers, min.diff.pct=0.0, logfc.threshold=1, only.pos=TRUE, min.pct=.3)
markersall  <- lapply(all, FindAllMarkers, min.diff.pct=0.0, logfc.threshold=0, only.pos=TRUE, min.pct=0, return.thresh=1)

write_xlsx(markersd, file.path(outpath, "markers-mindiff.xlsx"))
write_xlsx(markersfc, file.path(outpath, "markers-logfc.xlsx"))
write_xlsx(markersall, file.path(outpath, "markers-all.xlsx"))

# filter the comparison table to only keep M.nr that have at least 50 cells in both comp labels
all <- lapply(setNames(comps, comps)[1:7], LoadCompSet, tab=meta, rdspath=rdspath, filterbycc=50)

# run findAllMarkers on all groups
markersd  <- lapply(all[1:7], FindAllMarkers, min.diff.pct=0.3, logfc.threshold=0.6, only.pos=TRUE)
markersfc  <- lapply(all[1:7], FindAllMarkers, min.diff.pct=0.0, logfc.threshold=1, only.pos=TRUE, min.pct=.3)
markersall  <- lapply(all[1:7], FindAllMarkers, min.diff.pct=0.0, logfc.threshold=0, only.pos=TRUE, min.pct=0, return.thresh=1)

write_xlsx(markersd, file.path(outpath, "markers-mindiff-50cc.xlsx"))
write_xlsx(markersfc, file.path(outpath, "markers-logfc-50cc.xlsx"))
write_xlsx(markersall, file.path(outpath, "markers-all-50cc.xlsx"))
lapply(all[1:7],function(o) as.data.frame(table(o$M.nr, o$Grouping))) |>
	write_xlsx(file.path(outpath, "M-nr_50cc.xlsx"))



# Differential gene expression for CNV cluster clones M-02
obj <- LoadMandFilterCelltype(meta |> filter(M.nr == "M-02", TMA_ID==56 | TMA_ID==57), filter="Cancer Epithelial", minor=FALSE)
# add the cnv cluster metadata
obj <- AddMetadataFromFile(obj, file.path("infercnv", paste0("M-02", "-meta.Rds")))
obj$cnvclone <- NA_character_
obj$cnvclone[obj$cnv_clusters %in% c(0, 1, 5)] <- "Clone 1"
obj$cnvclone[obj$cnv_clusters %in% c(2, 3, 6)] <- "Clone 2"
Idents(obj) <- "cnvclone"
all <- FindMarkers(obj, ident.1="Clone 1", ident.2="Clone 2")
filtered <- all |> filter(p_val_adj < 0.01, abs(avg_log2FC) > .5, abs(pct.2 - pct.1) > .2 | (pct.1 > .3 & pct.2 > .3))
extra <- as.data.frame(table(obj$Donor.Block.ID, obj$cnvclone))
write_xlsx(list(all=all |> rownames_to_column("Gene"), filtered=filtered |> rownames_to_column("Gene"), extra=extra), file.path(diffoutpath, paste0("diff-cnvgroups-015_236-M-02.xlsx")))

# Differential gene expression for CNV cluster clones M-68
obj <- LoadMandFilterCelltype(meta |> filter(M.nr == "M-68"), filter="Cancer Epithelial", minor=FALSE)
# add the cnv cluster metadata
obj <- AddMetadataFromFile(obj, file.path("infercnv2", paste0("M-68", "-meta.Rds")))
Idents(obj) <- "cnv_clusters"
all <- FindAllMarkers(obj)
filtered <- all |> filter(p_val_adj < 0.01, abs(avg_log2FC) > .5, abs(pct.2 - pct.1) > .2 | (pct.1 > .3 & pct.2 > .3))
extra <- as.data.frame(table(obj$Donor.Block.ID))
write_xlsx(list(all=all, filtered=filtered, extra=extra), file.path(diffoutpath, paste0("diff-cnvclust-M-68.xlsx")))


# convert to long format for supplemental table
# combine with the older 50 cancer cell from diffexp2
# add a PASS columns to indicate the genes we use for ORA
de <- lapply(excel_sheets("markers-all-50cc.xlsx"), function(n) {
	read_xlsx("./diffexp2/markers-all-50cc.xlsx", sheet=n) |> 
		mutate(Comparison=paste(n, cluster),  .before=p_val) |>
		mutate(Pass=avg_log2FC>=1 & (pct.1+pct.2) >= .3) |>
		rename(Gene=gene) |>
		select(-cluster)
}) |> bind_rows() |> relocate(Gene, .before=p_val)


d1 <- read_xlsx("diff-cnvgroups-015_236-M-02.xlsx") |>
	mutate(Comparison="M-02 CNV cluster 0/1/5 vs. 2/3/6", .before="Gene") |>
	mutate(Pass=abs(avg_log2FC) >= 0.5 & (abs(pct.2 - pct.1) > .1 | (pct.1 > .2 & pct.2 > .2)))

d2 <- read_xlsx("diff-cnvclust-M-68.xlsx") |>
	mutate(Comparison=paste("M-68 CNV Cluster", cluster, "vs. rest")) |>
	mutate(Pass=abs(avg_log2FC) >= 0.5 & ( abs(pct.2 - pct.1) > .1 | (pct.1 > .2 & pct.2 > .2))) |>
	rename(Gene=gene) |> select(-cluster)

de |> bind_rows(d1) |> bind_rows(d2) |>
	write_xlsx("Supplemental-table-4.xlsx")



