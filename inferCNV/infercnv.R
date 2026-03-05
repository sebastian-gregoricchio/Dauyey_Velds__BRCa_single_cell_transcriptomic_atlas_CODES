library(tidyverse)
library(infercnv)

# TMA cores are grouped by M nr (patient)
# Celltype information is use to select a likely group for normal cells:

# wrapper function to run InferCNV which can be applied to all patient objects
RunInfer <- function(obj, name, outdir, metacol, normals, posfile, ...) {
	# get the count matrix
	m <- GetAssayData(obj, assay="Xenium", layer="counts")

	# get the metadata
	meta <- data.frame(
		cell=colnames(obj),
		anno=obj@meta.data[,metacol])
	dir.create(outdir, showWarnings=FALSE)
	mf <- file.path(outdir, paste0(name, "meta.txt"))
	write_delim(meta, file=mf, delim="\t", col_names=FALSE)

	infercnv_obj = CreateInfercnvObject(raw_counts_matrix=m,
		annotations_file=mf,
		delim="\t",
		gene_order_file=posfile,
		ref_group_names=normals)

	infercnv_obj = infercnv::run(infercnv_obj,
		cutoff=0.02,  # use 1 for smart-seq, 0.1 for 10x-genomics. changed to 0.02 for Xenium
		out_dir=outdir,  # dir is auto-created for storing outputs
		cluster_by_groups=TRUE,   # cluster
		denoise=TRUE,
		HMM=FALSE,
		...
	)
	# not need to save. infercnv saves as run.final.infercnv_obj in output dir
	#saveRDS(infercnv_obj, file.path(outdir, paste0(name, "-infercnv.Rds")))
}

# required functions
scriptpath="Dauyey_Velds__BRCa_single_cell_transcriptomic_atlas_CODES"
source(file.path(scriptpath, "functions", "functions.R"))
source(file.path(scriptpath, "functions", "seurat_helpers.R"))

# The analysis folder
ifpath <- "infercnv"

# load sample metadata
meta <- read_xlsx("./metadata_tma.xlsx")

rdspath <- "RDSfiles"

# load the gene position information
posfile <- file.path(ifpath, "https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt")

# the celltypes to use as normal reference
normals <- c("B-cells" ,"CAFs","Endothelial","Myeloid","Plasmablasts", "PVL", "T-cells")

# analysis mode switched to samples to avoid errors
# removing small groups might also helpt, but we don't use the subclustering
# and the hmm prediction so it might be good enough
meta |> filter(!is.na(M.nr)) |> group_by(M.nr) |> group_walk(function(rows, gr) {
	# load and combine the cores:
	objs <- rows |> group_by(`Donor Block ID`) |>
		GroupLoadRDS(rdspath)

	# merge, but do not process
	all <- merge(objs[[1]], objs[-1])
	all <- JoinLayers(all)

	# add the major cell type labels from SpaceXR
	all <- AddSpaceXRFirstType(all)

	# change type factor to chr to avoid error in infercnv when 0 cells pass
	all$spacexrmajor_first_type <- as.character(all$spacexrmajor_first_type)
	all <- subset(all, cells=which(!is.na(all$spacexrmajor_first_type) & all$nCount_Xenium > 110))


	out <- file.path(ifpath, gr$M.nr)
	RunInfer(all, gr$M.nr, out, "spacexrmajor_first_type", intersect(normals, unique(all$spacexrmajor_first_type)), posfile, analysis_mode="samples")
}, .keep=TRUE)

