library(spacexr)
library(tidyverse)
library(readxl)

## Restrict the number of cores the blas library can use or it will saturate the server!
library(RhpcBLASctl)
Sys.setenv(OMP_NUM_THREADS = "8", MKL_NUM_THREADS = "8",
           OPENBLAS_NUM_THREADS = "8", NUMEXPR_NUM_THREADS = "8", R_MAX_NUM_THREADS = "8")
blas_set_num_threads(8)
omp_set_num_threads(8)


# load helper functions
scriptpath="Dauyey_Velds__BRCa_single_cell_transcriptomic_atlas_CODES"
source(file.path(scriptpath, "functions", "functions.R"))


rdspath <- "./RDSfiles/"
outpath <- "./swarbrick/"

# read the core metadata
meta <- read_xlsx("./metadata_tma.xlsx")

# edit config
config <- list(
	normmethod = "LogNormalize",
	scalefactor=100,
	doscale = TRUE,
	nvfeatures = 0,
	pcadims = 30,
	nneighbors = 20,
	resolution = 1
	)

# We use the raw counts from the Swarbrick dataset as source
# See swarbrick_ref.R

swarbrick1 <- readRDS(file="swarbrick_raw.Rds")
class_df <- swarbrick1@meta.data |> select(celltype_major, celltype_minor) |>
	unique() |>
	remove_rownames() |>
	column_to_rownames("celltype_minor") |>
	rename(class=celltype_major)

####run swarbrick by donor Block ID and use major celltypes
ref_labels <- swarbrick1$celltype_major %>% as.factor()

dir.create(file.path(outpath, "RCTD_meta_by_dbi_major")
dir.create(file.path(outpath, "RCTD_meta_by_dbi_minor")

meta |> filter(!is.na(M.nr)) |> 
	group_by(`Donor Block ID`, Tissue_type) |> 
	group_walk(function(rows, gr) {
	print(gr)
	print(rows)
	M.nr <- unique(rows$M.nr)

	# skip if done
	if (file.exists(file.path(outpath, "RCTD_meta_by_dbi",
			paste0("precomp_swarbrick_RCTD_", M.nr,"_", gr$`Donor Block ID`,".rds")))) return()

	# load and combine the cores:
	obj <- rows |> GroupLoadRDS(rdspath)
	
	if (ncol(obj[[1]]) < 500) return()
	obj <- try(ProcessXenium(obj[[1]], config))
	# can fail on low cell numbers
	if(inherits(obj, "try-error")) return()

	common_genes <- intersect(rownames(swarbrick1), rownames(obj))

	ref.obj <- Reference(GetAssayData(swarbrick1, assay="RNA", layer="counts")[common_genes, ],
                     cell_types = ref_labels, min_UMI = 10, require_int = FALSE)

	test.obj <- SpatialRNA(coords = obj@meta.data[,c("x_centroid","y_centroid")] %>% as.data.frame(),
                       counts = GetAssayData(obj, assay = "Xenium", layer = "counts")[common_genes, ],
                       require_int = FALSE)

 	RCTD <- create.RCTD(
 		test.obj, 
  		ref.obj, 
  		UMI_min = 20, 
  		counts_MIN = 5, 
  		UMI_min_sigma = 100,
  		max_cores = 8, 
  	)
  
  	RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  	
	saveRDS(file=file.path(outpath, "RCTD_meta_by_dbi_major",
			paste0("precomp_swarbrick_RCTD_", M.nr,"_", gr$`Donor Block ID`,".rds")),
		RCTD)
	saveRDS(file=file.path(outpath, "RCTD_meta_by_dbi_major",
			paste0(M.nr, "_", gr$`Donor Block ID`, "_swarbrick_meta.Rds")),
		RCTD@results$results_df)
}, .keep=TRUE)


# and for the minor cell types
ref_labels <- swarbrick1$celltype_minor %>% as.factor()
meta |> filter(!is.na(M.nr), TMA=="TMA17") |>
	group_by(`Donor Block ID`, Tissue_type) |>
	group_walk(function(rows, gr) {
	print(gr)
	print(rows)
	M.nr <- unique(rows$M.nr)

	# skip if done
	if (file.exists(file.path(outpath, "RCTD_meta_by_dbi_minor",
			paste0("precomp_swarbrick_minor_RCTD_", M.nr,"_", gsub(" ", "_",gr$`Donor Block ID`),".rds")))) return()

	# load and combine the cores:
	obj <- rows |> GroupLoadRDS(rdspath)
	
	if (ncol(obj[[1]]) < 500) return()
	obj <- try(ProcessXenium(obj[[1]], config))
	# can fail on low cell numbers
	if(inherits(obj, "try-error")) return()
	

	common_genes <- intersect(rownames(swarbrick1), rownames(obj))

	ref.obj <- Reference(GetAssayData(swarbrick1, assay="RNA", layer="counts")[common_genes, ],
                     cell_types = ref_labels, min_UMI = 10, require_int = FALSE)

	test.obj <- SpatialRNA(coords = obj@meta.data[,c("x_centroid","y_centroid")] %>% as.data.frame(),
                       counts = GetAssayData(obj, assay = "Xenium", layer = "counts")[common_genes, ],
                       require_int = FALSE)

 	RCTD <- create.RCTD(
 		test.obj, 
  		ref.obj, 
  		UMI_min = 20, 
  		counts_MIN = 5, 
  		UMI_min_sigma = 100,
  		max_cores = 12, 
  		class_df = class_df # highly recommended 
  	)
  
  	RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  	
	saveRDS(file=file.path(outpath, "RCTD_meta_by_dbi_minor",
			paste0("precomp_swarbrick_minor_RCTD_", M.nr,"_", gsub(" ", "_",gr$`Donor Block ID`),".rds")),
		RCTD)
	saveRDS(file=file.path(outpath, "RCTD_meta_by_dbi_minor/",
			paste0(M.nr, "_", gsub(" ", "_",gr$`Donor Block ID`), "_swarbrick_minor_meta.Rds")),
		RCTD@results$results_df)

}, .keep=TRUE)

