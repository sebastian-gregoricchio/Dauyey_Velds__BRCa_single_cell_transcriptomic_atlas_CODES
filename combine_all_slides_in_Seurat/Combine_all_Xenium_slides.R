library(Matrix)
library(glmGamPoi)
library(SeuratObject)
library(Seurat)
library(SeuratWrappers) # (+ banksy installed)
library(SingleCellExperiment)
library(SpatialExperiment)
library(SummarizedExperiment)

## dev MAtrix version is required
# install.packages("Matrix", repos = "http://R-Forge.R-project.org") v >1.8.0


#### list of RDS files
rds_files_tma <- c(list.files("./core-RDSfiles/TMA13", full.names = TRUE),
                   list.files("./core-RDSfiles/TMA15", full.names = TRUE),
                   list.files("./core-RDSfiles/TMA17", full.names = TRUE))


# function definition to collect individual counts
read_counts <- function(x) {
  id <- stringr::str_extract(basename(x), "TMA\\d+_\\d+")

  object <- readRDS(x)

  mat <- SeuratObject::LayerData(object, assay = "Xenium", layer = "counts")
  colnames(mat) <- paste(colnames(mat), id, sep = "_")

  rm(object); gc() # to save memory space
  return(mat)
}


# function definition to collect individual metadata
read_metadata <- function(path) {
  id <- stringr::str_extract(basename(path), "TMA\\d+_\\d+")

  object <- readRDS(path)

  meta <- object@meta.data
  rm(object); gc()

  meta$tma_id  <- id
  meta$cell_id <- paste(id, rownames(meta), sep = "_")
  return(meta)
}


# collect all metadata
metadata_list <- purrr::map(rds_files_tma, read_metadata, .progress = TRUE)

# Combine metadata in one table
metadata <- dplyr::bind_rows(metadata_list)
stopifnot(!anyDuplicated(metadata$cell_id))
rownames(metadata) <- metadata$cell_id


# collect all raw counts in a list (dgCMatrix)
counts_list <- purrr::map(rds_files_tma, read_counts, .progress = TRUE)

# combine the dgCMatrices in a unique one
merged_counts <- SeuratObject::RowMergeSparseMatrices(mat1 = counts_list[[1]],
                                                      mat2 = counts_list[-1])


# Some checks
nrow(merged_counts)   # ---> should equal to the panel size (n = 5001)
anyDuplicated(colnames(merged_counts))  # -----> must be 0


#### load in Seurat
xenium <-
  Seurat::CreateSeuratObject(
    counts = merged_counts,
    meta.data = metadata,
    project = "Xenium",
    min.cells = 0,
    min.features = 0,
    assay = "Xenium")

xenium <- Seurat::NormalizeData(xenium,
                                normalization.method = "LogNormalize",
                                scale.factor = 100)


# Dimensional reduction (PCA and UMAP)
VariableFeatures(xenium) <- rownames(xenium)   # v1 panel; use FindVariableFeatures() on Prime 5K
xenium <- ScaleData(xenium)
xenium <- RunPCA(xenium, npcs = 50)

options(future.globals.maxSize = 8 * 1024^3)   # 8 GB
xenium <- FindNeighbors(xenium, dims = 1:50, k.param = 16)

xenium <- FindClusters(xenium, algorithm = 1)  # Louvain
xenium <- Seurat::RunUMAP(xenium, reduction = "pca", dims = 1:50, verbose = TRUE)



## Make cenrtroids per each TMA to avoid the overlaps
for (t in unique(xenium$TMA)) {
  cells <- colnames(xenium)[xenium$TMA == t]

  cents <- CreateCentroids(data.frame(x = xenium$x_centroid[cells],
                                      y = xenium$y_centroid[cells],
                                      cell = cells))

  xenium[[t]] <- CreateFOV(coords = cents,
                           type = "centroids",
                           assay = "Xenium")
}


# Centroids for display only
offsets <- c(TMA13 = 0, TMA15 = 15000, TMA17 = 30000)

xenium$x_centroid_display <- xenium$x_centroid + offsets[as.character(xenium$TMA)]


xenium[["all_TMA_display"]] <-
  CreateFOV(coords = CreateCentroids(data.frame(x = xenium$x_centroid_display,
                                                y = xenium$y_centroid)),
            type = "centroids",
            assay = "Xenium")


saveRDS(object = xenium,
        file =  "./xenium_all_TMAs_combined.Rds")
