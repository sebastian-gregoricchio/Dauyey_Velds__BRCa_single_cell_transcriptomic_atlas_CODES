library(Seurat)

# We use the raw counts from the Swarbrick dataset as source
# See swarbrick_ref.Rds
# Download and extract GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz

swpath <- "./"

# The files need to be renamed to be compatible with Seurat
file.rename("count_matrix_barcodes.tsv", "barcodes.tsv")
file.rename("count_matrix_genes", "genes.tsv")
file.rename("count_matrix_sparse.mtx", "matrix.mtx")

# read the metadata
m <- read.csv("metadata.csv", row.names=1)
# Create the Seurat Object
brc <- CreateSeuratObject(Read10X(data.dir="./", gene.column=1))
brc<- AddMetaData(brc, m)

saveRDS(brc, "swarbrick_raw.Rds")




