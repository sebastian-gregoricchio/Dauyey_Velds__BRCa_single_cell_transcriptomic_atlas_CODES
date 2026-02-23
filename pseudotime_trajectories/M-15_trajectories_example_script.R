## Load libraries
require(dplyr)
require(ggplot2)
require(Seurat)
require(slingshot)
require(SingleCellExperiment)
require(tradeSeq)

## Imported but not loaded libraries
# require(ggnewscale)
# require(ggtext)
# require(viridis)
# require(grid)
# require(patchwork)
# require(BiocParallel)
# require(msigdbr)
# require(clusterProfiler)
# require(purrr)


## -------------------- FUNCTIONS and Pre-defined PARAMETERS --------------------

# Params
patient_id <- "M-15"

axis.text.color <- "black"
line.color <- "black"
bkg.plot <- element_blank()

pt.size <- 3

expression_palette <- viridis::viridis(10, direction = 1)

umap_theme <- theme(aspect.ratio = 1,
                    axis.line = element_line(arrow = grid::arrow(type = "closed", ends = "last"), color = line.color),
                    axis.ticks = element_blank(),
                    axis.text = element_blank(),
                    axis.title = element_text(hjust = 0.90, color = axis.text.color),
                    panel.background = element_blank(),
                    plot.background = element_blank(),
                    legend.text = element_text(color = axis.text.color),
                    legend.background = bkg.plot,
                    legend.title = bkg.plot,
                    plot.title = ggtext::element_markdown(hjust = 0.5, color = axis.text.color),
                    plot.subtitle = ggtext::element_markdown(hjust = 0.5, color = axis.text.color))

colors_tissue_type <-
  c("Primary Breast" = "pink2",
    "Recurrent Breast" = "pink4",
    "Breast angiosarcoma" = "firebrick",
    "Malignant Pleural Effusion" = "purple",
    "Lymph Nodes" = "steelblue",
    "Ovarian" = "darkorange2",
    "Myometrium" = "orange",
    "Liver" = "chocolate4",
    "Neural" = "gold3",
    "Skin" = "forestgreen")


colors_cancer_cells_minor =
  c(# Cancer Epithelial (5 subtypes → gradient of firebrick)
    "Cancer Basal SC" = "forestgreen",
    "Cancer Cycling"  = "purple1",
    "Cancer Her2 SC"  = "steelblue4",
    "Cancer LumA SC"  = "darkorange2",
    "Cancer LumB SC"  = "gold3")


# Extract slingshot smooth curves and turn into arrow segments
add_curve_arrows = 
  function(p, sce, every = 10, length_cm = 0.15, linewidth = 0.8) {
    
    # Libraries
    require(Seurat)
    require(slingshot)
    require(SingleCellExperiment)
    require(ggplot2)
    require(grid)   
    
    
    crvs <- slingCurves(sce)
    
    list.paths <- list()
    k <- 0
    
    for (nm in names(crvs)) {
      S = as.data.frame(crvs[[nm]]$s[, c("umap_1","umap_2")])     # smoothed coords
      colnames(S) <- c("x","y")
      if (nrow(S) < 2) next
      k <- k+1
      S$Trajectory = nm
      list.paths[[k]] = S
    }
    
    paths = do.call(rbind, list.paths)
    
    # # take short segments every `every` points to avoid clutter
    # idx = seq(1, nrow(S)-1, by = every)
    # segs = data.frame(x = S$x[idx], y = S$y[idx],
    #                   xend = S$x[idx+1], yend = S$y[idx+1])
    p <-
      p +
      # geom_segment(data = segs,
      #              aes(x = x, y = y, xend = xend, yend = yend),
      #              arrow = arrow(length = unit(length_cm, "cm"), type = "closed"),
      #              lineend = "round", linewidth = linewidth) +
      ggnewscale::new_scale_color() +
      geom_path(
        data = paths,
        aes(x = x, y = y, color = Trajectory),
        linewidth = 1,
        arrow = arrow(type = "closed", length = unit(0.25, "cm"))
      )
    return(p)
  }


merged.ora.dotplot =
  function(named.ora.list,
           palette = viridis::viridis(option = "mako", n = 20, direction = +1, begin = 0.3),
           title = NULL,
           subtitle = NULL) {
    
    if (is.null(names(named.ora.list))) {
      stop("The list of ORA results provided does not have names. Use `names(named.ora.list)` = c(...)`.")
    }
    
    require(ggplot2)
    
    combined = clusterProfiler::merge_result(named.ora.list)
    
    combined_tb = as.data.frame(combined)
    
    dotplot_combo =
      enrichplot::dotplot(object = combined,
                          showCategory = 100) +
      ggtitle(label = title, subtitle = subtitle) +
      #viridis::scale_fill_viridis(option = "mako", direction = -1, begin = 0.3) +
      scale_fill_gradientn(colours = palette, name = "P<sub>adjusted</sub>") +
      xlab(NULL) +
      theme(plot.title = ggtext::element_markdown(hjust = 0.5),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5),
            legend.title = ggtext::element_markdown(),
            axis.ticks.y = element_blank()) +
      geom_text(data = combined_tb %>% mutate(Cluster_n = paste0(Cluster, "\n(", gsub(".*[/]", "", GeneRatio), ")")),
                aes(y = Description, x = Cluster_n,
                    label = Count),
                inherit.aes = FALSE)
    
    return(dotplot_combo)
  }


## ------------------------------------------------------------------------------


# Load the full table of the xenium
xenium <- readRDS(file = "./xenium_all_TMAs_combined.Rds")

xenium_pat <- subset(xenium, subset = M.nr == patient_id & Major_cell_annotation == "Cancer Epithelial" & grepl("Cancer", Minor_cell_annotation))

## Remove genes with score 0 for all cells and reloading data
cancer_pat_counts = Seurat::GetAssay(xenium_pat, assay = "Xenium")$counts
cancer_pat_counts = cancer_pat_counts[rowSums(cancer_pat_counts) > 0,]

xenium_pat =
  Seurat::CreateSeuratObject(counts = cancer_pat_counts,
                             meta.data = xenium_pat@meta.data,
                             project = "Xenium",
                             assay = gsub("-", "", paste0(patient_id, "_cancer")))

xenium_pat = Seurat::NormalizeData(object = xenium_pat)
xenium_pat = Seurat::ScaleData(object = xenium_pat)
xenium_pat = Seurat::FindVariableFeatures(object = xenium_pat)
xenium_pat = Seurat::RunPCA(object = xenium_pat)
xenium_pat = Seurat::RunUMAP(object = xenium_pat, reduction = "pca", dims = 1:50)



## ------------------------ Computed pseudotime trajectories -----------------------------
## Make Seurat clusters
clust.resolution <- 0.5
xenium_pat <- FindNeighbors(xenium_pat, dims = 1:50)
xenium_pat <- FindClusters(xenium_pat, resolution = clust.resolution, random.seed = 42, algorithm = 2)


## Compute trajectories
sce <- as.SingleCellExperiment(xenium_pat)
sce <- slingshot(data = sce,
                 dist.method = 'mnn', 
                 clusterLabels = paste0(gsub("-","", patient_id), "_cancer_snn_res.",clust.resolution),
                 reducedDim = 'UMAP',
                 # start.clus = "0", # if the start is known
                 # end.clus = "4",  # if the end is known
                 maxit = 50)

xenium_pat$pseudotime <- slingPseudotime(sce)[,1]


## -------------------------- Trajectories visualization ---------------------------------

## Make some general UMAPs
umap_tissue_origin <-
  Seurat::DimPlot(xenium_pat, group.by = "Tissue_type", raster=T, raster.dpi = c(800,800), pt.size = pt.size) + # geom_scattermore() for rasterization
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_color_manual(values = colors_tissue_type,
                     breaks = names(colors_tissue_type)) +
  ggtitle(label = paste0("**", patient_id, " Cancer cells | Tissue origin**"),
          subtitle = paste0("*n* = ", prettyNum(nrow(xenium_pat@meta.data), big.mark = ","))) +
  umap_theme


umap_cell_types_minor <-
  Seurat::DimPlot(xenium_pat, group.by = "Minor_cell_annotation", raster=T, raster.dpi = c(800,800), pt.size = pt.size) + # geom_scattermore() for rasterization
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_color_manual(values = colors_cancer_cells_minor,
                     breaks = names(colors_cancer_cells_minor)) +
  ggtitle(label = paste0("**", patient_id, " Cancer cells | cell type (minor)**"),
          subtitle = paste0("*n* = ", prettyNum(nrow(xenium_pat@meta.data), big.mark = ","))) +
  umap_theme


umap_sample_type <-
  Seurat::DimPlot(xenium_pat, group.by = "Donor.Block.ID", raster=T, raster.dpi = c(800,800), pt.size = pt.size) + # geom_scattermore() for rasterization
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle(label = paste0("**", patient_id, " Cancer cells | sample date**"),
          subtitle = paste0("*n* = ", prettyNum(nrow(xenium_pat@meta.data), big.mark = ","))) +
  umap_theme


## UMAP of Seurat clusters with trajectories
umap_seurat_clusters <-
  Seurat::DimPlot(xenium_pat, group.by = paste0(gsub("-","",patient_id),"_cancer_snn_res.",clust.resolution), raster=T, raster.dpi = c(800,800), pt.size = pt.size) +
  ggtitle( paste0("**", patient_id, " | slingshot pseudotime trajectories**")) +
  scale_color_manual(values = rainbow(length(unique(xenium_pat@meta.data[,paste0(gsub("-","",patient_id),"_cancer_snn_res.",clust.resolution)]))),
                     breaks = 0:(length(unique(xenium_pat@meta.data[,paste0(gsub("-","",patient_id),"_cancer_snn_res.",clust.resolution)]))-1)) +
  umap_theme


umap_trajectories <-
  add_curve_arrows(umap_seurat_clusters, sce, every = 8) +
  xlab("UMAP 1") +
  ylab("UMAP 2")


## Visualize combined UMAPs
patchwork::wrap_plots(umap_tissue_origin, umap_cell_types_minor,
                      umap_sample_type, umap_trajectories, nrow = 1)



## ------------------- Define the genes contribution to the lineages ---------------------

### Extract pseudotime and lineage weights for the lineage of interest
# - pt gives pseudotime for lineages (NAs for cells not assigned to it).
# - w gives soft assignment weights (0–1) for lineages. This is critical if cells share early trunk states or you have branching.

pt <- slingPseudotime(sce)     # matrix: cells x lineages
w <- slingCurveWeights(sce)    # matrix: cells x lineages

counts <- assay(sce, "logcounts")
keep <- rowSums(!is.na(pt)) > 0  # cells with at least one assigned lineage

pt2 <- pt
pt2[is.na(pt2)] <- 0



### Fit gene-wise smoothers along pseudotime for that lineage (tradeSeq)
# Defining multi-core parameters
set.seed(42)
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers = 15 # use 15 cores

# general additive model fitting (-- it takes some time --)
gam <- tradeSeq::fitGAM(counts = counts[, keep],
                        pseudotime = pt2[keep, , drop=FALSE],
                        cellWeights = w[keep, , drop=FALSE],
                        nknots = 6, # adjust (5–8 typical); more knots = more flexible curves
                        family = "gaussian", # because the xenium data are like intensities with a kind of normal distribution,
                        BPPARAM = BPPARAM)


# Test for genes associated with pseudotime in a lineage and rank them
assoc <- tradeSeq::associationTest(gam, lineages = TRUE)
assoc <- assoc[order(assoc$pvalue), , drop=FALSE]

# Test for global associations
assoc_global <- tradeSeq::associationTest(gam)
assoc_global$padj <- p.adjust(assoc_global$pvalue, "BH")
globally_sig <- unique(rownames(assoc_global)[assoc_global$padj < 0.01])
globally_sig <- globally_sig[!is.na(globally_sig)]



## Selecting the most associated genes per each lineage
list_lineage_association = list()
k=0

for (i in 1:length(grep("pvalue_",colnames(assoc)))) {
  k = k+1
  
  pval_lin = assoc[,colnames(assoc) == paste0("pvalue_",i), drop = F]
  pval_lin$gene = rownames(pval_lin)
  pval_lin$padj = p.adjust(pval_lin[,1], method = "BH")
  
  lin_genes =
    data.frame(gene = pval_lin$gene,
               pvalue = pval_lin[,1],
               padj_BH = pval_lin$padj,
               lineage = i) %>%
    dplyr::mutate(signif = padj_BH < 0.01)
  
  if (nrow(lin_genes) > 0) {list_lineage_association[[k]] = lin_genes}
}

if (length(list_lineage_association) > 1) {
  list_lineage_association_combined = do.call(rbind, list_lineage_association)
} else {
  list_lineage_association_combined = list_lineage_association
}


# Extract predicted "expression" smooth for associated genes
genes_smooth_list <- purrr::map(.x = globally_sig,
                               .f = function(x){
                                 tradeSeq::predictSmooth(models = gam, gene = x) %>%
                                   dplyr::mutate(gene = x)},
                               .progress = TRUE)

genes_smooth_combo <- do.call(rbind, genes_smooth_list)

genes_smooth_combo_assoc <-
  dplyr::left_join(x = genes_smooth_combo,
                   y = list_lineage_association_combined[,-5],
                   by = c("gene", "lineage"))




## Compute a lineage-specific effect size from the fitted curves
# lineage-specific range on the model scale

effect_tb <-
  genes_smooth_combo_assoc %>%
  dplyr::group_by(gene, lineage, padj_BH) %>%
  summarise(effect_range = max(yhat) - min(yhat), .groups="drop") %>% # scaling between min-max effect
  arrange(lineage, padj_BH, desc(effect_range))



n.top.genes <- 50

top_genes_lineage <- list()

for (i in 1:length(unique(effect_tb$lineage))) {
  top_genes_lineage[[i]] <- 
    (effect_tb %>%
       dplyr::filter(lineage == i,
                     padj_BH < 0.05))$gene[1:n.top.genes]
}

names(top_genes_lineage) <- paste0("lineage_", 1:length(unique(effect_tb$lineage)))

## Find overlap between lineages
overlap_top_genes <-
  venn::venn(x = top_genes_lineage, ggplot = TRUE, ilabels = "counts", zcolor = "style")
overlap_top_genes


# Get the unique list of genes "important"
top_genes_unique <- unique(c(unlist(top_genes_lineage, use.names = FALSE)))


## Plot the top genes expression (fitted) pattern per each trajectory
pattern_lines_top_genes <- 
  ggplot(data = genes_smooth_combo_assoc,
         aes(x = time,
             y = yhat,
             color = factor(lineage),
             group = gene)) +
  geom_line(show.legend = FALSE) +
  ylab("Fitted gene expression") +
  xlab("Pseudotime") +
  facet_wrap(~paste0("Trajectory ", lineage)) +
  ggpubr::theme_pubr()



# ---------------------------- OverRepresentation Analyses -------------------------------
## Collect signatures
hallamrks <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")

signatures_li <-
  msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP") %>%
  dplyr::filter(grepl("^LI_ESTROGENE", gs_name))


## Filter for Xenium 5001 genes panel
### The list of 5001 genes can be downloaded here --> https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://cdn.10xgenomics.com/raw/upload/v1715726653/software-support/Xenium-panels/5K_panel_files/XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv&ved=2ahUKEwi0rcr0iuiSAxX09gIHHeu-AoAQFnoECAwQAQ&usg=AOvVaw13J94Jx8EhTHZmh4Gj2eSH
fivek_panel = read.csv("./XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")[,1]
all_signatures = all_signatures %>% dplyr::filter(gene_symbol %in% fivek_panel)


## Perform ORA
ora_top_genes_lineage <-
  purrr::map(.x = 1:length(unique(effect_tb$lineage)),
             .f = function(x){
               clusterProfiler::enricher(gene = top_genes_lineage[[x]],
                                         pvalueCutoff = 0.01,
                                         pAdjustMethod = "BH",
                                         qvalueCutoff = 0.01,
                                         TERM2GENE = all_signatures %>% dplyr::select(gs_name, gene_symbol))
             })
names(ora_top_genes_lineage) <- paste0("Trajectory ", 1:length(unique(effect_tb$lineage)))


# Plot combined results
dotplot.ora_top_genes <- merged.ora.dotplot(ora_top_genes_lineage)



## Combine plots
patchwork::wrap_plots(umap_trajectories, pattern_lines_top_genes, dotplot.ora_top_genes, nrow = 1)






