library(patchwork)
requireNamespace("scales")
requireNamespace("stringr")

theme_image_black <- function(...) {
  theme_void(...) %+replace% 
    theme(
		panel.background = element_rect(fill ="black", color="black"),
				plot.background = element_rect(color  = 'black', fill ="black"),
				legend.key.spacing.y = unit(-2, 'mm')
    )
}

ImageDimPlotTMA <- function(obj, group.by, cols, size=1.5, ncol=NULL, nrow=NULL) {
	plots <- lapply(unique(obj$TMA_ID), function(tma) {
		sobj <- suppressWarnings(subset(obj, subset=(TMA_ID==tma)))
		title <- sobj@meta.data |> select(any_of(c("TMA", "TMA_ID", "M.nr", "Donor.Block.ID", "Donor Block ID", "Tissue_type"))) |>
			head(1) |> paste(collapse=" | ");
		ImageDimPlot(sobj, group.by=group.by, size=size) +
			scale_fill_manual(values=cols, drop=FALSE, na.value="gray30") +
			guides(fill=guide_legend(ncol=1, override.aes=list(size=2))) +
			ggtitle(title) +
			theme_image_black()
	})

	p2 <- wrap_plots(plots, ncol=ncol, nrow=nrow) & theme(plot.background = element_rect(color  = 'black', fill ="black"))

	p2
}

ImageDimPlotTMAGuide <- function(obj, group.by, flip_xy, cols=NULL, size=1.5, na.value="gray30", na.legend=TRUE, ncol=NULL, nrow=NULL, ...) {
	plots <- obj@meta.data |> group_by(TMA, TMA_ID) |> group_map(function(rows, gr) {
			
		sobj <- suppressWarnings(subset(obj, subset=(TMA_ID==gr$TMA_ID & TMA==gr$TMA)))
		
		title <- sobj@meta.data |> select(any_of(c("TMA", "TMA_ID", "M.nr", "Donor.Block.ID", "Donor Block ID", "Tissue_type"))) |>
			mutate(across(everything(), as.character)) |>
			head(1) |> paste(collapse=" | ");

		stopifnot("Unknown metadata column"=group.by %in% colnames(sobj@meta.data))

		if (is.null(cols)) {
			l <- if (is.factor(sobj@meta.data[,group.by])) {
				levels(sobj@meta.data[,group.by])
			} else {
				stringr::str_sort(unique(obj@meta.data[,group.by]), numeric=TRUE)
			}
			cols <- setNames(scales::hue_pal()(length(l)), l)
		}

		core <- ImageDimPlot(sobj, group.by=group.by, size=size, flip_xy=flip_xy, ...) +
			scale_fill_manual(values=cols, drop=FALSE, na.value=na.value) +
			guides(fill="none") +
			ggtitle(title) +
			scale_x_reverse() +
			theme(panel.background = element_rect(fill ="black", color="black"),
				plot.background = element_rect(color  = 'black', fill ="black"),
				legend.key.spacing.y = unit(-2, 'mm'),
				panel.grid=element_blank()
			)
				

		ylim <- obj@meta.data |> select(all_of(c("TMA_ID", group.by))) |>  
			dplyr::count(TMA_ID) |> pull(n) |> max()
		co <- floor(ylim / 100)

		legend <- sobj@meta.data |>  select(all_of(group.by)) |> dplyr::count(.data[[group.by]]) |>
			filter(!is.na(.data[[group.by]]) | na.legend) |>
			ggplot(aes(x=1, y=n, fill=.data[[group.by]], label=ifelse(n > co, paste(.data[[group.by]], n, sep=" : "), ""))) +
			geom_col() +
			scale_fill_manual(values=cols, guide="none", na.value=na.value) +
			geom_text(size = 3, position = position_stack(vjust = 0.5), color="white") +
			ylim(0,ylim) +
			ggtitle(group.by) +
			theme_void() + 
			theme(plot.background = element_rect(color  = 'black', fill ="black"), title=element_text(color="white"))

		core + legend + plot_layout(widths=c(.85,.15))
	})

	p2 <- wrap_plots(plots, ncol=ncol, nrow=nrow) & theme(plot.background = element_rect(color  = 'black', fill ="black"))

	p2
}


ImageFeaturePlotTMAGuide <- function(obj, features, flip_xy, cols=NULL, size=size, ncol=NULL, nrow=NULL,min.cutoff, max.cutoff) {
	plots <- obj@meta.data |> group_by(TMA, TMA_ID) |> group_map(function(rows, gr) {
			
		sobj <- suppressWarnings(subset(obj, subset=(TMA_ID==gr$TMA_ID & TMA==gr$TMA)))
		
		title <- sobj@meta.data |> select(any_of(c("TMA", "TMA_ID", "M.nr", "Donor.Block.ID", "Donor Block ID", "Tissue_type"))) |>
			mutate(across(everything(), as.character)) |>
			head(1) |> paste(collapse=" | ");

		stopifnot("Unknown metadata column"=features %in% colnames(sobj@meta.data))

		if (is.null(cols)) {
			l <- if (is.factor(sobj@meta.data[,features])) {
				levels(sobj@meta.data[,features])
			} else {
				stringr::str_sort(unique(obj@meta.data[,features]), numeric=TRUE)
			}
			cols <- setNames(scales::hue_pal()(length(l)), l)
		}

		core <- ImageFeaturePlot(sobj, features=features, size=size, min.cutoff=min.cutoff, max.cutoff=max.cutoff) +
			coord_flip() +
			#scale_fill_manual(values=cols, drop=FALSE, na.value="gray30") +
			#guides(fill="none") +
			#ggtitle(title) +
			scale_x_reverse() +
			theme(panel.background = element_rect(fill ="black", color="black"),
				plot.background = element_rect(color  = 'black', fill ="black"),
				legend.key.spacing.y = unit(-2, 'mm'))
				

	})

	p2 <- wrap_plots(plots, ncol=ncol, nrow=nrow) & theme(plot.background = element_rect(color  = 'black', fill ="black"))

	p2
}

