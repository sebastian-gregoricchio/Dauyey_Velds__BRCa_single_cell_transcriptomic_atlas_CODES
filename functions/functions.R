# Commonly uses helper functions in the analysis steps

# example: TMA13_A24_60033_0043800_10_11_M-29_E16-60069_A1.Rds
# Rds files are create from the provided rows of the metadata
SampleRdsName <- function(rows, path=NULL) {
	# replace the TMA in the slide with the TMA ID
	stopifnot("Cannot make file name for multiple TMA's"=length(unique(rows$TMA)) == 1)
	stopifnot("Cannot make file name for M.nr"=length(unique(rows$M.nr)) == 1)

	tma <- sub("^TMA", unique(rows$TMA), unique(rows$Slidename))
	xenid <- unique(rows$Xenium_slide)
	mnr <- unique(rows$M.nr)
	bid <- unique(rows$`Donor Block ID`)

	fname <- paste(tma, xenid, paste(rows$TMA_ID, collapse="_"), mnr, sep="_") |>
		paste(ifelse(length(bid) == 1, gsub(" ","_", bid), "multi"), sep="_") |>
		paste0(".Rds")

	if (!is.null(path)) {
		file.path(path, fname)
	} else {
		fname
	}
}

# Load the spacexr metadata output by M.nr and update
# the provided Seurat metadata
AddSpaceXRFirstType <- function(obj, prefix="spacexrmajor", minor=FALSE) {
	M.nr <- unique(obj$M.nr)
	print(M.nr)
	stopifnot(length(M.nr) == 1)

	spacexrpath = ifelse(minor,"./swarbrick/RCTD_meta_by_dbi_minor/", "./swarbrick/RCTD_meta_by_dbi/")
	space_files <- list.files(path=spacexrpath, pattern=paste0("^", M.nr,".*_meta.Rds"), full.names=TRUE)
	print(space_files)
	if (length(space_files) == 0) {
		return(AddMetaData(obj, rep(NA, ncol(obj)), col.name=paste(prefix, "first_type", sep="_")))
	}

	obj_space <- lapply(space_files, readRDS) |> bind_rows()
	if (sum(rownames(obj_space) %in% colnames(obj)) == 0)
		return(AddMetaData(obj, rep(NA, ncol(obj)), col.name=paste(prefix, "first_type", sep="_")))

	print(head(obj_space))
	# convert to factor 
	AddMetaData(obj, metadata=obj_space |> rename_with(function(n) paste(prefix, n, sep="_")))
}

# Save (selected) metadata columns to file
SaveMetadata <- function(obj, file, cols=colnames(obj@meta.data)) {
	nc <- setdiff(cols, colnames(obj@meta.data))
	if (length(nc) > 0) {
		stop("Column(s): '", nc, "' not in metadata")
	}

	m <- obj@meta.data[,cols, drop=FALSE]
	saveRDS(m, file)
}

# Add 
AddMetadataFromFile <- function(obj, file, overwrite=FALSE, cols=NULL) {
	tab <- readRDS(file)

	cols <- if (!is.null(cols) & length(cols) > 0) {
		nc <- setdiff(cols, colnames(tab))
		if (length(nc) > 0) {
			stop("Column(s): (", paste(nc, sep=", "), ") not in metadata")
		}
	} else {
			colnames(tab)
	}

	ic <- intersect(colnames(tab), colnames(obj@meta.data))
	if (length(ic) > 0 & !overwrite) {
		stop("Overlapping column names: (", paste(ic, sep=", "), "). Use overwrite=TRUE to force")
	}

	if (nrow(tab) != ncol(obj))
		warning("Different number of rows in Seurat object and metadata file")

	AddMetaData(obj, metadata=tab[,cols, drop=FALSE])
}

