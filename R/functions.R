#' Preprocess Seurat Object
#'
#' Performs standard pre-processing workflow for scRNA-seq data
#'
#' @param assay Assay to use
#' @param scale Perform linear transformation 'Scaling'
#' @param normalize Perform normalization
#' @param features Identify highly variable features
#' @param legacy_settings Use legacy settings
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' panc8[["RNA"]] <- seurat_preprocess(panc8[["RNA"]])
#'
seurat_preprocess <- function(assay, scale = TRUE, normalize = TRUE, features = NULL, ...) {
  # Normalize data
  
  if (normalize) {
    assay <- Seurat::NormalizeData(assay, verbose = FALSE, ...)
  }
  
  # Filter out only variable genes
  assay <- Seurat::FindVariableFeatures(assay, selection.method = "vst", verbose = FALSE, ...)
  
  # Regress out unwanted sources of variation
  if (scale) {
    assay <- Seurat::ScaleData(assay, features = rownames(assay), ...)
  }
  
  return(assay)
}

#' Run Seurat Pipeline
#'
#' This functions allows you to Preprocess, Cluster and Reduce Dimensions for a single seurat object.
#'
#' @param seu A Seurat object
#' @param assay Assay of interest in Seurat object
#' @param resolution Resolution for clustering cells. Default set to 0.6.
#' @param reduction Dimensional reduction object seu
#' @param organism Organism
#' @param ... Extra parameters passed to seurat_process
#'
#' @return
#' @export
#'
#' @examples
#'
#' processed_seu <- seurat_process(panc8)
#'
seurat_process <- function(seu, assay = "RNA", resolution = seq(0.2, 1.0, by = 0.2), reduction = "pca", organism = "human", ...) {
  assays <- names(seu@assays)
  
  assays <- assays[assays %in% c("RNA", "transcript")]
  
  for (assay in assays) {
    seu[[assay]] <- seurat_preprocess(seu[[assay]], scale = TRUE, ...)
  }
  
  # PCA
  seu <- seurat_reduce_dimensions(seu, check_duplicates = FALSE, reduction = reduction, ...)
  
  seu <- seurat_cluster(seu = seu, resolution = resolution, reduction = reduction, ...)
  
  seu <- find_all_markers(seu, seurat_assay = "RNA")
  
  # if (feature == "RNA"){
  #   enriched_seu <- tryCatch(getEnrichedPathways(seu), error = function(e) e)
  #   enrichr_available <- !any(class(enriched_seu) == "error")
  #   if(enrichr_available){
  #     seu <- enriched_seu
  #   }
  # }
  
  # annotate low read count category in seurat metadata
  seu <- add_read_count_col(seu)
  
  # annotate cell cycle scoring to seurat objects
  seu <- annotate_cell_cycle(seu, organism = organism, ...)
  
  # annotate mitochondrial percentage in seurat metadata
  seu <- add_percent_mito(seu, organism = organism)
  
  return(seu)
}

#' Dimensional Reduction
#'
#' Run PCA, TSNE and UMAP on a seurat object
#' perplexity should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation
#'
#' @param seu A Seurat object
#' @param assay Assay of interest to be run on the seurat object
#' @param reduction Set dimensional reduction object
#' @param legacy_settings Use legacy settings
#' @param ... Extra parameters passed to seurat_reduce_dimensions
#'
#' @return
#' @export
#'
#' @examples
seurat_reduce_dimensions <- function(seu, assay = "RNA", reduction = "pca", legacy_settings = FALSE, ...) {
  if ("integrated" %in% names(seu@assays)) {
    assay <- "integrated"
  } else {
    assay <- "RNA"
  }
  
  num_samples <- dim(seu)[[2]]
  
  if (num_samples < 50) {
    npcs <- num_samples - 1
  } else {
    npcs <- 50
  }
  
  if (legacy_settings) {
    message("using legacy settings")
    seu <- Seurat::RunPCA(seu, assay = assay, features = rownames(seu))
  } else {
    # seu <- Seurat::RunPCA(object = seu, do.print = FALSE, npcs = npcs, ...)
    seu <- Seurat::RunPCA(object = seu, assay = assay, features = Seurat::VariableFeatures(object = seu), do.print = FALSE, npcs = npcs, ...)
  }
  
  if (reduction == "harmony") {
    seu <- harmony::RunHarmony(seu, "batch")
  }
  
  if ((ncol(seu) - 1) > 3 * 30) {
    seu <- Seurat::RunTSNE(object = seu, assay = assay, reduction = reduction, dims = 1:30, check_duplicates = FALSE)
    seu <- Seurat::RunUMAP(object = seu, assay = assay, reduction = reduction, dims = 1:30, check_duplicates = FALSE)
  }
  
  return(seu)
}

#' Run Louvain Clustering at Multiple Resolutions
#'
#' @param seu A seurat object
#' @param resolution Clustering resolution
#' @param custom_clust
#' @param reduction Set dimensional reduction object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
seurat_cluster <- function(seu = seu, resolution = 0.6, custom_clust = NULL, reduction = "pca", algorithm = 1, ...) {
  message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] Clustering Cells..."))
  seu <- FindNeighbors(object = seu, dims = 1:30, reduction = reduction)
  
  if (length(resolution) > 1) {
    for (i in resolution) {
      message(paste0("clustering at ", i, " resolution"))
      seu <- Seurat::FindClusters(object = seu, resolution = i, algorithm = algorithm, ...)
    }
  } else if (length(resolution) == 1) {
    message(paste0("clustering at ", resolution, " resolution"))
    seu <- Seurat::FindClusters(object = seu, resolution = resolution, algorithm = algorithm, ...)
  }
  
  if (!is.null(custom_clust)) {
    seu <- Seurat::StashIdent(object = seu, save.name = "old.ident")
    clusters <- tibble::tibble("sample_id" = rownames(seu[[]])) %>%
      tibble::rownames_to_column("order") %>%
      dplyr::inner_join(custom_clust, by = "sample_id") %>%
      dplyr::pull(cluster) %>%
      identity()
    
    Idents(object = seu) <- clusters
    
    
    return(seu)
  }
  
  return(seu)
}

#' Find All Markers
#'
#' Find all markers at a range of resolutions
#'
#' @param seu A seurat object.
#' @param metavar A metadata variable to group by.
#' @param seurat_assay Assay to use, Default "RNA".
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' markers_stashed_seu <- find_all_markers(panc8)
#' marker_genes <- Misc(markers_stashed_seu, "markers")
#' str(marker_genes)
find_all_markers <- function(seu, metavar = NULL, seurat_assay = "RNA", ...) {
  if (is.null(metavar)) {
    resolutions <- colnames(seu[[]])[grepl(paste0(seurat_assay, "_snn_res."), colnames(seu[[]]))]
    
    cluster_index <- grepl(paste0(seurat_assay, "_snn_res."), colnames(seu[[]]))
    
    if (!any(cluster_index)) {
      warning("no clusters found in metadata. runnings seurat_cluster")
      seu <- seurat_cluster(seu, resolution = seq(0.2, 1.0, by = 0.2))
    }
    
    clusters <- seu[[]][, cluster_index]
    
    cluster_levels <- purrr::map_int(clusters, ~ length(unique(.x)))
    cluster_levels <- cluster_levels[cluster_levels > 1]
    
    clusters <- dplyr::select(clusters, dplyr::one_of(names(cluster_levels)))
    metavar <- names(clusters)
  }
  
  new_markers <- purrr::map(metavar, stash_marker_features, seu, seurat_assay = seurat_assay, ...)
  names(new_markers) <- metavar
  
  old_markers <- seu@misc$markers[!names(seu@misc$markers) %in% names(new_markers)]
  
  seu@misc$markers <- c(old_markers, new_markers)
  
  return(seu)
}

#' Enframe seurat markers
#'
#' @param marker_table
#'
#' @return
#' @export
#'
#' @examples
enframe_markers <- function(marker_table) {
  marker_table %>%
    dplyr::select(Gene.Name, Cluster) %>%
    dplyr::mutate(rn = row_number()) %>%
    tidyr::pivot_wider(names_from = Cluster, values_from = Gene.Name) %>%
    dplyr::select(-rn)
}

#' Stash Marker Genes in a Seurat Object
#'
#' Marker Genes will be stored in slot `@misc$markers`
#'
#' @param metavar A metadata variable to group by
#' @param seu A seurat object
#' @param seurat_assay An assay to use
#' @param top_n Use top n genes, Default "200"
#' @param p_val_cutoff p value cut-off, Default value is "0.5"
#'
#' @return
#'
#' @examples
#'
#' seu <- stash_marker_features(metavar = "batch", seu, seurat_assay = "RNA")
#'
stash_marker_features <- function(metavar, seu, seurat_assay, top_n = 200, p_val_cutoff = 0.5) {
  message(paste0("stashing presto markers for ", metavar))
  
  markers <- list()
  
  Idents(seu) <- seu@meta.data[[metavar]]
  
  DefaultAssay(seu) <- seurat_assay
  
  # markers$presto <-
  #   seu |>
  #   JoinLayers() |>
  #   Seurat::FindAllMarkers() %>%
  #   tibble::rownames_to_column("feature") |>
  #     dplyr::group_by(cluster) %>%
  #     dplyr::filter(p_val_adj < p_val_cutoff) %>%
  #     dplyr::top_n(n = top_n, wt = avg_log2FC) %>%
  #     dplyr::arrange(cluster, desc(avg_log2FC)) %>%
  #     dplyr::select(Gene.Name = feature, Average.Log.Fold.Change = avg_log2FC, Adjusted.pvalue = p_val_adj, Cluster = cluster) |>
  #   identity()
  
  markers$presto <-
    presto::wilcoxauc(seu, metavar, seurat_assay = seurat_assay) %>%
    dplyr::group_by(group) %>%
    dplyr::filter(padj < p_val_cutoff) %>%
    dplyr::top_n(n = top_n, wt = logFC) %>%
    dplyr::arrange(group, desc(logFC)) %>%
    dplyr::select(Gene.Name = feature, Average.Log.Fold.Change = logFC, Adjusted.pvalue = padj, avgExpr, Cluster = group)
  
  return(markers)
}

#' Annotate Low Read Count Category
#'
#' Add a Read Count Categorical Variable to Seurat Object (based on nCount_RNA)
#'
#' @param seu A seurat object
#' @param thresh Set a threshold for low read count
#'
#' @return
#' @export
#'
#' @examples
add_read_count_col <- function(seu, thresh = 1e5) {
  rc <- seu[["nCount_gene"]] < thresh
  
  seu <- Seurat::AddMetaData(
    object = seu,
    metadata = rc,
    col.name = "low_read_count"
  )
}

#' Annotate percent mitochondrial reads per cell
#'
#'  Add a Read Count Categorical Variable to Seurat Object (based on nCount_RNA)
#'
#' @param seu A seurat object
#'
#' @param seurat_assay gene
#'
#' @return
#' @export
#'
#' @examples
#' add_percent_mito(panc8)
#' add_percent_mito(baron2016singlecell, organism = "mouse")
add_percent_mito <- function(seu, organism = organism, seurat_assay = "gene") {
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-", assay = seurat_assay)
  
  return(seu)
}

#' Annotate Cell Cycle
#'
#' Annotate Cell Cycle for Gene and Transcript Seurat Objects
#'
#' @param seu A seurat object
#' @param organism organism. Default "human"
#'
#' @return
#' @export
#'
#' @examples
#' annotate_cell_cycle(panc8, organism = "human")
#' annotate_cell_cycle(baron2016singlecell, organism = "mouse")
annotate_cell_cycle <- function(seu, organism = "human", ...) {
  # setdefaultassay to "gene"
  # Seurat::DefaultAssay(seu) <- "gene"
  
  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
  
  if (organism == "mouse") {
    s_genes <-
      dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in% s_genes) %>%
      dplyr::pull(MGI.symbol)
    
    g2m_genes <-
      dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in% g2m_genes) %>%
      dplyr::pull(MGI.symbol)
  }
  
  seu <- CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
}