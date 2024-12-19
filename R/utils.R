#' @export
get_border_boundaries <- function(sobj, fovs) {
  box::use(
    Seurat[...]
  )

  # Images() returns the slide name --> Only one value
  slide_name <- Images(sobj)[1]
  # Returns the first tissue name in a vector
  # Apparently all the same value
  tissue_name <- sobj@meta.data$Run_Tissue_name[1]

  # `id` returns the cell id
  # `fov` must be part of the given FOVs (only take chosen ones)
  # Tissue must be identical?
  # Returns all cells that are part of the FOVs
  cells_in_boundary <- sobj$id[
    (sobj$fov %in% fovs) &
      (sobj$Run_Tissue_name == tissue_name)
  ]

  # Returns all centroids for the slide
  centroids <- sobj@images[[slide_name]]$centroids
  # Renames all cells according to the slide
  centroids@cells <- paste(slide_name, centroids@cells, sep = "_")

  # Returns 4 coordinates to build a boundary box for the cells of interest
  border_boundaries <- apply(
    centroids@coords[centroids@cells %in% cells_in_boundary, ],
    2,
    range
  )

  border_boundaries
}

#' @export
plot_cells <- function(sobj, fovs) {
  box::use(
    ggplot2[xlim, ylim],
    Seurat[...]
  )

  # Get the border boundaries for ours FOVs
  fov_boundaries <- get_border_boundaries(sobj, fovs)
  slide_name <- Images(sobj)[1]

  ImageDimPlot(
    sobj,
    fov = slide_name,
    border.color = "black"
  ) +
    xlim(fov_boundaries[, 2]) +
    ylim(fov_boundaries[, 1])
}

#' @export
plot_gene_expression_markers <- function(sobj, fovs, genes) {
  box::use(
    ggplot2[xlim, ylim],
    Seurat[ImageDimPlot]
  )

  fov_boundaries <- get_border_boundaries(sobj, fovs)
  slide_name <- Images(sobj)[1]

  ImageDimPlot(
    sojb,
    fov = slide_name,
    border.color = "black",
    alpha = 0.3, # Reduce alpha of cell fills to improve molecule visualization
    molecules = genes,
    mols.size = 0.8,
    nmols = 100000, # Set the total number of molecules to visualize
    axes = FALSE
  ) +
    xlim(fov_boundaries[, 2]) +
    ylim(fov_boundaries[, 1])
}
