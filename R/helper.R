#' Merge Two Seurat Objects
#'
#' This function merges two Seurat objects, adds unique cell identifiers for
#' each object, and saves the resulting merged Seurat object as an RDS file.
#'
#' @param sobj_1 A Seurat object to be merged.
#' @param sobj_2 A Seurat object to be merged.
#' @param ids A character vector of length 2 containing unique cell identifiers
#'   for each Seurat object.
#' @param save_dir A string specifying the directory where the merged Seurat
#'   object should be saved.
#'
#' @return A Seurat object containing the merged data.
#'
#' @export
merge_seurat_objects <- function(sobj_1,
                                 sobj_2,
                                 ids,
                                 save_dir) {
  # Merges the two Seurat objects and adds the provided cell identifiers
  sobj_liver_merged <- merge(
    sobj_1,
    y = sobj_2,
    add.cell.ids = ids
  )

  # Saves the merged Seurat object
  saveRDS(sobj_liver_merged, file = file.path(save_dir, "sobj_merged_liver.Rds"))

  sobj_liver_merged
}

#' Get the Border Boundaries of Cells in Given FOVs
#'
#' This function calculates the boundaries (4 coordinates) of the cells within
#' the specified Fields of View (FOVs) from a Seurat object, based on the tissue
#' and FOV information.
#'
#' @param sobj A Seurat object containing spatially resolved transcriptomics data.
#' @param fovs A vector of FOV identifiers to specify which regions of the tissue to consider.
#'
#' @return A 2x2 matrix containing the minimum and maximum coordinates for the
#'   x and y axes (bounding box) of the cells in the specified FOVs.
#'
#' @export
get_border_boundaries <- function(sobj, fovs) {
  box::use(
    Seurat[Images]
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

#' Plot Cells in a Given FOV
#'
#' This function generates a 2D plot of cells in a specific Field of View (FOV)
#' from a Seurat object, displaying the boundary box around the specified cells
#' and focusing on the area defined by the given FOVs.
#'
#' @param sobj A Seurat object containing spatial data.
#' @param fovs A vector of FOV identifiers specifying which regions to visualize.
#'
#' @return A ggplot object that visualizes the cells within the boundaries of
#'   the specified FOVs, with the boundary box drawn around the cells of
#'   interest.
#'
#' @export
plot_cells_in_fov <- function(sobj, fovs) {
  box::use(
    ggplot2[xlim, ylim],
    Seurat[ImageDimPlot, Images]
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

#' Plot Gene Expression Markers in a Given FOV
#'
#' This function generates a 2D plot of gene expression markers within a given
#' Field of View (FOV) in a Seurat object. The plot shows cells colored based on
#' the expression of specified genes, with boundary boxes highlighting the cells
#' in the selected FOV.
#'
#' @param sobj A Seurat object containing spatial transcriptomics data.
#' @param fovs A vector of FOV identifiers specifying which regions to visualize.
#' @param genes A character vector containing the names of the genes to be visualized in the plot.
#'
#' @return A ggplot object displaying the gene expression markers in the
#'   specified FOVs, with boundary boxes drawn around the cells of interest and
#'   the gene expression levels visualized.
#'
#' @export
plot_gene_expression_markers <- function(sobj, fovs, genes) {
  box::use(
    ggplot2[xlim, ylim],
    Seurat[ImageDimPlot, Images]
  )

  fov_boundaries <- get_border_boundaries(sobj, fovs)
  slide_name <- Images(sobj)[1]

  ImageDimPlot(
    sobj,
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
