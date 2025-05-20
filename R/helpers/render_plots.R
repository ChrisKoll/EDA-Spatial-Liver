#' Custom categorical color palette
#'
#' @description
#' A predefined color palette optimized for visualizing categorical data in
#' plots. The palette is based on personal branding colors, starting with a
#' dark blue `#1D2951`.
#'
#' @details
#' This palette contains 12 distinct colors and is intended for use in functions
#' that visualize multiple groups.
#'
#' @keywords internal
custom_categorical <- c(
  "#1D2951", # Space Cadet (dark blue)
  "#FFD166", # Warm Yellow
  "#4BA3C3", # Sky Blue
  "#F25C54", # Coral Red
  "#007C91", # Deep Teal
  "#EF798A", # Dusty Rose
  "#6DCDB8", # Mint Green
  "#FF9F1C", # Vivid Orange
  "#6B4C9A", # Violet Indigo
  "#A4C400", # Lime Green
  "#8E7DBE", # Lavender Purple
  "#2E8B57" # Sea Green
)

#' Get spatial boundaries for FOVs
#'
#' @description
#' `get_spatial_boundaries()` returns the x and y spatial limits (bounding box)
#' of cells within specified fields of view (FOVs) in a `Seurat` spatial object.
#'
#' @details
#' The function extracts centroid coordinates for all cells within the given
#' FOVs and computes the spatial extent (minimum and maximum) in both x and y
#' dimensions.
#'
#' @param sobj A `Seurat` object containing spatial transcriptomics data.
#' @param fovs A character or numeric vector of FOV identifiers to include in
#'   the boundary calculation.
#'
#' @return A `2 x 2` matrix of numeric values representing the spatial range:
#'
#'   * Row 1: `y`-axis min and max.
#'   * Row 2: `x`-axis min and max.
#'
#' @keywords internal
get_spatial_boundaries <- function(sobj, fovs) {
  box::use(
    Seurat[Images]
  )


  # Images() returns the slide name --> Only one value
  slide_name <- Images(sobj)[1]
  # Returns the first tissue name in a vector
  tissue_name <- sobj@meta.data$Run_Tissue_name[1]

  # `id` returns the cell id
  # `fov` must be part of the given FOVs (only take chosen ones)
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
  spatial_boundaries <- apply(
    centroids@coords[centroids@cells %in% cells_in_boundary, ],
    2,
    range
  )

  spatial_boundaries
}

#' Render cell map for FOVs
#'
#' @description
#' `render_cell_map()` plots the spatial distribution of all cells in the
#' specified FOVs from a `Seurat` object. It includes a bounding box to zoom the
#' plot to the relevant region of interest.
#'
#' @details
#' The function uses `Seurat::ImageDimPlot()` to render spatial positions of
#' cells, focusing the viewport on the selected FOVs. The bounding box is
#' automatically calculated based on cell centroids.
#'
#' @param sobj A `Seurat` object containing spatial data.
#' @param fovs A character or numeric vector of FOV identifiers specifying
#'   which regions to visualize.
#'
#' @return A `ggplot` object displaying the spatial cell layout with a bounding
#'   box over the selected FOVs.
#'
#' @export
render_cell_map <- function(sobj, fovs) {
  box::use(
    ggplot2[labs, xlim, ylim],
    Seurat[ImageDimPlot, Images]
  )

  # Get the border boundaries for ours FOVs
  fov_boundaries <- get_spatial_boundaries(sobj, fovs)
  slide_name <- Images(sobj)[1]

  ImageDimPlot(
    sobj,
    fov = slide_name,
    border.color = "#000000"
  ) +
    xlim(fov_boundaries[, 2]) +
    ylim(fov_boundaries[, 1]) +
    labs(fill = "Cell Types")
}

#' Render transcript map for FOVs
#'
#' @description
#' `render_transcript_map()` visualizes transcripts in the specified FOVs from
#' a `Seurat` spatial object. The function overlays molecule-level gene
#' expression data onto a spatial cell layout, with a bounding box highlighting
#' the region of interest.
#'
#' @details
#' This plot is generated using `Seurat::ImageDimPlot()` with additional
#' arguments to control molecule visualization. It is particularly useful for
#' spatial marker analysis across specific tissue regions.
#'
#' @param sobj A `Seurat` object containing spatial transcriptomics data.
#' @param fovs A character or numeric vector of FOV identifiers specifying
#'   which regions to visualize.
#' @param genes A character vector of gene names to be visualized as spatial
#'   markers.
#'
#' @return A `ggplot` object showing gene expression in the selected FOVs,
#'   with molecule-level detail and bounding box framing the target area.
#'
#' @export
render_transcript_map <- function(sobj, fovs, genes) {
  box::use(
    ggplot2[guides, xlim, ylim],
    Seurat[ImageDimPlot, Images]
  )

  fov_boundaries <- get_spatial_boundaries(sobj, fovs)
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
    ylim(fov_boundaries[, 1]) +
    guides(fill = "none") # Remvoes cell type legend
}

#' Get custom colors for plotting
#'
#' @description
#' Selects a subset of the `custom_categorical` palette based on the number of
#' groups in the input data. If the number of groups exceeds the palette length,
#' an error is thrown.
#'
#' @param data A `data.frame` or `tibble` containing the group column.
#' @param group A `character` string indicating the name of the grouping column.
#'
#' @return A `character` vector of hexadecimal color codes.
#'
#' @keywords internal
get_custom_colors <- function(data, group) {
  box::use(
    cli[cli_abort],
    dplyr[n_distinct],
    rlang[caller_env]
  )

  # Determine number of groups
  n_groups <- n_distinct(data[[group]])

  # Number of groups to visualize must fit custom palette
  # Currently 12 custom colors are supported
  if (n_groups > length(custom_categorical)) {
    cli_abort(
      c(
        "Groups must not exceed {.val {length(custom_categorical)}}.",
        "i" = "You provided {.val {n_groups}} groups to visualize."
      ),
      call = caller_env()
    )
  }

  # Assign custom colors according to number of groups
  custom_colors <- custom_categorical[1:n_groups]

  custom_colors
}

#' Render density plot
#'
#' @description
#' Creates a kernel density plot for one or more distributions grouped by a
#' categorical variable.
#'
#' @details
#' The input data must contain two required columns:
#'
#' - `value`: Numeric values to be visualized.
#' - group: A categorical variable indicating group membership.
#'
#' The plot uses a semi-transparent fill to allow overlapping distributions to
#' remain visible.
#'
#' @param data A `data.frame` or `tibble` containing `value` and grouping column.
#' @param group An unquoted column name used for grouping the distributions.
#' @param legend `character` string for the legend title. Defaults to `"Legend"`.
#' @param title `character` string for the plot title. Defaults to `"Density Plot"`.
#' @param x_lab `character` string for the x-axis label. Defaults to `"Value"`.
#' @param y_lab `character` string for the y-axis label. Defaults to `"Density"`.
#'
#' @return A `ggplot` object representing the density plot.
#'
#' @export
render_density_plot <- function(data,
                                group,
                                legend = "Legend",
                                title = "Density Plot",
                                x_lab = "Value",
                                y_lab = "Density") {
  box::use(
    ggplot2,
    rlang[as_string, ensym]
  )

  # Convert to symbol
  group <- ensym(group)

  # Get custom colors for number of groups to plot
  custom_colors <- get_custom_colors(data, as_string(group))

  # Plot density
  plot <- data |>
    ggplot2$ggplot(ggplot2$aes(x = value, fill = !!group)) +
    ggplot2$geom_density(alpha = 0.5) + # 0.5 to make all distributions visible
    ggplot2$labs(title = title, x = x_lab, y = y_lab) +
    ggplot2$scale_fill_manual(values = custom_colors, name = legend) +
    ggplot2$theme_minimal()

  plot
}

#' Render violin plot
#'
#' @description
#' Produces a violin plot overlaid with a box plot to display the distribution
#' and summary statistics of one or more groups.
#'
#' @details
#' The input data must contain the following columns:
#'
#' - `value`: Numeric values for each observation.
#' - group: Categorical variable indicating group membership.
#'
#' The function visualizes the shape and spread of each group using a violin
#' plot, with a centered boxplot showing the interquartile range.
#'
#' @param data A `data.frame` or `tibble` containing the required columns.
#' @param group An unquoted column name used for grouping.
#' @param legend `character` string for the legend title. Defaults to `"Legend"`.
#' @param title `character` string for the plot title. Defaults to `"Violin Plot"`.
#' @param x_lab `character` string for the x-axis label. Defaults to `"Group"`.
#' @param y_lab `character` string for the y-axis label. Defaults to `"Value"`.
#'
#' @return A `ggplot` object representing the violin plot.
#'
#' @export
render_violin_plot <- function(data,
                               group,
                               legend = "Legend",
                               title = "Violin Plot",
                               x_lab = "Group",
                               y_lab = "Value") {
  box::use(
    ggplot2,
    rlang[as_string, ensym]
  )

  # Convert to symbol
  group <- ensym(group)

  # Get custom colors for number of groups to plot
  custom_colors <- get_custom_colors(data, as_string(group))

  # Plot violins
  plot <- data |>
    ggplot2$ggplot(ggplot2$aes(x = !!group, y = value, fill = !!group)) +
    ggplot2$geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + # Draw quartiles
    ggplot2$geom_boxplot(width = 0.1, fill = "#ffffff") + # Color box plot
    ggplot2$labs(title = title, x = x_lab, y = y_lab) +
    ggplot2$scale_fill_manual(values = custom_colors, name = legend) +
    ggplot2$theme_minimal()

  plot
}

#' Render ECDF plot
#'
#' @description
#' Generates an empirical cumulative distribution function (ECDF) plot to
#' compare multiple distributions using step lines.
#'
#' @details
#' The input data must contain:
#'
#' - `value`: Numeric values to be plotted.
#' - group: Categorical grouping variable.
#'
#' Each group's ECDF is plotted as a separate line to show the cumulative probability.
#'
#' @param data A `data.frame` or `tibble` with columns `value` and `group`.
#' @param group Unquoted name of the column used for grouping.
#' @param legend `character` string for the legend title. Defaults to `"Legend"`.
#' @param title `character` string for the plot title. Defaults to `"ECDF Plot"`.
#' @param x_lab `character` string for the x-axis label. Defaults to `"Value"`.
#' @param y_lab `character` string for the y-axis label. Defaults to `"Cumulative Probability"`.
#'
#' @return A `ggplot` object representing the ECDF plot.
#'
#' @export
render_ecdf_plot <- function(data,
                             group,
                             legend = "Legend",
                             title = "ECDF Plot",
                             x_lab = "Value",
                             y_lab = "Cumulative Probability") {
  box::use(
    ggplot2,
    rlang[as_string, ensym]
  )

  # Convert to symbol
  group <- ensym(group)

  # Get custom colors for number of groups to plot
  custom_colors <- get_custom_colors(data, as_string(group))

  # Plot cumulative proportion
  plot <- data |>
    ggplot2$ggplot(ggplot2$aes(x = value, color = !!group)) +
    ggplot2$stat_ecdf(geom = "step") +
    ggplot2$labs(title = title, x = x_lab, y = y_lab) +
    ggplot2$scale_color_manual(values = custom_colors, name = legend) +
    ggplot2$theme_minimal()

  plot
}

#' Render UMAP plot
#'
#' @description
#' Creates a scatter plot of UMAP embeddings colored by a categorical variable.
#'
#' @details
#' The input data must contain:
#'
#' - `umap1`, `umap2`: Numeric coordinates from UMAP dimensionality reduction.
#' - color: A categorical variable to group points.
#'
#' If the number of unique groups exceeds the palette size, the legend is
#' disabled.
#'
#' @param data A `data.frame` or `tibble` with UMAP coordinates and a color column.
#' @param color Unquoted name of the grouping column.
#' @param legend `character` string for the legend title. Defaults to `"Legend"`.
#' @param title `character` string for the plot title. Defaults to `"UMAP Plot"`.
#' @param x_lab `character` string for the x-axis label. Defaults to `"UMAP1"`.
#' @param y_lab `character` string for the y-axis label. Defaults to `"UMAP2"`.
#'
#' @return A `ggplot` object visualizing the UMAP results.
#'
#' @export
render_umap_plot <- function(data,
                             color,
                             legend = "Legend",
                             title = "UMAP Plot",
                             x_lab = "UMAP1",
                             y_lab = "UMAP2") {
  box::use(
    cli[cli_abort],
    dplyr[n_distinct],
    ggplot2,
    rlang[as_string, caller_env, ensym],
  )

  # Convert to symbol
  color <- ensym(color)

  # Determine number of groups
  n_groups <- n_distinct(data[[as_string(color)]])

  # Plot UMAP
  plot <- data |>
    ggplot2$ggplot(ggplot2$aes(x = umap1, y = umap2, color = !!color)) +
    ggplot2$geom_point() +
    ggplot2$labs(title = title, x = x_lab, y = y_lab) +
    ggplot2$theme_minimal()

  # Number of groups can be greater than available custom colors
  # Currently 12 custom colors are supported
  if (n_groups <= length(custom_categorical)) {
    # Assign custom colors according to number of groups
    custom_colors <- custom_categorical[1:n_groups]

    plot <- plot + ggplot2$scale_color_manual(values = custom_colors, name = legend)
  } else {
    plot <- plot + ggplot2$guides(color = "none")
  }

  plot
}

#' Render bar plot
#'
#' @description
#' Draws a horizontal bar plot where bars are colored by a categorical variable.
#'
#' @details
#' The input data must contain:
#'
#' - `value`: A numeric variable to define bar lengths.
#' - `group`: A factor or character variable mapped to the y-axis.
#' - color: A categorical variable used to fill the bars.
#'
#' The function supports up to 12 unique fill groups based on the custom color
#' palette.
#'
#' @param data A `data.frame` or `tibble` containing columns `value`, `group`, and `color`.
#' @param color Unquoted column name used for fill color.
#' @param legend `character` string for the legend title. Defaults to `"Legend"`.
#' @param title `character` string for the plot title. Defaults to `"Bar Plot"`.
#' @param x_lab `character` string for the x-axis label. Defaults to `"Value"`.
#' @param y_lab `character` string for the y-axis label. Defaults to `"Group"`.
#'
#' @return A `ggplot` object representing the bar plot.
#'
#' @export
render_bar_plot <- function(data,
                            color,
                            legend = "Legend",
                            title = "Bar Plot",
                            x_lab = "Value",
                            y_lab = "Group") {
  box::use(
    ggplot2,
    rlang[as_string, ensym]
  )

  # Convert to symbol
  group <- ensym(group)

  # Get custom colors for number of groups to plot
  custom_colors <- get_custom_colors(data, as_string(group))

  plot <- data |>
    ggplot2$ggplot(ggplot2$aes(x = value, y = group, fill = !!color)) +
    ggplot2$geom_col(position = "identity", alpha = 0.75) +
    ggplot2$labs(title = title, x = x_lab, y = y_lab) +
    ggplot2$scale_fill_manual(values = custom_colors, name = legend) +
    ggplot2$theme_minimal()
}
