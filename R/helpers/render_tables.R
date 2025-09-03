#' Render styled table
#'
#' @description
#' `render_table()` returns a styled HTML table using the `gt` package.
#' The output is useful for rendering clean and readable tables in R
#' Markdown documents or Shiny applications.
#'
#' @details
#' This function applies a consistent table theme:
#'
#' * Thick top, heading, and bottom borders to frame the table.
#' * Bold column labels.
#' * No horizontal lines between rows.
#' * Left-aligned content for all columns.
#'
#' @param data A `data.frame` or `tibble` to be rendered as a `gt` table.
#'
#' @return A `gt_tbl` object
#'
#' @export
render_table <- function(data) {
  box::use(
    gt
  )

  table <- data |>
    gt$gt() |>
    gt$tab_options(
      table.border.top.width = 2, # Thick top border
      heading.border.bottom.width = 2, # Thick bottom border for heading
      table.border.bottom.width = 2, # Thick bottom border
      table_body.hlines.width = 0 # No horizontal lines between rows
    ) |>
    gt$tab_style(
      style = gt$cell_text(weight = "bold"),
      locations = gt$cells_column_labels(gt$everything()) # Bold column titles
    ) |>
    gt$cols_align(
      align = "left",
      columns = gt$everything() # Left-align all columns
    )

  table
}
