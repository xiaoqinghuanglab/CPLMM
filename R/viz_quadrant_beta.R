#' Quadrant plot of beta coefficients from Wald results
#'
#' Plots \code{beta_x_col} (x-axis) vs \code{beta_y_col} (y-axis) with axes at 0.
#' Points are split by FDR significance, significant ones highlighted and (optionally)
#' annotated. Corner labels show counts per quadrant.
#'
#' @param wald_df Data frame of Wald results.
#' @param beta_x_col Column for x beta (default "Beta 1").
#' @param beta_y_col Column for y beta (default "Beta 3").
#' @param fdr_col Column for FDR-adjusted p-values (default "Adjusted P-value (FDR)_1";
#'   if missing, tries common alternatives like "Adjusted P-value 1" or "FDR").
#' @param protein_col Column with protein identifiers. Default "Protein".
#' @param fdr_threshold Numeric FDR cutoff. Default 0.05.
#' @param annotate Logical; annotate significant proteins? Default TRUE.
#' @param style_config Optional list of ggplot2 theme overrides.
#' @param export Logical; write plot files? Default FALSE.
#' @param export_dir Output directory. Default "Figures".
#' @param export_name Base filename (no extension). Default "quadrant_plot".
#' @param export_formats Character vector of extensions. Default c("pdf","svg").
#'
#' @return A ggplot object.
#' @export
plot_quadrant_beta <- function(
    wald_df,
    beta_x_col = "Beta 1",
    beta_y_col = "Beta 3",
    fdr_col    = "Adjusted P-value (FDR)_1",
    protein_col = "Protein",
    fdr_threshold = 0.05,
    annotate = TRUE,
    style_config = NULL,
    export = FALSE,
    export_dir = "Figures",
    export_name = "quadrant_plot",
    export_formats = c("pdf", "svg")
) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("rlang", quietly = TRUE)

  df <- tibble::as_tibble(wald_df)

  # Resolve columns (be forgiving about FDR name)
  if (!all(c(beta_x_col, beta_y_col, protein_col) %in% names(df))) {
    miss <- setdiff(c(beta_x_col, beta_y_col, protein_col), names(df))
    stop("plot_quadrant_beta(): missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!fdr_col %in% names(df)) {
    candidates <- c("Adjusted P-value 1", "FDR", "Adjusted P-value", "Adjusted P-value 2", "FDR_1")
    alt <- intersect(candidates, names(df))
    if (length(alt) == 0) {
      stop("plot_quadrant_beta(): FDR column '", fdr_col, "' not found and no common alternatives present.", call. = FALSE)
    }
    fdr_col <- alt[1]
  }

  # Numeric vectors & keep finite rows
  x <- suppressWarnings(as.numeric(df[[beta_x_col]]))
  y <- suppressWarnings(as.numeric(df[[beta_y_col]]))
  fdr <- suppressWarnings(as.numeric(df[[fdr_col]]))
  prot <- as.character(df[[protein_col]])

  keep <- is.finite(x) & is.finite(y) & is.finite(fdr)
  dfp <- tibble::tibble(
    !!beta_x_col := x[keep],
    !!beta_y_col := y[keep],
    !!fdr_col    := fdr[keep],
    !!protein_col := prot[keep]
  )

  if (!nrow(dfp)) {
    stop("plot_quadrant_beta(): no finite rows to plot.", call. = FALSE)
  }

  dfp$Significant <- dfp[[fdr_col]] < fdr_threshold

  # Quadrant counts
  q1 <- sum(dfp[[beta_x_col]] > 0 & dfp[[beta_y_col]] > 0)
  q2 <- sum(dfp[[beta_x_col]] < 0 & dfp[[beta_y_col]] > 0)
  q3 <- sum(dfp[[beta_x_col]] < 0 & dfp[[beta_y_col]] < 0)
  q4 <- sum(dfp[[beta_x_col]] > 0 & dfp[[beta_y_col]] < 0)

  # Limits with 10% padding (replicates Python min/max * 1.1 behavior)
  x_min <- min(dfp[[beta_x_col]], na.rm = TRUE)
  x_max <- max(dfp[[beta_x_col]], na.rm = TRUE)
  y_min <- min(dfp[[beta_y_col]], na.rm = TRUE)
  y_max <- max(dfp[[beta_y_col]], na.rm = TRUE)
  x_lim <- c(x_min * 1.1, x_max * 1.1)
  y_lim <- c(y_min * 1.1, y_max * 1.1)

  # Corner label positions (same heuristic as Python)
  lab_pos <- list(
    q1 = c(x = x_max * 0.7, y = y_max * 0.9),
    q2 = c(x = x_min * 0.7, y = y_max * 0.9),
    q3 = c(x = x_min * 0.7, y = y_min * 0.9),
    q4 = c(x = x_max * 0.7, y = y_min * 0.9)
  )

  # Base plot
  p <- ggplot2::ggplot() +
    # Non-significant
    ggplot2::geom_point(
      data = dfp[!dfp$Significant, , drop = FALSE],
      ggplot2::aes(x = .data[[beta_x_col]], y = .data[[beta_y_col]]),
      color = "black", alpha = 0.3, size = 1.6, stroke = 0
    ) +
    # Significant
    ggplot2::geom_point(
      data = dfp[dfp$Significant, , drop = FALSE],
      ggplot2::aes(x = .data[[beta_x_col]], y = .data[[beta_y_col]]),
      color = "#B24745", alpha = 0.8, size = 2.0, stroke = 0
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::labs(
      x = if (grepl("1", beta_x_col)) paste0(beta_x_col, " (Before Onset)") else beta_x_col,
      y = if (grepl("3", beta_y_col)) paste0(beta_y_col, " (After Onset)")  else beta_y_col,
      title = "Quadrant Plot of Beta Coefficients"
    ) +
    ggplot2::coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "#DDDDDD", linewidth = 0.5),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 10),
      axis.text  = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6)
    )

  # Annotations for significant proteins
  if (isTRUE(annotate) && any(dfp$Significant)) {
    df_lab <- dfp[dfp$Significant, , drop = FALSE]
    p <- p + ggplot2::geom_text(
      data = df_lab,
      ggplot2::aes(x = .data[[beta_x_col]], y = .data[[beta_y_col]], label = .data[[protein_col]]),
      size = 2.8, hjust = 1, vjust = 0, color = "#B24745"
    )
  }

  # Quadrant count labels
  p <- p +
    ggplot2::annotate("text", x = lab_pos$q1["x"], y = lab_pos$q1["y"], label = paste0("n=", q1), size = 3.2) +
    ggplot2::annotate("text", x = lab_pos$q2["x"], y = lab_pos$q2["y"], label = paste0("n=", q2), size = 3.2) +
    ggplot2::annotate("text", x = lab_pos$q3["x"], y = lab_pos$q3["y"], label = paste0("n=", q3), size = 3.2) +
    ggplot2::annotate("text", x = lab_pos$q4["x"], y = lab_pos$q4["y"], label = paste0("n=", q4), size = 3.2)

  # Optional theme overrides
  if (!is.null(style_config) && length(style_config)) {
    p <- p + do.call(ggplot2::theme, style_config)
  }

  if (isTRUE(export)) {
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    for (fmt in export_formats) {
      outfile <- file.path(export_dir, paste0(export_name, ".", fmt))
      ggplot2::ggsave(outfile, p, width = 6, height = 6, dpi = 300)
    }
  }

  p
}
