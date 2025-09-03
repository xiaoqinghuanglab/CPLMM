#' Volcano plot for Wald test results
#'
#' Builds a volcano using signed \eqn{\log_2} fold-change between post- and pre-onset
#' slopes (\code{beta_after_col / beta_before_col}) and \eqn{-\log_{10}(p)}.
#' Points are colored by significance: \code{p < pval_threshold} \emph{and}
#' \code{FDR < fdr_threshold}. You can optionally annotate significant proteins
#' (or a custom list) and export in multiple formats.
#'
#' @param wald_df Data frame with Wald results (from \code{compute_wald_test()} or similar).
#' @param pval_col Column with raw p-values. Default \code{"P-value 1"}.
#' @param fdr_col Column with FDR-adjusted p-values. Default
#'   \code{"Adjusted P-value (FDR)_1"}. If not found, common alternatives are tried
#'   (e.g. \code{"Adjusted P-value 1"}, \code{"FDR"}).
#' @param beta_before_col Column with slope before onset (β1). Default \code{"Beta 1"}.
#' @param beta_after_col Column with slope after onset (β3). Default \code{"Beta 3"}.
#' @param protein_col Column with protein/gene names. Default \code{"Protein"}.
#' @param pval_threshold Raw p-value cutoff for significance. Default \code{0.05}.
#' @param fdr_threshold FDR cutoff for significance. Default \code{0.05}.
#' @param annotate Logical; annotate points? Default \code{TRUE}.
#' @param annotate_list Optional character vector of proteins to annotate; if \code{NULL}
#'   and \code{annotate=TRUE}, annotates all points labeled "Significant".
#' @param style_config Optional list of \pkg{ggplot2} theme overrides (e.g.,
#'   \code{list(legend.position = "bottom")}).
#' @param export Logical; save the plot? Default \code{FALSE}.
#' @param export_dir Directory to save into. Default \code{"Figures"}.
#' @param export_name Base filename (no extension). Default \code{"wald_volcano"}.
#' @param export_formats Character vector of extensions, e.g. \code{c("pdf","svg")}.
#'
#' @return A \pkg{ggplot2} object.
#' @export
plot_wald_volcano <- function(
    wald_df,
    pval_col = "P-value 1",
    fdr_col  = "Adjusted P-value (FDR)_1",
    beta_before_col = "Beta 1",
    beta_after_col  = "Beta 3",
    protein_col = "Protein",
    pval_threshold = 0.05,
    fdr_threshold  = 0.05,
    annotate = TRUE,
    annotate_list = NULL,
    style_config = NULL,
    export = FALSE,
    export_dir = "Figures",
    export_name = "wald_volcano",
    export_formats = c("pdf","svg")
) {
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  df <- tibble::as_tibble(wald_df)

  # Resolve columns (be forgiving about fdr col naming)
  if (!pval_col %in% names(df)) {
    stop("plot_wald_volcano(): column '", pval_col, "' not found.", call. = FALSE)
  }
  if (!fdr_col %in% names(df)) {
    candidates <- c("Adjusted P-value 1", "Adjusted P-value 2", "FDR", "FDR_1", "Adjusted P-value (FDR)_1")
    alt <- intersect(candidates, names(df))
    if (length(alt) == 0) {
      stop("plot_wald_volcano(): FDR column '", fdr_col,
           "' not found, and no common alternatives present.", call. = FALSE)
    }
    fdr_col <- alt[1]
  }
  need <- c(beta_before_col, beta_after_col, protein_col)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("plot_wald_volcano(): missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # Derived metrics
  df$neglog10p <- -log10(suppressWarnings(as.numeric(df[[pval_col]])))
  fc <- suppressWarnings(as.numeric(df[[beta_after_col]]) / as.numeric(df[[beta_before_col]]))
  df$log2FoldChange <- sign(fc) * log2(abs(fc))
  df$Significance <- ifelse(df[[pval_col]] < pval_threshold & df[[fdr_col]] < fdr_threshold,
                            "Significant", "Not Significant")

  # Keep only finite values for plotting
  df <- df[is.finite(df$neglog10p) & is.finite(df$log2FoldChange), , drop = FALSE]

  # Base plot
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$log2FoldChange, y = .data$neglog10p, fill = .data$Significance)
  ) +
    ggplot2::geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "steelblue", linewidth = 0.5) +
    ggplot2::geom_point(shape = 21, color = "black", stroke = 0.3, alpha = 0.7, size = 1.8) +
    ggplot2::scale_fill_manual(values = c("Not Significant" = "gray", "Significant" = "#B24745")) +
    ggplot2::labs(
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      title = "Volcano Plot for Wald Test",
      fill = "Significance"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "#DDDDDD", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 10),
      axis.text  = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.background = ggplot2::element_rect(fill = "white", color = "grey80")
    )

  # Optional annotations
  if (isTRUE(annotate)) {
    if (!is.null(annotate_list)) {
      lab_df <- df[df[[protein_col]] %in% annotate_list, , drop = FALSE]
    } else {
      lab_df <- df[df$Significance == "Significant", , drop = FALSE]
    }
    if (nrow(lab_df)) {
      p <- p + ggplot2::geom_text(
        data = lab_df,
        ggplot2::aes(label = .data[[protein_col]]),
        size = 2.8, vjust = -0.2, hjust = 1, color = "black"
      )
    }
  }

  # Style overrides, if provided
  if (!is.null(style_config) && length(style_config)) {
    p <- p + do.call(ggplot2::theme, style_config)
  }

  # Export if requested
  if (isTRUE(export)) {
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    for (fmt in export_formats) {
      outfile <- file.path(export_dir, paste0(export_name, ".", fmt))
      ggplot2::ggsave(outfile, p, width = 7, height = 5, dpi = 300)
    }
  }

  p
}
