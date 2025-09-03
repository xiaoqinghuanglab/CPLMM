#' Mann–Whitney (Wilcoxon rank-sum) comparison for two groups across genes
#'
#' Given a long-format expression table with columns `Gene`, `Expression`, `Source`,
#' performs a two-sided Wilcoxon rank-sum test (aka Mann–Whitney U) for each gene
#' comparing `group2` vs `group1`. It reports group means, their difference
#' (`Delta_mean = mean(group2) - mean(group1)`), test statistic, raw p-value,
#' BH/FDR-adjusted p-value, and a significance flag at `alpha`.
#'
#' @param combined_expr A data.frame/tibble with columns `Gene`, `Expression`, `Source`.
#' @param gene_list Character vector of gene/protein names to test.
#' @param group1 Character; label of the first group in `Source`.
#' @param group2 Character; label of the second group in `Source`.
#' @param alpha Numeric; FDR threshold for `Significant`. Default 0.05.
#' @param rank_by Character; one of `"FDR"`, `"p_value"`, or `"Delta_mean"`.
#'   Sorting is ascending for p/FDR, descending for Delta_mean (mirrors Python).
#'
#' @return A tibble with columns:
#'   `Gene`, `<group1>_mean`, `<group2>_mean`, `Delta_mean`,
#'   `U_statistic`, `p_value`, `FDR`, `Significant`, and sorted by `rank_by`.
#' @export
compare_groups_mannwhitney <- function(
    combined_expr,
    gene_list,
    group1,
    group2,
    alpha = 0.05,
    rank_by = "FDR"
) {
  stopifnot(all(c("Gene","Expression","Source") %in% names(combined_expr)))
  df <- combined_expr

  rows <- list()
  k <- 0L

  col1 <- paste0(group1, "_mean")
  col2 <- paste0(group2, "_mean")

  for (g in gene_list) {
    expr1 <- df$Expression[df$Gene == g & df$Source == group1]
    expr2 <- df$Expression[df$Gene == g & df$Source == group2]
    expr1 <- expr1[!is.na(expr1)]
    expr2 <- expr2[!is.na(expr2)]

    if (length(expr1) > 0L && length(expr2) > 0L) {
      wt <- tryCatch(
        stats::wilcox.test(expr1, expr2, alternative = "two.sided", exact = FALSE),
        error = function(e) NULL
      )
      if (!is.null(wt)) {
        k <- k + 1L
        m1 <- mean(expr1)
        m2 <- mean(expr2)
        # Build a single-row data.frame with dynamic column names
        row <- data.frame(
          Gene = g,
          check.names = FALSE, stringsAsFactors = FALSE
        )
        row[[col1]] <- m1
        row[[col2]] <- m2
        row[["Delta_mean"]]  <- m2 - m1
        row[["U_statistic"]] <- as.numeric(wt$statistic)
        row[["p_value"]]     <- wt$p.value

        # Reorder columns for consistency
        row <- row[, c("Gene", col1, col2, "Delta_mean", "U_statistic", "p_value"), drop = FALSE]
        rows[[k]] <- row
      }
    }
  }

  if (k == 0L) {
    out <- data.frame(Gene = character(), check.names = FALSE, stringsAsFactors = FALSE)
    out[[col1]] <- numeric()
    out[[col2]] <- numeric()
    out[["Delta_mean"]]  <- numeric()
    out[["U_statistic"]] <- numeric()
    out[["p_value"]]     <- numeric()
    out[["FDR"]]         <- numeric()
    out[["Significant"]] <- logical()
    return(tibble::as_tibble(out))
  }

  out <- do.call(rbind, rows)

  # FDR (BH)
  out$FDR <- stats::p.adjust(out$p_value, method = "BH")
  out$Significant <- out$FDR < alpha

  # sorting
  if (rank_by %in% names(out)) {
    if (rank_by == "Delta_mean") {
      out <- out[order(-out$Delta_mean), , drop = FALSE]
    } else {
      out <- out[order(out[[rank_by]]), , drop = FALSE]
    }
  } else {
    out <- out[order(out$FDR), , drop = FALSE]
  }

  rownames(out) <- NULL
  tibble::as_tibble(out)
}
