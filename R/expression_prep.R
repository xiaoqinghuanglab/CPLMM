#' Prepare long-format expression data across strata and diagnostic categories
#'
#' Takes three cohort-specific data frames (normal-only, status-change,
#' abnormal-only) plus the full data frame with a diagnosis/category column.
#' Selects a subset of genes/proteins, pivots them to long format, coerces
#' values to numeric, tags the source, and binds all results together.
#'
#' @param df_normal_only Data frame for normal-only subjects.
#' @param df_status_change Data frame for status-change subjects.
#' @param df_abnormal_only Data frame for abnormal-only subjects.
#' @param df_all Full data frame containing the `category_col`.
#' @param subset_genes Character vector of gene/protein column names to include.
#' @param category_col Column name in `df_all` with diagnosis categories.
#'   Default: `"CATEGORY"`.
#' @param categories Character vector of category labels to include from `df_all`.
#'   Default: `c("SCD","MCI","AD Dementia","FTD Dementia")`.
#' @param normal_label Label to use for the normal-only cohort in the output.
#'   Default: `"Normal_only"`.
#' @param status_label Label to use for the status-change cohort in the output.
#'   Default: `"Status_change"`.
#' @param abnormal_label Label to use for the abnormal-only cohort in the output.
#'   Default: `"Abnormal_only"`.
#'
#' @return A tibble with columns `Gene`, `Expression`, `Source`, where
#'   `Gene` is an ordered factor with levels matching `subset_genes`.
#' @export
prepare_combined_expression <- function(
    df_normal_only,
    df_status_change,
    df_abnormal_only,
    df_all,
    subset_genes,
    category_col = "CATEGORY",
    categories = c("SCD", "MCI", "AD Dementia", "FTD Dementia"),
    normal_label  = "Normal_only",
    status_label  = "Status_change",
    abnormal_label = "Abnormal_only"
) {
  # deps
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("rlang", quietly = TRUE)

  # basic checks similar to pandas KeyError if missing columns
  needed <- unique(subset_genes)
  miss_n <- setdiff(needed, names(df_normal_only))
  miss_s <- setdiff(needed, names(df_status_change))
  miss_a <- setdiff(needed, names(df_abnormal_only))
  miss_all <- setdiff(c(needed, category_col), names(df_all))
  if (length(miss_n) || length(miss_s) || length(miss_a) || length(miss_all)) {
    stop("Missing required columns in inputs:\n",
         if (length(miss_n))  paste0(" - df_normal_only: ", paste(miss_n, collapse = ", "), "\n"),
         if (length(miss_s))  paste0(" - df_status_change: ", paste(miss_s, collapse = ", "), "\n"),
         if (length(miss_a))  paste0(" - df_abnormal_only: ", paste(miss_a, collapse = ", "), "\n"),
         if (length(miss_all))paste0(" - df_all: ", paste(miss_all, collapse = ", "), "\n"),
         call. = FALSE)
  }

  melt_subset <- function(df, label) {
    df |>
      dplyr::select(dplyr::all_of(subset_genes)) |>
      tidyr::pivot_longer(dplyr::everything(), names_to = "Gene", values_to = "Expression") |>
      dplyr::mutate(
        Expression = suppressWarnings(as.numeric(.data$Expression)),
        Source = label
      )
  }

  # core groups
  melted_normal   <- melt_subset(df_normal_only, normal_label)
  melted_status   <- melt_subset(df_status_change, status_label)
  melted_abnormal <- melt_subset(df_abnormal_only, abnormal_label)

  # diagnostic categories from df_all
  category_melted <- lapply(categories, function(cat) {
    df_all |>
      dplyr::filter(.data[[category_col]] == cat) |>
      dplyr::select(dplyr::all_of(subset_genes)) |>
      tidyr::pivot_longer(dplyr::everything(), names_to = "Gene", values_to = "Expression") |>
      dplyr::mutate(
        Expression = suppressWarnings(as.numeric(.data$Expression)),
        Source = gsub(" ", "_", cat, fixed = TRUE)
      )
  })

  combined_expr <- dplyr::bind_rows(
    melted_normal, melted_status, melted_abnormal, dplyr::bind_rows(category_melted)
  ) |>
    dplyr::mutate(Gene = factor(.data$Gene, levels = subset_genes, ordered = TRUE))

  tibble::as_tibble(combined_expr)
}
