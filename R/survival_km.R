#' KM plot for biomarker threshold crossing (with log-rank and risk table)
#'
#' Builds event data via \code{compute_event_df()} for three groups:
#' Status-Change (all \code{wd_df}), Normal-only (\code{normal_df}), and MCI
#' (subset of \code{wd_df} where \code{CATEGORY == "MCI"}). Fits KM curves,
#' draws a vertical line at onset (0), annotates the Normal vs MCI log-rank
#' p-value, and shows an at-risk table. Internally, times are shifted so they're
#' non-negative (required by \pkg{survival}), but the x-axis shows the original
#' requested time points.
#'
#' @param biomarker_name Character; biomarker column (e.g. "NFL").
#' @param threshold Numeric; threshold for event definition.
#' @param wd_df Data frame of Status-Change subjects (must contain SUBID, PROCEDURE_AGE, ONSET_AGE, CATEGORY, biomarker).
#' @param normal_df Data frame of Normal-only subjects.
#' @param abnm_df Data frame of Abnormal-only subjects (kept for parity; not plotted by default).
#' @param threshold_dict Optional; ignored (kept for parity with Python signature).
#' @param save Logical; save the combined figure as SVG? Default FALSE.
#' @param save_path Directory to save into if \code{save=TRUE}. Default "Figures".
#' @param jama_palette Named character vector of colors; names should include
#'   "Status-Change", "Normal", "MCI". Defaults to JAMA-like palette.
#' @param time_points Numeric vector of x tick labels (in years to onset) to show.
#'   Default: \code{seq(-10, 4, by = 2)} (matches Python's -10..+4).
#'
#' @return A patchwork/ggplot object with the KM and risk-table stacked.
#' @export
plot_km_with_threshold <- function(
    biomarker_name,
    threshold,
    wd_df,
    normal_df,
    abnm_df,
    threshold_dict = NULL,
    save = FALSE,
    save_path = "Figures",
    jama_palette = NULL,
    time_points = seq(-10, 4, by = 2)
) {
  # deps
  requireNamespace("survival", quietly = TRUE)
  requireNamespace("survminer", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("patchwork", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # ---- 1) Build event data (Status-Change, Normal, MCI) ----
  if (!("CATEGORY" %in% names(wd_df))) {
    stop("wd_df must contain a CATEGORY column (to extract MCI subgroup).", call. = FALSE)
  }
  mci_ids <- unique(wd_df$SUBID[wd_df$CATEGORY == "MCI"])
  mci_df  <- wd_df[wd_df$SUBID %in% mci_ids, , drop = FALSE]

  ev_normal <- compute_event_df(normal_df, threshold, biomarker_name, "Normal")
  ev_mci    <- compute_event_df(mci_df,    threshold, biomarker_name, "MCI")
  ev_status <- compute_event_df(wd_df,     threshold, biomarker_name, "Status-Change")

  event_df <- dplyr::bind_rows(ev_status, ev_normal, ev_mci)
  if (!nrow(event_df)) stop("No event data could be computed.", call. = FALSE)

  # enforce group order to match palette order
  event_df$group <- factor(event_df$group, levels = c("Status-Change","Normal","MCI"))

  # ---- 2) Time shifting for survival (R requires non-negative time) ----
  # We shift all times so min time >= 0, but label axis with original values.
  min_time_all <- suppressWarnings(min(c(event_df$time, time_points), na.rm = TRUE))
  offset <- if (is.finite(min_time_all) && min_time_all < 0) -min_time_all else 0
  event_df$time_shifted <- event_df$time + offset
  ticks_shifted <- time_points + offset

  # ---- 3) Fit KM ----
  sf <- survival::survfit(survival::Surv(time_shifted, event) ~ group, data = event_df)

  # default palette
  if (is.null(jama_palette)) {
    jama_palette <- c("Status-Change" = "#DF8F44", "Normal" = "#00A1D5", "MCI" = "#B24745")
  }
  # survminer expects a vector of colors in the order of strata
  strata_names <- names(sf$strata)
  # e.g., "group=Normal" -> "Normal"
  strata_levels <- sub("^group=", "", strata_names)
  pal_vec <- unname(jama_palette[strata_levels])

  # break spacing for risk table ticks if evenly spaced
  equal_step <- FALSE
  if (length(time_points) >= 2) {
    diffs <- diff(time_points)
    equal_step <- all(abs(diffs - diffs[1]) < 1e-9)
  }
  break_by <- if (equal_step) abs(diff(time_points)[1]) else NULL

  g <- survminer::ggsurvplot(
    sf,
    data = event_df,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text.col = TRUE,
    break.time.by = break_by,    # if NULL, survminer picks breaks
    xlim = range(ticks_shifted, na.rm = TRUE),
    palette = pal_vec,
    legend.title = "",
    legend = "top",
    ggtheme = ggplot2::theme_minimal(base_size = 12)
  )

  # ---- 4) Customize axes, onset line, titles, labels ----
  x_min <- min(ticks_shifted, na.rm = TRUE)
  x_max <- max(ticks_shifted, na.rm = TRUE)
  x_pos <- x_min + 0.65 * (x_max - x_min)
  y_pos <- 0.1

  g$plot <- g$plot +
    ggplot2::geom_vline(xintercept = offset, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    ggplot2::scale_x_continuous(breaks = ticks_shifted, labels = time_points, limits = c(x_min, x_max)) +
    ggplot2::labs(
      title = sprintf("%s \u2265 %.2f", biomarker_name, threshold),  # ≥
      x = "Years to Onset",
      y = "Proportion Not Yet Crossed"
    ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(linetype = "dotted", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )

  # Sync x scale for risk table and label ticks with original (pre-shift) labels
  if (!is.null(g$table)) {
    g$table <- g$table +
      ggplot2::scale_x_continuous(breaks = ticks_shifted, labels = time_points, limits = c(x_min, x_max)) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  # ---- 5) Log-rank test (Normal vs MCI) and annotate p-value ----
  has_norm <- any(event_df$group == "Normal")
  has_mci  <- any(event_df$group == "MCI")
  if (has_norm && has_mci) {
    sub_nm <- event_df[event_df$group %in% c("Normal","MCI"), , drop = FALSE]
    sub_nm$grp2 <- droplevels(sub_nm$group)
    sd <- survival::survdiff(survival::Surv(time_shifted, event) ~ grp2, data = sub_nm)
    # df = number of groups - 1 = 1
    pval <- stats::pchisq(sd$chisq, df = (length(sd$n) - 1), lower.tail = FALSE)
    g$plot <- g$plot +
      ggplot2::annotate("text", x = x_pos, y = y_pos,
                        label = sprintf("p = %.3e", pval),
                        size = 3.8, hjust = 0)
  }

  # ---- 6) Stack KM + risk table and optionally save ----
  combined <- if (is.null(g$table)) {
    g$plot
  } else {
    g$plot / g$table + patchwork::plot_layout(heights = c(3, 1))
  }

  if (isTRUE(save)) {
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    outfile <- file.path(save_path, sprintf("KM_%s_JAMA.svg", biomarker_name))
    ggplot2::ggsave(outfile, combined, width = 6, height = 8, dpi = 300)
  }

  combined
}
