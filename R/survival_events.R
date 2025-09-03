#' Compute event/censoring times from longitudinal biomarker threshold crossing
#'
#' For each subject, sorts by \code{PROCEDURE_AGE}, takes the first value of
#' \code{onset_source} as the onset reference, and checks whether the biomarker
#' ever crosses \code{threshold}. If yes, event time = age at the first crossing
#' minus onset; if not, censor at last observed \code{PROCEDURE_AGE} minus onset.
#'
#' @param df A data.frame/tibble of longitudinal visits with columns
#'   \code{SUBID}, \code{PROCEDURE_AGE}, the biomarker column, and \code{onset_source}.
#' @param threshold Numeric threshold for event definition.
#' @param biomarker Character; biomarker column name.
#' @param group_label Character; label stored in the \code{group} column of output.
#' @param onset_source Character; column giving onset-age reference. Default \code{"ONSET_AGE"}.
#'
#' @return A tibble with columns \code{SUBID}, \code{time}, \code{event} (0/1), \code{group}.
#' @export
compute_event_df <- function(df,
                             threshold,
                             biomarker,
                             group_label,
                             onset_source = "ONSET_AGE") {
  # basic checks
  need <- c("SUBID", "PROCEDURE_AGE", biomarker, onset_source)
  miss <- setdiff(need, names(df))
  if (length(miss)) {
    stop("compute_event_df(): missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  # split-apply in base to avoid version-specific dplyr pitfalls
  split_list <- split(df, df$SUBID)
  out_rows <- vector("list", length(split_list))
  i <- 0L

  for (subid in names(split_list)) {
    g <- split_list[[subid]]
    # order by PROCEDURE_AGE
    ord <- order(g$PROCEDURE_AGE)
    g <- g[ord, , drop = FALSE]

    onset_age <- g[[onset_source]][1]
    # guard against all-NA onset
    if (is.na(onset_age)) {
      # if onset missing entirely, skip or set event=0, time=NA; we'll choose NA time + event 0
      time <- NA_real_
      event <- 0L
    } else {
      # indices where biomarker crosses threshold
      idx <- which(suppressWarnings(as.numeric(g[[biomarker]]) >= threshold))
      if (length(idx) > 0L) {
        time <- g$PROCEDURE_AGE[min(idx)] - onset_age
        event <- 1L
      } else {
        time <- max(g$PROCEDURE_AGE, na.rm = TRUE) - onset_age
        event <- 0L
      }
    }

    i <- i + 1L
    out_rows[[i]] <- data.frame(
      SUBID = subid,
      time = as.numeric(time),
      event = as.integer(event),
      group = group_label,
      stringsAsFactors = FALSE
    )
  }

  tibble::as_tibble(do.call(rbind, out_rows))
}
