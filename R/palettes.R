#' JAMA palette for diagnostic categories (Python parity)
#'
#' Returns a named character vector of hex colors matching your Python
#' \code{JAMA_PALETTE} (keys use underscores).
#'
#' @return Named character vector: Normal_only, SCD, MCI, AD_Dementia, FTD_Dementia.
#' @examples
#' pal <- jama_palette()
#' pal["MCI"]
#' @export
jama_palette <- function() {
  c(
    Normal_only  = "#00A1D5",
    SCD          = "#DF8F44",
    MCI          = "#B24745",
    AD_Dementia  = "#79AF97",
    FTD_Dementia = "#6A6599"
  )
}

#' ggplot2 scales using the JAMA diagnostic palette
#' @inheritParams ggplot2::scale_color_manual
#' @export
scale_color_jama_diagnosis <- function(...) {
  ggplot2::scale_color_manual(values = jama_palette(), ...)
}

#' @rdname scale_color_jama_diagnosis
#' @export
scale_fill_jama_diagnosis <- function(...) {
  ggplot2::scale_fill_manual(values = jama_palette(), ...)
}
