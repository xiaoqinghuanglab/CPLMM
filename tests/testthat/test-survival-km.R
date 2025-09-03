test_that("plot_km_with_threshold returns a patchwork object", {
  set.seed(1)
  mk <- function(n_id, cat) {
    n <- n_id * 3
    df <- data.frame(
      SUBID = rep(paste0(cat, 1:n_id), each = 3),
      PROCEDURE_AGE = rep(60:62, times = n_id),
      ONSET_AGE = rep(59, n),
      CATEGORY = cat,
      NFL = rnorm(n, if (cat == "MCI") 15 else 8, 3)
    )
    df
  }
  wd <- rbind(mk(5, "MCI"), mk(5, "Other"))
  norm <- mk(4, "Normal")
  abnm <- mk(3, "Abn")  # not used in plot, passed for parity

  p <- plot_km_with_threshold(
    biomarker_name = "NFL",
    threshold = 12,
    wd_df = wd,
    normal_df = norm,
    abnm_df = abnm,
    time_points = seq(-2, 4, by = 2),
    save = FALSE
  )
  expect_s3_class(p, "gg")
})
