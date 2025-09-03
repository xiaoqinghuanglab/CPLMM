test_that("plot_top_pathways_bar returns a ggplot", {
  set.seed(1)
  df <- tibble::tibble(
    Cleaned_Pathway = rep(paste("Path", 1:10), each = 6),
    Gene            = sample(paste0("G", 1:40), 60, replace = TRUE),
    LogQValue       = runif(60, 1, 6),
    BioCategory_Manual = rep(c("Neuro","Immune","Glia","Metabolism","Other"), length.out = 60)
  )
  p <- plot_top_pathways_bar(df, top_n = 8, annotate = TRUE)
  expect_s3_class(p, "ggplot")
})
