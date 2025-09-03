test_that("plot_pathway_bubble returns a ggplot", {
  df <- tibble::tibble(
    Cleaned_Pathway = rep(c("Axon guidance","Synapse","Immune"), each = 6),
    BioCategory_Manual = rep(c("Neuro","Neuro","Immune"), each = 6),
    Source = rep(c("DAVID","Reactome"), times = 9),
    LogQValue = runif(18, 1, 6),
    Gene = paste0("G", sample(LETTERS, 18, TRUE))
  )
  p <- plot_pathway_bubble(df, export = FALSE)
  expect_s3_class(p, "ggplot")
})
