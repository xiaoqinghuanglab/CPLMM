test_that("plot_pathway_gene_heatmap returns a ggplot", {
  set.seed(1)
  df <- tibble::tibble(
    Cleaned_Pathway = rep(c("Axon guidance","Synapse","Immune signaling","Myelination"), each = 5),
    BioCategory_Manual = rep(c("Neuro","Neuro","Immune","Glia"), each = 5),
    Gene = sample(paste0("G", 1:10), size = 20, replace = TRUE)
  )
  p <- plot_pathway_gene_heatmap(df, export = FALSE)
  expect_s3_class(p, "ggplot")
})
