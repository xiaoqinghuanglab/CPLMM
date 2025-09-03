test_that("plot_expression_boxplot returns a ggplot", {
  set.seed(42)
  genes <- c("G1","G2","G3")
  df <- tibble::tibble(
    Gene = rep(genes, each = 60),
    Expression = rnorm(180, 0, 1),
    Source = rep(rep(c("Normal_only","SCD","MCI"), each = 20), times = 3)
  )
  p <- plot_expression_boxplot(
    combined_expr = df,
    gene_order = genes,
    rotation = 45,
    export = FALSE
  )
  expect_s3_class(p, "ggplot")
})
