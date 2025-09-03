test_that("compare_groups_mannwhitney computes stats and sorts", {
  set.seed(123)
  genes <- c("G1","G2","G3")
  df <- tibble::tibble(
    Gene = rep(genes, each = 20),
    Expression = c(rnorm(20, 0, 1),  # G1 group spread
                   rnorm(10, 0.0, 1), rnorm(10, 1.0, 1),  # G2: group2 shifted up
                   rnorm(20, 0, 1)),  # G3 no strong shift
    Source = rep(rep(c("A","B"), each = 10), times = 3)
  )

  res <- compare_groups_mannwhitney(
    combined_expr = df,
    gene_list = genes,
    group1 = "A",
    group2 = "B",
    alpha = 0.05,
    rank_by = "FDR"
  )

  expect_true(all(c("Gene","A_mean","B_mean","Delta_mean","U_statistic","p_value","FDR","Significant") %in% names(res)))
  expect_equal(nrow(res), 3L)
  # With G2 having a mean shift, it should usually be near the top (smallest FDR)
  expect_true(res$Gene[1] %in% c("G2"))
})
