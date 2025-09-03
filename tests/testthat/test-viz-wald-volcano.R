test_that("plot_wald_volcano returns a ggplot", {
  df <- tibble::tibble(
    Protein = c("A","B","C","D"),
    `Beta 1` = c(0.2, 0.1, -0.1, 0.05),
    `Beta 3` = c(-0.2, 0.3, -0.05, 0.20),
    `P-value 1` = c(1e-5, 0.02, 0.5, 0.001),
    `Adjusted P-value 1` = c(2e-5, 0.04, 0.6, 0.005)
  )
  p <- plot_wald_volcano(
    wald_df = df,
    pval_col = "P-value 1",
    fdr_col = "Adjusted P-value 1",  # allow our local column name
    beta_before_col = "Beta 1",
    beta_after_col  = "Beta 3",
    protein_col = "Protein",
    export = FALSE
  )
  expect_s3_class(p, "ggplot")
})
