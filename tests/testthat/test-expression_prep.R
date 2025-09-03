test_that("prepare_combined_expression pivots, coerces, and labels correctly", {
  genes <- c("G1","G2")
  df_all <- tibble::tibble(
    CATEGORY = c("Normal","SCD","MCI","AD Dementia"),
    G1 = c("1.0","2.0","3.0","4.0"),
    G2 = c("5.0","6.0","7.0","8.0")
  )
  df_normal_only   <- df_all[1, , drop = FALSE]
  df_status_change <- df_all[2, , drop = FALSE]
  df_abnormal_only <- df_all[4, , drop = FALSE]

  out <- prepare_combined_expression(
    df_normal_only = df_normal_only,
    df_status_change = df_status_change,
    df_abnormal_only = df_abnormal_only,
    df_all = df_all,
    subset_genes = genes,
    category_col = "CATEGORY",
    categories = c("SCD","MCI"),
    normal_label = "Normal_only",
    status_label = "Status_change",
    abnormal_label = "Abnormal_only"
  )

  expect_true(all(c("Gene","Expression","Source") %in% names(out)))
  expect_s3_class(out$Gene, "factor")
  expect_true(is.ordered(out$Gene))
  expect_equal(levels(out$Gene), genes)
  expect_type(out$Expression, "double")

  # Sources should include cohort labels + categories with underscores
  expect_true(all(c("Normal_only","Status_change","Abnormal_only","SCD","MCI") %in% unique(out$Source)))
})
