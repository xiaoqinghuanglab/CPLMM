test_that("compute_wald_test works and ranks correctly", {
  df <- tibble::tibble(
    Protein = c("A","B","C"),
    `Beta 1` = c( 0.30,  0.10, 0.00),
    `SE Beta 1` = c(0.10, 0.10, 0.10),
    `Beta 3` = c(-0.10,  0.08, 0.00),
    `SE Beta 3` = c(0.10, 0.10, 0.20),
    `Beta 2` = c( 0.02,  0.01, 0.00),
    `SE Beta 2` = c(0.10, 0.10, 0.10),
    `Beta 4` = c(-0.01,  0.02, 0.00),
    `SE Beta 4` = c(0.10, 0.10, 0.20)
  )

  res <- compute_wald_test(df, adjust_p = TRUE, rank_by = 1, alpha = 0.05)

  expect_true(all(c("Protein","Wald Statistic 1","P-value 1","Adjusted P-value 1","Rank") %in% names(res)))
  expect_equal(nrow(res), 3L)

  # Because Protein A has the largest |Beta1 - Beta3| relative to SEs, it should usually rank #1 by test 1
  expect_equal(res$Protein[1], "A")

  # Significant flags exist and are logical
  expect_type(res$`Significant 1`, "logical")
  expect_type(res$`Significant 2`, "logical")

  # Ranking by test 2 yields a valid ordering and the same row count
  res2 <- compute_wald_test(df, adjust_p = FALSE, rank_by = 2)
  expect_equal(nrow(res2), 3L)
  expect_true(all(order(res2$`P-value 2`) == seq_len(3)))
})
