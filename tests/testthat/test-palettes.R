test_that("jama_palette returns expected names and hex codes", {
  pal <- jama_palette()
  expect_named(pal, c("Normal_only","SCD","MCI","AD_Dementia","FTD_Dementia"))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", unname(pal))))
  expect_equal(pal[["MCI"]], "#B24745")
})
