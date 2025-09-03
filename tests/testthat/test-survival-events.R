test_that("compute_event_df marks first crossing and censors otherwise", {
  df <- tibble::tibble(
    SUBID = c("S1","S1","S1","S2","S2"),
    PROCEDURE_AGE = c(60, 61, 62, 70, 71),
    ONSET_AGE = c(59, 59, 59, 68, 68),
    NFL = c(10, 15, 20, 5, 6)   # threshold at 12
  )
  # S1 crosses at age 61 (first age >= 12); S2 never crosses, censor at last age 71
  ev <- compute_event_df(df, threshold = 12, biomarker = "NFL", group_label = "Status_change")

  expect_equal(sort(ev$SUBID), c("S1","S2"))
  expect_equal(ev$event[ev$SUBID == "S1"], 1L)
  expect_equal(ev$event[ev$SUBID == "S2"], 0L)

  # time is relative to ONSET_AGE (first row per subject); S1: 61-59 = 2, S2 censor: 71-68 = 3
  expect_equal(ev$time[ev$SUBID == "S1"], 2)
  expect_equal(ev$time[ev$SUBID == "S2"], 3)
  expect_true(all(ev$group == "Status_change"))
})
