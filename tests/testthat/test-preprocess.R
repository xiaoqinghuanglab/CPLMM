test_that("time transforms work", {
  df <- tibble::tibble(age = c(60, 62), onset_age = c(58, 60))
  out1 <- calculate_years_since_onset(df)
  expect_true("years_since_onset" %in% names(out1))
  expect_equal(out1$years_since_onset, c(2, 2))

  out2 <- add_piecewise_age(df)
  expect_true("piecewise_age" %in% names(out2))
  expect_equal(out2$piecewise_age, c(2, 2))
})

test_that("set_onset_age_for_normals updates only normal-only subjects", {
  df <- tibble::tibble(
    subject_id = c(1,1,2,2,3,3),
    status_raw = c("Normal","Normal","Normal","Abnormal","Normal","Normal"),
    age = c(60, 62, 61, 63, 55, 58)
  )
  out <- set_onset_age_for_normals(df)
  # subject 1 & 3 are normal-only; subject 2 is mixed
  expect_true(all(out$decage_mutated[out$subject_id == 1] == 62))
  expect_true(all(out$decage_mutated[out$subject_id == 3] == 58))
  expect_true(all(is.na(out$decage_mutated[out$subject_id == 2])))
})

test_that("correct_status_by_onset flips statuses according to age vs onset", {
  df <- tibble::tibble(
    age = c(60, 65),
    onset_age = c(62, 62),
    status_cleaned = c("Abnormal", "Normal")
  )
  out <- correct_status_by_onset(df)
  expect_equal(out$status_cleaned, c("Normal","Abnormal"))
})

test_that("enforce_unidirectional_status_change sets all later to Abnormal", {
  df <- tibble::tibble(
    subject_id = c(1,1,1,1),
    procedure_date = as.Date("2020-01-01") + c(0, 10, 20, 30),
    status_cleaned = c("Normal","Normal","Abnormal","Normal")
  )
  out <- enforce_unidirectional_status_change(df)
  expect_equal(out$status_cleaned, c("Normal","Normal","Abnormal","Abnormal"))
})

test_that("identify_status_change_subjects keeps only changers", {
  df <- tibble::tibble(
    subject_id = c(1,1,2,2,3,3),
    procedure_date = as.Date("2020-01-01") + c(0,5,0,5,0,5),
    status_cleaned = c("Normal","Abnormal", "Normal","Normal", "Abnormal","Abnormal")
  )
  keep <- identify_status_change_subjects(df)
  expect_true(all(unique(keep$subject_id) == 1))
})

test_that("get_status_groups splits normal-only and abnormal-only correctly", {
  df <- tibble::tibble(
    subject_id = c(1,1,2,2,3,3),
    status_cleaned = c("Normal","Normal", "Abnormal","Abnormal", "Normal","Abnormal")
  )
  parts <- get_status_groups(df)
  expect_true(setequal(unique(parts$df_normal$subject_id), 1))
  expect_true(setequal(unique(parts$df_abnormal$subject_id), 2))
})
