#' Calculate years since disease onset
#'
#' @param df A data.frame/tibble.
#' @param age_col Column name for age. Default: `"age"`.
#' @param onset_age_col Column name for onset age. Default: `"onset_age"`.
#' @param new_col Name of the new column to create. Default: `"years_since_onset"`.
#'
#' @return A tibble with `new_col` added as `age - onset_age`.
#' @export
calculate_years_since_onset <- function(df,
                                        age_col = "age",
                                        onset_age_col = "onset_age",
                                        new_col = "years_since_onset") {
  df <- tibble::as_tibble(df)
  stopifnot(all(c(age_col, onset_age_col) %in% names(df)))
  df[[new_col]] <- df[[age_col]] - df[[onset_age_col]]
  df
}

#' Create piecewise age (relative to onset)
#'
#' Mirrors the Python behavior: `piecewise_age = age - onset_age`.
#'
#' @inheritParams calculate_years_since_onset
#' @param new_col Name of the new column. Default: `"piecewise_age"`.
#'
#' @return A tibble with `new_col` added.
#' @export
add_piecewise_age <- function(df,
                              age_col = "age",
                              onset_age_col = "onset_age",
                              new_col = "piecewise_age") {
  calculate_years_since_onset(df, age_col = age_col, onset_age_col = onset_age_col, new_col = new_col)
}

#' Set onset age for subjects who remain Normal across all visits
#'
#' For subjects with `status_col == "Normal"` at every visit, set `mutated_col`
#' equal to the subject-specific maximum of `age_col` (in-place column creation if missing).
#'
#' @param df A data.frame/tibble.
#' @param subject_id_col Subject ID column. Default: `"subject_id"`.
#' @param status_col Status column (raw). Default: `"status_raw"`.
#' @param age_col Age column. Default: `"age"`.
#' @param mutated_col Column to write max-age into. Default: `"decage_mutated"`.
#'
#' @return A tibble with `mutated_col` updated for normal-only subjects.
#' @export
set_onset_age_for_normals <- function(df,
                                      subject_id_col = "subject_id",
                                      status_col = "status_raw",
                                      age_col = "age",
                                      mutated_col = "decage_mutated") {
  df <- tibble::as_tibble(df)
  stopifnot(all(c(subject_id_col, status_col, age_col) %in% names(df)))

  # Ensure target column exists
  if (!mutated_col %in% names(df)) df[[mutated_col]] <- NA_real_

  by_subj <- dplyr::group_by(df, .data[[subject_id_col]])
  flags <- dplyr::summarise(by_subj, all_normal = all(.data[[status_col]] == "Normal"), .groups = "drop")
  normal_ids <- flags %>% dplyr::filter(.data$all_normal) %>% dplyr::pull(.data[[subject_id_col]])

  if (length(normal_ids)) {
    # compute max age per subject only for normal-only IDs
    max_age_map <- df %>%
      dplyr::filter(.data[[subject_id_col]] %in% normal_ids) %>%
      dplyr::group_by(.data[[subject_id_col]]) %>%
      dplyr::summarise(max_age = max(.data[[age_col]], na.rm = TRUE), .groups = "drop")
    df <- df %>%
      dplyr::left_join(max_age_map, by = rlang::set_names(subject_id_col, subject_id_col)) %>%
      dplyr::mutate(
        !!mutated_col := dplyr::if_else(.data[[subject_id_col]] %in% normal_ids,
                                        .data$max_age, .data[[mutated_col]])
      ) %>%
      dplyr::mutate(max_age = NULL)
  }

  df
}

#' Correct status labels using onset age logic
#'
#' If `onset_age > age` and status is `"Abnormal"`, set to `"Normal"`.
#' If `onset_age < age` and status is `"Normal"`, set to `"Abnormal"`.
#'
#' @param df A data.frame/tibble.
#' @param age_col Age column. Default: `"age"`.
#' @param onset_age_col Onset age column. Default: `"onset_age"`.
#' @param status_col Cleaned status column to modify. Default: `"status_cleaned"`.
#'
#' @return A tibble with corrected `status_col`.
#' @export
correct_status_by_onset <- function(df,
                                    age_col = "age",
                                    onset_age_col = "onset_age",
                                    status_col = "status_cleaned") {
  df <- tibble::as_tibble(df)
  stopifnot(all(c(age_col, onset_age_col, status_col) %in% names(df)))

  too_early <- (df[[onset_age_col]] > df[[age_col]]) & (df[[status_col]] == "Abnormal")
  too_late  <- (df[[onset_age_col]] < df[[age_col]]) & (df[[status_col]] == "Normal")

  df[[status_col]][too_early] <- "Normal"
  df[[status_col]][too_late]  <- "Abnormal"
  df
}

#' Enforce unidirectional status change: once Abnormal, always Abnormal
#'
#' Sorts by subject and date, and sets all visits *at or after* the first
#' `"Abnormal"` status to `"Abnormal"`.
#'
#' @param df A data.frame/tibble.
#' @param subject_id_col Subject ID column. Default: `"subject_id"`.
#' @param status_col Status column to modify. Default: `"status_cleaned"`.
#' @param date_col Date/time column used for ordering. Default: `"procedure_date"`.
#'
#' @return A tibble with monotonized `status_col`.
#' @export
enforce_unidirectional_status_change <- function(df,
                                                 subject_id_col = "subject_id",
                                                 status_col = "status_cleaned",
                                                 date_col = "procedure_date") {
  df <- tibble::as_tibble(df)
  stopifnot(all(c(subject_id_col, status_col, date_col) %in% names(df)))

  df %>%
    dplyr::arrange(.data[[subject_id_col]], .data[[date_col]]) %>%
    dplyr::group_by(.data[[subject_id_col]]) %>%
    dplyr::mutate(
      .abn_seen = dplyr::cumany(.data[[status_col]] == "Abnormal"),
      !!status_col := dplyr::if_else(.abn_seen, "Abnormal", .data[[status_col]])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(.abn_seen = NULL)
}

#' Keep only subjects who change from Normal to Abnormal
#'
#' Returns subjects who have at least one `"Normal"` visit and at least one
#' non-Normal visit, preserving order by subject and date.
#'
#' @param df A data.frame/tibble.
#' @param subject_id_col Subject ID column. Default: `"subject_id"`.
#' @param status_col Status column (cleaned). Default: `"status_cleaned"`.
#' @param date_col Date/time column. Default: `"procedure_date"`.
#'
#' @return A tibble containing only status-change subjects.
#' @export
identify_status_change_subjects <- function(df,
                                            subject_id_col = "subject_id",
                                            status_col = "status_cleaned",
                                            date_col = "procedure_date") {
  df <- tibble::as_tibble(df)
  stopifnot(all(c(subject_id_col, status_col, date_col) %in% names(df)))

  df %>%
    dplyr::arrange(.data[[subject_id_col]], .data[[date_col]]) %>%
    dplyr::group_by(.data[[subject_id_col]]) %>%
    dplyr::filter(any(.data[[status_col]] == "Normal") & any(.data[[status_col]] != "Normal")) %>%
    dplyr::ungroup()
}

#' Split into normal-only and abnormal-only subject groups
#'
#' Normal-only: all visits have `"Normal"`.
#' Abnormal-only: no visits have `"Normal"`.
#'
#' @param df A data.frame/tibble.
#' @param subject_id_col Subject ID column. Default: `"subject_id"`.
#' @param status_col Status column (cleaned). Default: `"status_cleaned"`.
#'
#' @return A named list with elements `df_normal` and `df_abnormal`.
#' @export
get_status_groups <- function(df,
                              subject_id_col = "subject_id",
                              status_col = "status_cleaned") {
  df <- tibble::as_tibble(df)
  stopifnot(all(c(subject_id_col, status_col) %in% names(df)))

  flags <- df %>%
    dplyr::group_by(.data[[subject_id_col]]) %>%
    dplyr::summarise(
      all_normal = all(.data[[status_col]] == "Normal"),
      none_normal = all(.data[[status_col]] != "Normal"),
      .groups = "drop"
    )

  normal_ids   <- flags %>% dplyr::filter(.data$all_normal)  %>% dplyr::pull(.data[[subject_id_col]])
  abnormal_ids <- flags %>% dplyr::filter(.data$none_normal) %>% dplyr::pull(.data[[subject_id_col]])

  list(
    df_normal   = df %>% dplyr::filter(.data[[subject_id_col]] %in% normal_ids),
    df_abnormal = df %>% dplyr::filter(.data[[subject_id_col]] %in% abnormal_ids)
  )
}
