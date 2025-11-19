

# Core pipeline utilities for anaerobic digester forecasting. ---------

#' Add seasonal sine/cosine terms.
add_seasonal_terms <- function(df, time_col) {
  df %>%
    dplyr::mutate(
      .doy = lubridate::yday(.data[[time_col]]),
      S_sin_doy = sin(2 * pi * .doy / 365.25),
      S_cos_doy = cos(2 * pi * .doy / 365.25)
    ) %>%
    dplyr::select(-.doy)
}


#' Aggregate numeric columns across a set of digesters.
aggregate_digesters <- function(df,
                                id_col,
                                time_col,
                                ids = c("D1", "D2"),
                                aggregated_id = "D1D2") {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  df %>%
    dplyr::filter(.data[[id_col]] %in% ids) %>%
    dplyr::group_by(.data[[time_col]]) %>%
    dplyr::summarise(dplyr::across(all_of(num_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::mutate(!!id_col := aggregated_id) %>%
    dplyr::relocate(all_of(c(id_col, time_col)))
}


#' Create leakage-safe targets for a given outcome/horizon.
#' make_forecast_targets() returns a version of df where each row has an additional
#' column .y_target, which is the outcome observed horizon_days in the future
#' for that same entity. If use_smoothed = TRUE, that outcome is first turned
#' into a trailing moving average.
make_forecast_targets <- function(df,
                         outcome_col,
                         id_col,
                         time_col,
                         horizon_days,
                         use_smoothed = TRUE,
                         smooth_window = 6) {
  df1 <- df %>%
    dplyr::arrange(.data[[id_col]], .data[[time_col]]) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(
      .y_raw = .data[[outcome_col]],
      .y_smooth = slider::slide_dbl(.y_raw, mean, .before = smooth_window - 1, .complete = TRUE)
    ) %>%
    dplyr::ungroup()

  future_vals <- df1 %>%
    dplyr::transmute(
      !!id_col := .data[[id_col]],
      !!time_col := .data[[time_col]] - lubridate::days(horizon_days),
      .y_target = if (use_smoothed) .y_smooth else .y_raw
    )

  df1 %>%
    dplyr::left_join(future_vals, by = c(id_col, time_col)) %>%
    tidyr::drop_na(.y_target)
}
