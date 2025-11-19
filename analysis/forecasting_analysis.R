
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(plsmod)
  library(slider)
  library(tidymodels)
  library(rsample)
  library(glue)
  library(digest)
  library(modeltime)
  library(timetk)
})

tidymodels_prefer()
source("R/data_functions.R")

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
cfg_data_path     <- "./data/mdata.species.clr_all.csv"
cfg_id_col        <- "Digester"
cfg_time_col      <- "Sample_Date"
cfg_outcome_col   <- "A"

cfg_primary_digesters <- c("D1", "D2")

# Forecast settings
cfg_forecast_horizons   <- c(1, 7, 14, 21)
cfg_assess_width_days   <- 14
cfg_lookback_days       <- 56
cfg_summary_window      <- 6   # feature rolling window
cfg_holdout_samples     <- 8L  # final samples held out
cfg_feature_sets_to_use <- c("fs_S", "fs_M")
cfg_debug_nsplits       <- 5   # Inf for all

# ------------------------------------------------------------------------------
# Data Loading & Feature Definitions
# ------------------------------------------------------------------------------
df_raw <- read_csv(cfg_data_path, show_col_types = FALSE) %>%
  mutate(!!cfg_time_col := as_date(.data[[cfg_time_col]]))

# Define column groups
process_cols <- c("T1", "A", "B", "V2", "P", "V3", "T2", "C", "H", "C3", "V", "AA", "O", "TP", "T3", "G", "H2", "C2")
seasonal_cols <- c("BOM_min_temp", "BOM_max_temp", "BOM_3pm_temp", "BOM_3pm_pressure", "S_sin_doy", "S_cos_doy")
microbial_otu_cols <- grep("_otu$", names(df_raw), value = TRUE)
microbial_comm_cols <- c("su__sugar_degraders_fermenters", "aa__amino_acid_degraders", "c4__butyrate_and_valerate_degraders",
                         "fa__long_chain_fatty_acid_degraders", "pro__propionate_degraders", "ac__acetoclastic_methanogens",
                         "h2__hydrogenotrophic_methanogens")
diversity_cols <- c("uniqueOTUs", "Shannon", "Simpson", "invSimpson", "PielouEvenness", "Soerensen", "SoerensenAD1AD2",
                    "BrayCurtisAD1AD2", "JaccardAD1AD2", "Rolling_Jaccard6days", "Rolling_Bray6days",
                    "Rate_of_change_Jaccard_per_day", "Rate_of_change_Bray_per_day")

cfg_feature_sets <- list(
  fs_S   = seasonal_cols,
  fs_M   = c(microbial_otu_cols, microbial_comm_cols, diversity_cols),
  fs_SP  = union(seasonal_cols, process_cols),
  fs_SPM = union(seasonal_cols, union(process_cols, c(microbial_otu_cols, microbial_comm_cols, diversity_cols))),
  fs_P   = process_cols
)
feature_sets <- cfg_feature_sets[cfg_feature_sets_to_use]

# ------------------------------------------------------------------------------
# Output Setup
# ------------------------------------------------------------------------------
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

metrics_set <- yardstick::metric_set(yardstick::rmse, yardstick::mae, yardstick::rsq_trad)
ctrl_resamples <- tune::control_resamples(save_pred = TRUE, verbose = FALSE)

# ------------------------------------------------------------------------------
# Pre-Loop Data Preparation
# ------------------------------------------------------------------------------
# 1. Aggregation and Seasonality
df_primary <- df_raw %>%
  filter(.data[[cfg_id_col]] %in% cfg_primary_digesters) %>%
  add_seasonal_terms(cfg_time_col)

df_agg <- aggregate_digesters(df_primary, cfg_id_col, cfg_time_col) %>%
  arrange(.data[[cfg_id_col]], .data[[cfg_time_col]]) %>%
  group_by(.data[[cfg_id_col]]) %>%
  mutate(target_last_value = .data[[cfg_outcome_col]]) %>% # Current value (lag0)
  ungroup()

# 2. Generate ALL targets upfront (Moving target generation out of the loop)
# We create columns like 'A_lead1', 'A_lead7', 'A_lead14' etc.
df_master <- df_agg
target_map <- list()

for(h in cfg_forecast_horizons) {
  col_name <- glue("{cfg_outcome_col}_lead{h}")
  target_map[[as.character(h)]] <- col_name

  df_master <- df_master %>%
    group_by(.data[[cfg_id_col]]) %>%
    mutate(!!col_name := lead(.data[[cfg_outcome_col]], n = h)) %>%
    ungroup()
}

# ------------------------------------------------------------------------------
# Model Specifications
# ------------------------------------------------------------------------------
# Shared specs
lm_spec <- parsnip::linear_reg() %>% set_engine("lm")
lasso_spec <- parsnip::linear_reg(penalty = tune(), mixture = 1) %>% set_engine("glmnet")
lasso_grid <- grid_regular(penalty(range = c(-3, 0), trans = scales::log10_trans()), levels = 7)

pls_spec <- parsnip::pls(mode = "regression", num_comp = tune(), predictor_prop = tune()) %>% set_engine("mixOmics")
pls_grid <- grid_space_filling(num_comp(range = c(1L, 6L)), predictor_prop(range = c(0.2, 1)), size = 10)

arima_spec <- modeltime::arima_reg() %>% set_engine("auto_arima")

# ------------------------------------------------------------------------------
# 1. Define Base Recipes (OUTSIDE LOOP)
# ------------------------------------------------------------------------------
# We use 'df_master' to establish column names and types.
# We do NOT set an outcome yet.
base_recipes <- sapply(names(feature_sets), function(x) {

  recipes::recipe(df_master) %>%
    # 1. Assign roles manually (No formula used)
    recipes::update_role(!!sym(cfg_time_col), new_role = "time_index") %>%
    recipes::update_role(!!sym(cfg_id_col), new_role = "id") %>%
    recipes::update_role(all_of(feature_sets[[x]]), new_role = "predictor") %>%

    # 2. Feature Engineering (Identical for all horizons)
    timetk::step_timeseries_signature(!!sym(cfg_time_col)) %>%
    timetk::step_fourier(!!sym(cfg_time_col), period = 7, K = 1) %>%
    recipes::step_lag(all_of(feature_sets[[x]]), lag = c(1, 3, 7, 14)) %>%

    # Rolling features
    timetk::step_slidify_augment(
      all_of(feature_sets[[x]]), period = cfg_summary_window, .f = ~mean(.x, na.rm=TRUE),
      align = "right", partial = TRUE, prefix = glue("roll{cfg_summary_window}_mean_")
    ) %>%
    timetk::step_slidify_augment(
      all_of(feature_sets[[x]]), period = cfg_summary_window, .f = ~sd(.x, na.rm=TRUE),
      align = "right", partial = TRUE, prefix = glue("roll{cfg_summary_window}_sd_")
    ) %>%
    timetk::step_slidify_augment(
      all_of(feature_sets[[x]]), period = cfg_summary_window,
      .f = ~{val<-suppressWarnings(max(.x, na.rm=TRUE)); if(is.finite(val)) val else NA_real_},
      align = "right", partial = TRUE, prefix = glue("roll{cfg_summary_window}_max_")
    ) %>%

    # 3. Cleanup and Normalization
    recipes::step_rm(ends_with(".lbl")) %>%
    recipes::step_zv(all_predictors()) %>%
    recipes::step_impute_mean(all_numeric_predictors()) %>%
    recipes::step_impute_mode(all_nominal_predictors()) %>%
    recipes::step_corr(all_predictors(), threshold = 0.9) %>%
    recipes::step_normalize(all_predictors())

}, simplify = FALSE, USE.NAMES = TRUE)

# extract_best_and_predict:
# Takes a row from a workflow_set, finalizes it, fits on Train, predicts Test
fit_and_predict_holdout <- function(wflow_id, result, workflow, train_data, test_data) {

  if(nrow(test_data) == 0) return(tibble())

  # 1. Determine Best Configuration (if tunable) or just use the workflow (if not)
  # Note: 'result' is the tuning result object
  best_config <- tryCatch(select_best(result, metric = "rmse"), error = function(e) NULL)

  final_wf <- if (!is.null(best_config)) {
    finalize_workflow(workflow, best_config)
  } else {
    workflow
  }

  # 2. Fit on ALL training data
  fit_final <- parsnip::fit(final_wf, data = train_data)

  # 3. Predict on Holdout
  preds <- predict(fit_final, new_data = test_data)

  # 4. Bind with ID/Date/Truth
  bind_cols(
    test_data %>% select(any_of(c("Digester", "Sample_Date", ".y_target"))),
    preds
  ) %>%
    mutate(model = wflow_id)
}

# ------------------------------------------------------------------------------
# Main Forecasting Loop
# ------------------------------------------------------------------------------
all_metrics_list <- list()
all_holdout_preds <- list()
all_tuning_logs <- list()

message(glue("Starting forecast loop. Horizons: {paste(cfg_forecast_horizons, collapse=', ')}"))

for (h in cfg_forecast_horizons) {

  # --- 1. Setup Data & Target ---
  curr_target <- target_map[[as.character(h)]]
  message(glue("--- Horizon: {h} | Target: {curr_target} ---"))

  df_h <- df_master %>%
    filter(!is.na(.data[[curr_target]])) %>%
    mutate(.y_target = .data[[curr_target]])

  # Split
  df_h <- df_h %>%
    arrange(.data[[cfg_id_col]], .data[[cfg_time_col]]) %>%
    group_by(.data[[cfg_id_col]]) %>%
    mutate(is_holdout = row_number() > pmax(n() - cfg_holdout_samples, 0L)) %>%
    ungroup()

  df_train <- df_h %>% filter(!is_holdout) %>% select(-is_holdout)
  df_test  <- df_h %>% filter(is_holdout)  %>% select(-is_holdout)

  if(nrow(df_train) == 0) stop("Empty training set")

  # Resamples
  resamples_sliding <- sliding_period(
    df_train, index = !!sym(cfg_time_col), period = "day",
    lookback = cfg_lookback_days, assess_start = 1, assess_stop = cfg_assess_width_days,
    step = 1, complete = TRUE
  )

  if(cfg_debug_nsplits < nrow(resamples_sliding)) {
    resamples_sliding <- resamples_sliding[1:cfg_debug_nsplits, ]
  }

  # --- 2. Recipe Finalization ---

  # A. The ML Recipes (from Base Recipes)
  recs_ml <- base_recipes %>%
    map(~ .x %>%
          recipes::update_role(!!sym(curr_target), new_role = "outcome") %>%
          recipes::step_naomit(all_outcomes())
    )

  # B. The Baseline Recipes (Created specifically for this target)
  # ARIMA: Needs Date + Outcome only
  rec_arima <- recipe(df_train) %>%
    update_role(!!sym(cfg_time_col), new_role = "time_index") %>%
    update_role(!!sym(curr_target), new_role = "outcome") %>%
    step_naomit(all_outcomes())

  # Last Value: Needs LastVal + Outcome + ID/Date
  rec_last <- recipe(df_train) %>%
    update_role(!!sym(curr_target), new_role = "outcome") %>%
    update_role(target_last_value, new_role = "predictor") %>%
    update_role(!!sym(cfg_time_col), new_role = "time_index") %>%
    update_role(!!sym(cfg_id_col), new_role = "id") %>%
    step_impute_mean(all_predictors())

  # --- 3. Build Unified Workflow Set ---

  # We construct the set in blocks to ensure correct pairings
  # Block 1: ML Models (All ML Recipes x LM, Lasso, PLS)
  wset_ml <- workflow_set(
    preproc = recs_ml,
    models  = list(lm = lm_spec, lasso = lasso_spec, pls = pls_spec),
    cross   = TRUE
  )

  # Block 2: ARIMA (ARIMA Recipe x ARIMA Spec)
  wset_arima <- workflow_set(
    preproc = list(simple = rec_arima),
    models  = list(arima = arima_spec),
    cross   = TRUE
  )

  # Block 3: Last Value (LastVal Recipe x LM Spec)
  wset_last <- workflow_set(
    preproc = list(baseline = rec_last),
    models  = list(lastobs = lm_spec),
    cross   = TRUE
  )

  # Combine and configure grids
  all_workflows <- bind_rows(wset_ml, wset_arima, wset_last) %>%
    # Attach specific grids to specific model types using Regex on wflow_id
    option_add(grid = lasso_grid, id = "lasso") %>%
    option_add(grid = pls_grid,   id = "pls")

  # --- 4. Execute (The "One Ring" Rule) ---
  # workflow_map handles standard resampling (LM/ARIMA) AND tuning (Lasso/PLS) simultaneously

  message("   -> Fitting all models...")
  results_master <- all_workflows %>%
    workflow_map(
      fn = "tune_grid",
      resamples = resamples_sliding,
      metrics = metrics_set,
      control = ctrl_resamples,
      verbose = TRUE
    )

  # --- 5. Extract Results (Operating on the Master Object) ---

  # A. Tuning Logs (Detailed candidates)
  horizon_log <- results_master %>%
    mutate(metrics = map(result, ~collect_metrics(.x, summarize = FALSE))) %>%
    select(wflow_id, metrics) %>%
    unnest(metrics) %>%
    mutate(horizon = h)

  all_tuning_logs[[as.character(h)]] <- horizon_log

  # B. Assessment Metrics (Best candidate per model)
  horizon_assess <- results_master %>%
    mutate(best_metrics = map(result, function(res) {
      # select_best works for tuned models; for non-tuned, it effectively selects the only result
      tryCatch({
        best <- select_best(res, metric = "rmse")
        collect_metrics(res, summarize = FALSE) %>% semi_join(best, by = ".config")
      }, error = function(e) collect_metrics(res, summarize = FALSE)) # Fallback for simple models
    })) %>%
    select(model = wflow_id, best_metrics) %>%
    unnest(best_metrics) %>%
    mutate(horizon = h, data_role = "assessment")

  # C. Holdout Predictions
  message("   -> Predicting holdout...")

  horizon_holdout <- results_master %>%
    mutate(preds = pmap(list(wflow_id, result, info), function(id, res, info) {
      fit_and_predict_holdout(id, res, info$workflow[[1]], df_train, df_test)
    })) %>%
    select(preds) %>%
    unnest(preds) %>%
    mutate(horizon = h)

  # Calculate metrics for holdout
  if(nrow(horizon_holdout) > 0) {
    horizon_holdout_metrics <- horizon_holdout %>%
      group_by(model) %>%
      metrics_set(truth = .y_target, estimate = .pred) %>%
      mutate(horizon = h, data_role = "holdout", .config = "holdout", id = "Holdout")

    # Append everything
    all_metrics_list[[as.character(h)]] <- bind_rows(horizon_assess, horizon_holdout_metrics)
    all_holdout_preds[[as.character(h)]] <- horizon_holdout
  } else {
    all_metrics_list[[as.character(h)]] <- horizon_assess
  }
}

# ------------------------------------------------------------------------------
# Aggregation and Saving
# ------------------------------------------------------------------------------
final_metrics <- bind_rows(all_metrics_list)
final_preds   <- bind_rows(all_holdout_preds)
final_logs    <- bind_rows(all_tuning_logs)

write_csv(final_metrics, "outputs/tables/summary_metrics_all_horizons.csv")
write_csv(final_preds, "outputs/tables/holdout_predictions.csv")
# Detailed tuning results for every fold, every parameter candidate, every horizon
write_csv(final_logs, "outputs/tables/full_tuning_results_log.csv")

message("Pipeline complete. Files saved to outputs/tables/")

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------
# 1. Assessment MAE Boxplot
final_metrics %>%
  filter(.metric == "mae", data_role == "assessment") %>%
  ggplot(aes(y = .estimate, x = model, fill = model)) +
  geom_boxplot() + theme_light(base_size = 14) +
  labs(y = "MAE (Assessment Folds)", title = "Model Performance Distribution") +
  theme(axis.text.x = element_blank(), legend.position = "bottom")

# 2. Performance by Horizon
final_metrics %>%
  filter(.metric == "rmse", data_role == "assessment") %>%
  ggplot(aes(x = model, y = .estimate, colour = model)) +
  geom_jitter(width = 0.2, alpha = 0.4, show.legend = FALSE) +
  stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
  facet_wrap(~horizon, scales = "free_y") + theme_bw() +
  labs(title = "RMSE by Horizon", y = "RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Holdout Predictions
if(nrow(final_preds) > 0) {
  final_preds %>%
    ggplot() +
    geom_line(aes(x = !!sym(cfg_time_col), y = .y_target, group = interaction(horizon, !!sym(cfg_id_col))), colour = "black", linewidth = 0.8) +
    geom_line(aes(x = !!sym(cfg_time_col), y = .pred, colour = model, group = interaction(model, horizon)), alpha = 0.7) +
    facet_wrap(~horizon, scales = "free_y") + theme_minimal() +
    labs(title = "Holdout Predictions", y = "Outcome")
}
