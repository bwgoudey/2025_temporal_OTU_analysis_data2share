
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
source("./R/data_functions.R")

# ---------------
library(doFuture)
library(future)

n_cores=parallel::detectCores() - 1
registerDoFuture()
plan(multisession, workers = n_cores)

message(glue("Parallel processing enabled with {n_cores} cores."))
#-----------


# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
cfg_debug = 0
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
cfg_holdout_samples     <- 10  # final samples held out
cfg_feature_sets_to_use <- c("fs_P", "fs_S", "fs_Mcomm", "fs_SPM", "fs_SP")
if(cfg_debug)
  cfg_feature_sets_to_use = c("fs_P")
cfg_debug_nsplits       <- ifelse(cfg_debug, 5, Inf)   # Inf for all

#model tuning
cfg_models_subset <- character()
if(cfg_debug)
  cfg_models_subset <- c(str_c(cfg_feature_sets_to_use, "_lasso"), "arima_arima", "lastval_lastobs")
cfg_pls_levels = ifelse(cfg_debug, 3, 5)
cfg_lasso_levels = ifelse(cfg_debug, 3, 7)


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
microbial_diversity_cols <- c("uniqueOTUs", "Shannon", "Simpson", "invSimpson", "PielouEvenness", "Soerensen", "SoerensenAD1AD2",
                              "BrayCurtisAD1AD2", "JaccardAD1AD2", "Rolling_Jaccard6days", "Rolling_Bray6days",
                              "Rate_of_change_Jaccard_per_day", "Rate_of_change_Bray_per_day")

cfg_feature_sets <- list(
  fs_S   = seasonal_cols,
  fs_Mcomm   = c( microbial_comm_cols),
  fs_Mdiv   = c( microbial_diversity_cols),
  fs_M   = c(microbial_otu_cols, microbial_comm_cols, microbial_diversity_cols),
  fs_SP  = union(seasonal_cols, process_cols),
  fs_SPM = union(seasonal_cols, union(process_cols, c(microbial_comm_cols, microbial_diversity_cols))),
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
make_model_specs <- function() {
  list(
    lm = parsnip::linear_reg() %>% set_engine("lm"),
    lasso = parsnip::linear_reg(penalty = tune(), mixture = 1) %>% set_engine("glmnet"),
    pls = parsnip::pls(mode = "regression", num_comp = tune(), predictor_prop = tune()) %>% set_engine("mixOmics"),
    arima = modeltime::arima_reg() %>% set_engine("auto_arima"),
    lastobs = parsnip::linear_reg() %>% set_engine("lm")
  )
}

lasso_grid <- grid_regular(penalty(range = c(-3, 0), trans = scales::log10_trans()), levels = cfg_lasso_levels)
pls_grid <- grid_space_filling(num_comp(range = c(1L, 6L)), predictor_prop(range = c(0.2, 1)), size = cfg_pls_levels)
model_specs <- make_model_specs()

# ------------------------------------------------------------------------------
# 1. Define Base Recipes (OUTSIDE LOOP)
# ------------------------------------------------------------------------------
# We use 'df_master' to establish column names and types.
# We do NOT set an outcome yet.
base_recipes <- sapply(feature_sets, function(fs_cols) {

  recipes::recipe(df_master) %>%
    # 1. Assign roles manually (No formula used)
    recipes::update_role(!!sym(cfg_time_col), new_role = "time_index") %>%
    recipes::update_role(!!sym(cfg_id_col), new_role = "id") %>%
    recipes::update_role(all_of(fs_cols), new_role = "predictor") %>%

    # 2. Feature Engineering (Identical for all horizons)
    timetk::step_timeseries_signature(!!sym(cfg_time_col)) %>%
    timetk::step_fourier(!!sym(cfg_time_col), period = 7, K = 1) %>%
    recipes::step_lag(all_of(fs_cols), lag = c(1, 3, 7, 14)) %>%

    # Rolling features
    timetk::step_slidify_augment(
      all_of(fs_cols), period = cfg_summary_window, .f = ~mean(.x, na.rm=TRUE),
      align = "right", partial = TRUE, prefix = glue("roll{cfg_summary_window}_mean_")
    ) %>%
    timetk::step_slidify_augment(
      all_of(fs_cols), period = cfg_summary_window, .f = ~sd(.x, na.rm=TRUE),
      align = "right", partial = TRUE, prefix = glue("roll{cfg_summary_window}_sd_")
    ) %>%
    timetk::step_slidify_augment(
      all_of(fs_cols), period = cfg_summary_window,
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

# ------------------------------------------------------------------------------
# Baseline Recipe Helpers
# ------------------------------------------------------------------------------
make_arima_recipe <- function(df_train, target, time_col, id_col) {
  recipes::recipe(stats::as.formula(sprintf("%s ~ %s + %s", target, time_col, id_col)), data = df_train) %>%
    recipes::update_role(!!rlang::sym(target), new_role = "outcome") %>%
    recipes::update_role(!!rlang::sym(id_col), new_role = "id") %>%
    recipes::add_role(!!rlang::sym(time_col), new_role = "time_index") %>%
    recipes::step_naomit(recipes::all_outcomes())
}

make_lastval_recipe <- function(df_train, target, time_col, id_col) {
  recipes::recipe(stats::as.formula(sprintf("%s ~ target_last_value + %s + %s", target, time_col, id_col)), data = df_train) %>%
    recipes::update_role(!!rlang::sym(target), new_role = "outcome") %>%
    recipes::update_role(target_last_value, new_role = "predictor") %>%
    recipes::update_role(!!rlang::sym(time_col), new_role = "time_index") %>%
    recipes::update_role(!!rlang::sym(id_col), new_role = "id") %>%
    recipes::step_naomit(recipes::all_outcomes())
}

# ------------------------------------------------------------------------------
# Main Forecasting Loop
# ------------------------------------------------------------------------------
all_metrics_list <- list()
all_holdout_preds <- list()
all_assessment_preds <- list()
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
    mutate(.row_id = row_number()) %>%
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

  # If we've turned on the debugging flag, then only use a subset of the time data for
  # training and evaluation purposes.
  # NB: I'm not actually sure this saves all that much time but I'm doing it anyway
  if(cfg_debug_nsplits < nrow(resamples_sliding)) {
    resamples_class <- class(resamples_sliding)
    resamples_sliding <- resamples_sliding[1:cfg_debug_nsplits, ]
    class(resamples_sliding) <- resamples_class

  }

  # --- 2. Recipe Finalization ---

  # A. The ML Recipes (from Base Recipes)
  recs_ml <- base_recipes %>%
    map(~ .x %>%
          recipes::update_role(!!sym(curr_target), new_role = "outcome") %>%
          recipes::step_naomit(all_outcomes())
    )

  rec_arima <- make_arima_recipe(df_train, curr_target, cfg_time_col, cfg_id_col)
  rec_last <- make_lastval_recipe(df_train, curr_target, cfg_time_col, cfg_id_col)
  recipes_for_workflows <- c(recs_ml, list(lastval = rec_last, arima = rec_arima))

  # --- 3. Build Unified Workflow Set -----------
  all_workflows <- workflow_set(
    preproc = recipes_for_workflows,
    models = model_specs,
    cross = TRUE
  ) %>%
    filter(
      !(grepl("^arima_", wflow_id)   & !grepl("_arima$", wflow_id)),
      !(grepl("_arima$", wflow_id)   & !grepl("^arima_", wflow_id)),
      !(grepl("^lastval_", wflow_id) & !grepl("_lastobs$", wflow_id)),
      !(grepl("_lastobs$", wflow_id) & !grepl("^lastval_", wflow_id))
    ) %>%
    add_tune_grid_to_workflows("lasso", lasso_grid) %>%
    add_tune_grid_to_workflows("pls", pls_grid)

  # --- 4. Execute the workfllows -----------
  # workflow_map handles standard resampling (LM/ARIMA) AND tuning (Lasso/PLS) simultaneously

  #If debugging, we may only want to run a subset of models
  if(exists("cfg_models_subset") && length(cfg_models_subset))
    all_workflows <- all_workflows %>% filter(wflow_id %in% cfg_models_subset)

  message("   -> Fitting all models...")
  results_master <- all_workflows %>%
    workflow_map(
      fn = "tune_grid",
      resamples = resamples_sliding,
      metrics = metrics_set,
      control = ctrl_resamples,
      verbose = TRUE
    )

  # --- 5. Extract Results  ---

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

  # C. Assessment Predictions (per resample)
  horizon_assessment_preds <- results_master %>%
    mutate(preds = map(result, ~collect_predictions(.x, summarize = FALSE))) %>%
    select(model = wflow_id, preds) %>%
    unnest(preds) %>%
    mutate(horizon = h, data_role = "assessment") %>%
    rename(.row_id = .row) %>%
    left_join(
      df_train %>% select(.row_id, !!sym(cfg_time_col), !!sym(cfg_id_col), .y_target),
      by = ".row_id"
    ) %>%
    select(-.row_id)

  all_assessment_preds[[as.character(h)]] <- horizon_assessment_preds

  # D. Holdout Predictions
  message("   -> Predicting holdout...")

  horizon_holdout <- results_master %>%
    mutate(preds = pmap(list(wflow_id, result, info), function(id, res, info) {
      fit_and_predict_holdout(id, res, info$workflow[[1]], df_train, df_test)
    })) %>%
    select(preds) %>%
    unnest(preds) %>%
    mutate(horizon = h, data_role = "holdout")

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
final_holdout_preds   <- bind_rows(all_holdout_preds)
final_assessment_preds <- bind_rows(all_assessment_preds)
final_logs    <- bind_rows(all_tuning_logs)

write_csv(final_metrics, "outputs/tables/summary_metrics_all_horizons.csv")
write_csv(final_holdout_preds, "outputs/tables/holdout_predictions.csv")
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
  filter(!grepl("_lm", model)) %>%
  filter(.metric == "rmse") %>%
  ggplot(aes(x = model, y = .estimate, colour = model)) +
  geom_jitter(width = 0.2, alpha = 0.4, show.legend = FALSE) +
  geom_boxplot(aes(fill=model)) +
  facet_grid(data_role~horizon, scales = "free_y") + theme_bw() +
  labs(title = "RMSE by Horizon", y = "RMSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Assessment + Holdout Predictions (timeline)
if(nrow(final_holdout_preds) + nrow(final_assessment_preds) > 0) {
  actual_series <- df_master %>%
    arrange(.data[[cfg_time_col]]) %>%
    select(all_of(cfg_time_col), actual = all_of(cfg_outcome_col)) %>%
    distinct()

  prediction_plot_data <- bind_rows(final_assessment_preds, final_holdout_preds) %>%
    mutate(Sample_Date1 = Sample_Date+horizon,
           horizon = factor(horizon))

  prediction_plot_data %>%
    ggplot() +
    geom_line(data = actual_series, aes(x = Sample_Date, y = actual), colour = "black", linewidth = 0.8) +
    geom_point(data = actual_series, aes(x = Sample_Date, y = actual), colour = "black", alpha = 0.4, size = 1) +
    geom_point(aes(x = Sample_Date1, y = .pred, colour = model, shape = horizon, alpha = data_role), size = 2) +
    theme_minimal(base_size = 13) +
    scale_alpha_manual(values = c(assessment = 0.5, holdout = 1), name = "Data Role") +
    labs(title = "Assessment & Holdout Predictions", y = "Outcome", x = cfg_time_col, colour = "Model", shape = "Horizon")
}
