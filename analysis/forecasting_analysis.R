
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(slider)
  library(tidymodels)
  library(rsample)
  library(glue)
  library(digest)
  library(modeltime)
})

tidymodels_prefer()


source("R/data_functions.R")

DATA_PATH <- "./data/mdata.species.clr_all.csv"
df_raw <- readr::read_csv(DATA_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    Sample_Date := lubridate::as_date(.data[["Sample_Date"]])
  )

# Set up the Features ------------------
nm <- names(df_raw)
process_cols=c("T1", "A", "B", "V2", "P", "V3", "T2", "C", "H", "C3", "V", "AA", "O", "TP", "T3", "G", "H2", "C2")
seasonal_cols <- c("BOM_min_temp" , "BOM_max_temp"    , "BOM_3pm_temp" ,  "BOM_3pm_pressure" , "S_sin_doy", "S_cos_doy")
microbial_cols <- grep("_otu$|^(Shannon|uniqueOTUs|PielouEvenness|Simpson|invSimpson|Rolling_)", names(df), value = TRUE)
microbial_community_cols =  c("su__sugar_degraders_fermenters", "aa__amino_acid_degraders", "c4__butyrate_and_valerate_degraders",  "fa__long_chain_fatty_acid_degraders",  "pro__propionate_degraders" ,  "ac__acetoclastic_methanogens" ,  "h2__hydrogenotrophic_methanogens")
diversity_cols=c("uniqueOTUs", "Shannon", "Simpson","invSimpson", "PielouEvenness", "Soerensen", "SoerensenAD1AD2", "BrayCurtisAD1AD2", "JaccardAD1AD2", "Rolling_Jaccard6days", "Rolling_Bray6days","Rate_of_change_Jaccard_per_day", "Rate_of_change_Bray_per_day")

FEATURES<-list(
  S = seasonal_cols,
  M_indv = microbial_cols,
  M_comm = microbial_community_cols,
  M_div = diversity_cols,
  M = c(microbial_cols, microbial_community_cols, diversity_cols),
  P = process_cols
)

FEATURE_SETS <- list(
  fs_S  = FEATURES$S,
  fs_M  = FEATURES$M,
  fs_SP = union(FEATURES$S, FEATURES$P),
  fs_SPM = union(FEATURES$S, union(FEATURES$P, FEATURES$M)),
  fs_P = FEATURES$P
)

ensure_output_dirs <- function() {
  dir.create("outputs", showWarnings = FALSE)
  dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
  dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
  dir.create("outputs/logs", recursive = TRUE, showWarnings = FALSE)
}


ID_COL <- "Digester"
TIME_COL <- "Sample_Date"
PRIMARY_DIGESTERS <- c("D1", "D2")
FORECAST_HORIZONS <- c(0, 7, 14, 21)
ASSESS_WIDTH_DAYS <- 14 # Make predictions for points in th enext two weeks
SCHEME_LOOKBACKS <- list(rolling = Inf, sliding = 56)
SMOOTH_WINDOW <- 6
PRIMARY_TARGETS <- c("A")


df_primary <- df_raw %>%
  dplyr::filter(.data[[ID_COL]] %in% PRIMARY_DIGESTERS) %>%
  add_seasonal_terms(TIME_COL)

df_aggregate <- aggregate_digesters(df_primary, ID_COL, TIME_COL)


df<-df_primary

debug_target <- "A"
debug_horizon <- 0
debug_use_smoothed <- FALSE
debug_scheme_name <- "rolling"
feature_sets <- FEATURE_SETS[1:2] # keep seasonal + operator baseline for speed
debug_nsplits <- 5
summary_window=6


df_debug <- df_aggregate
outcome_col <- debug_target
horizon_days <- debug_horizon
use_smoothed <- debug_use_smoothed
scheme <- "rolling"
lookback_days <- 56
assess_width_days <- ASSESS_WIDTH_DAYS
smooth_window <- SMOOTH_WINDOW
id_col <- ID_COL
time_col <- TIME_COL


message(
  sprintf(
    "debug_simple -> target: %s | horizon: %s | smoothed: %s | feature_sets: %s",
    debug_target,
    debug_horizon,
    debug_use_smoothed,
    paste(names(feature_sets), collapse = ", ")
  )
)

# Run the analysis after this point --------------------


# Setup the data given the config

df_ml <- df_debug %>%
  add_seasonal_terms(time_col) %>%
  make_forecast_targets(
    outcome_col = outcome_col,
    id_col = id_col,
    time_col = time_col,
    horizon_days = horizon_days,
    use_smoothed = use_smoothed,
    smooth_window = smooth_window)



if (!"target_last_value" %in% names(df_ml)) {
  df_ml <- df_ml %>%
    dplyr::arrange(.data[[id_col]], .data[[time_col]]) %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::mutate(target_last_value = .data[[outcome_col]]) %>%
    dplyr::ungroup()
}

message(
  sprintf(
    "run_forecast_experiment: outcome %s, horizon %s, id %s, time %s",
    outcome_col,
    horizon_days,
    id_col,
    time_col
  )
)

if (!nrow(df_ml)) {
  stop("Prepared modelling frame is empty; check outcome column and preprocessing steps.")
}

# Set up our croiss-validation
# Here we are using a sliding window approach (where window size = lookback)
# Rolling-orin approach by setting lookback_days = Inf
# - might need to set complete = FALSE for this?
message(sprintf("build samples from windows: lookback: %s, forecast_days: %s", lookback_days, assess_width_days))
resamples_sliding <- rsample::sliding_period(
  df_ml,
  index = !!rlang::sym(time_col), # Date or POSIXct column
  period = "day",
  lookback = lookback_days, # number of days in the analysis window
  assess_start = 1,
  assess_stop = assess_width_days, #number of days in the forecasting window
  step = 1,
  complete = TRUE
)
rset_class <- class(resamples_sliding)


#For debug,overight our definiiton and restrict to only a few windows
max_debug_splits <- min(debug_nsplits, nrow(resamples_sliding))
if (max_debug_splits < nrow(resamples_sliding)) {
  message(sprintf("Debug run: restricting to first %s of %s splits", max_debug_splits, nrow(resamples_sliding)))
  resamples_sliding <- resamples_sliding[seq_len(max_debug_splits), ]
  class(resamples_sliding) <- rset_class
}

# Set up temporal features ----------------

message("Preparing single recipe + models")
predictor_recipes <- sapply(names(feature_sets), function(x) {
  recipes::recipe(
    stats::as.formula(sprintf("%s ~ %s + %s + %s", outcome_col, paste(feature_sets[[x]], collapse = "+"), time_col, id_col)),
    data = df_ml
  ) %>%
    recipes::update_role(!!rlang::sym(time_col), new_role = "time_index", old_role = "predictor") %>%
    recipes::update_role(!!rlang::sym(id_col), new_role = "id", old_role = "predictor") %>%
    timetk::step_timeseries_signature(!!rlang::sym(time_col)) %>%
    timetk::step_fourier(!!rlang::sym(time_col), period = 7, K = 1) %>%
    recipes::step_lag(tidyselect::all_of(feature_sets[[x]]), lag = c(1, 3, 7, 14)) %>%
    # rolling mean / sd / max over ~2 weeks (@ 3 samples/week)
    timetk::step_slidify_augment(
      tidyselect::all_of(feature_sets[[x]]), period = summary_window,
      .f = ~mean(.x, na.rm = TRUE),
      align = "right", partial = TRUE, prefix = sprintf("roll%d_mean_",summary_window)
    ) %>%
    timetk::step_slidify_augment(
      tidyselect::all_of(feature_sets[[x]]), period = summary_window,
      .f = ~stats::sd(.x, na.rm = TRUE),
      align = "right", partial = TRUE, prefix = sprintf("roll%d_sd_",summary_window)
    ) %>%
    timetk::step_slidify_augment(
      tidyselect::all_of(feature_sets[[x]]), period = summary_window,
      .f = ~{
        val <- suppressWarnings(max(.x, na.rm = TRUE))
        if (is.finite(val)) val else NA_real_
      },
      align = "right", partial = TRUE, prefix = sprintf("roll%d_max_",summary_window)
    ) %>%
    recipes::step_rm(tidyselect::ends_with(".lbl")) %>% # remove ordered label columns
    recipes::step_impute_mean(all_numeric_predictors()) %>%
    recipes::step_impute_mode(all_nominal_predictors()) %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    recipes::step_corr(recipes::all_predictors(), threshold = 0.9) %>%
    recipes::step_normalize(recipes::all_predictors())
}, simplify = FALSE, USE.NAMES = TRUE)

# Define metrics
metrics_set <- yardstick::metric_set(yardstick::rmse, yardstick::mae, yardstick::rsq_trad)
ctrl_resamples <- tune::control_resamples(save_pred = TRUE)

fit_univariate_baseline <- function(spec, label) {
  purrr::imap_dfr(resamples_sliding$splits, function(split_obj, idx) {
    split_id <- resamples_sliding$id[[idx]]
    train_df <- rsample::analysis(split_obj)
    test_df <- rsample::assessment(split_obj)

    if (!nrow(train_df) || !nrow(test_df)) {
      return(NULL)
    }

    fit <- tryCatch(
      parsnip::fit(
        spec,
        stats::as.formula(sprintf(".y_target ~ %s", time_col)),
        data = train_df
      ),
      error = function(e) {
        message(sprintf("%s failed on %s: %s", label, split_id, e$message))
        return(NULL)
      }
    )

    if (is.null(fit)) {
      return(NULL)
    }

    preds <- predict(fit, new_data = test_df)
    metric_df <- metrics_set(
      data = tibble::tibble(truth = test_df$.y_target, .pred = preds$.pred),
      truth = truth,
      estimate = .pred
    ) %>%
      dplyr::mutate(
        model = label,
        id = split_id,
        .config = "baseline",
        penalty = NA_real_
      )

    metric_df
  })
}

# Define the models
models <- list()
models[["lm"]] <- parsnip::linear_reg() %>% set_engine("lm")

models[["lasso"]] <- linear_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet") %>%
  set_mode("regression")

lasso_grid <- grid_regular(penalty(range = c(-3, 0), trans = transform_log10()), levels = 7)

arima_spec <- modeltime::arima_reg() %>% parsnip::set_engine("auto_arima")
naive_spec <- modeltime::naive_reg() %>% parsnip::set_engine("naive")

rec_lastval <- recipes::recipe(
  stats::as.formula(sprintf("%s ~ target_last_value + %s + %s", outcome_col, time_col, id_col)),
  data = df_ml
) %>%
  recipes::update_role(!!rlang::sym(time_col), new_role = "time_index", old_role = "predictor") %>%
  recipes::update_role(!!rlang::sym(id_col), new_role = "id", old_role = "predictor") %>%
  recipes::step_impute_mean(recipes::all_predictors()) %>%
  recipes::step_normalize(recipes::all_predictors())

wf_lastval <- workflows::workflow() %>%
  workflows::add_recipe(rec_lastval) %>%
  workflows::add_model(parsnip::linear_reg() %>% parsnip::set_engine("lm"))

res_lastval <- tune::fit_resamples(
  wf_lastval,
  resamples = resamples_sliding,
  metrics = metrics_set,
  control = ctrl_resamples
)

wset <- workflow_set(
  preproc = predictor_recipes,
  models = c(models)
)

results <- wset %>%
  workflow_map(
    fn = "tune_grid", # This is the key
    resamples = resamples_sliding,
    grid = lasso_grid,
    metrics = metrics_set,
    control = ctrl_resamples,
    verbose = TRUE # Good for debugging
  )

select_best_safe <- purrr::possibly(function(x) tune::select_best(x, metric = "rmse"), otherwise = NULL)

## - Clean up after model fit (best config only for tuned members)
perf_outer <- results %>%
  dplyr::mutate(
    model = wflow_id,
    best_config = purrr::map(result, select_best_safe),
    metrics = purrr::map2(result, best_config, function(res_obj, best_tbl) {
      metrics_tbl <- tune::collect_metrics(res_obj, summarize = FALSE)
      if (!is.null(best_tbl) && ".config" %in% names(metrics_tbl)) {
        metrics_tbl <- dplyr::semi_join(metrics_tbl, best_tbl, by = ".config")
      }
      metrics_tbl
    })
  ) %>%
  dplyr::select(model, metrics) %>%
  tidyr::unnest(metrics) %>%
  dplyr::ungroup()

baseline_metrics <- dplyr::bind_rows(
  fit_univariate_baseline(arima_spec, "Model_ARIMA"),
  fit_univariate_baseline(naive_spec, "Model_Naive"),
  tune::collect_metrics(res_lastval, summarize = FALSE) %>% dplyr::mutate(model = "Model_LastObs")
)

perf_outer <- dplyr::bind_rows(perf_outer, baseline_metrics)
perf_outer <- perf_outer %>%
  dplyr::mutate(
    horizon = horizon_days,
    data_role = "assessment"
  )

message("Metrics compiled from held-out assessment folds; use summarize = FALSE for per-split diagnostics.")

str(perf_outer %>% filter(.metric == "mae") %>% group_by(model) %>% summarise(n()))

perf_outer %>%
  filter(.metric == "mae") %>%
  ggplot(aes(y = .estimate, x = model, fill = model)) +
  geom_boxplot() +
  theme_light(base_size = 20) +
  ylab("mae") +
  theme(axis.text.x = element_blank())

plot_metric_by_horizon <- function(metrics_df, metric_name = "rmse") {
  metrics_df %>%
    dplyr::filter(.metric == metric_name) %>%
    ggplot2::ggplot(ggplot2::aes(x = model, y = .estimate, colour = model)) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.4, show.legend = FALSE) +
    ggplot2::stat_summary(fun = median, geom = "point", size = 2) +
    ggplot2::labs(
      title = sprintf("%s performance by horizon", toupper(metric_name)),
      y = metric_name,
      x = "Model"
    ) +
    ggplot2::facet_wrap(~horizon, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

metric_plot <- plot_metric_by_horizon(perf_outer, metric_name = "rmse")
print(metric_plot)

debug_result <- list(
    wset = wset,
    resamples = resamples_sliding,
    results = results,
    metrics = perf_outer,
    baselines = list(
      arima = "Model_ARIMA",
      naive = "Model_Naive",
      last_obs = res_lastval
    ),
    plots = list(
      metric = metric_plot
    )
)

print(head(debug_result$metrics))
























#
#
# Manual debug ---------------------
# BG: Not sure if this works (after changes on Nov This might not work anymore)19/11/25) but
# gives an idea of how to simply set this up.

#--- 1. Define what you want to debug ---
wflow_id_to_debug <- "Model1_S_lm" # The workflow you want to check
wflow_id_to_debug <- wset$wflow_id[[2]] # The workflow you want to check
split_to_debug <- resamples_sliding$splits[[1]] # The first data split

# --- 2. Get the raw data from the split ---
train_data <- analysis(split_to_debug)
test_data  <- assessment(split_to_debug)

# Get the specific workflow from your workflow_set
# This is the correct replacement for `pull_workflow()`
single_wflow <- extract_workflow(wset, wflow_id_to_debug)

# Get the recipe from that single workflow
# This is the correct replacement for `pull_workflow_preprocessor()`
recipe_obj <- extract_preprocessor(single_wflow)

# Get the parsnip model spec from that single workflow
# This is the correct replacement for `pull_workflow_model()`
model_obj <- extract_spec_parsnip(single_wflow) %>% parsnip::set_args(penalty=1)

# --- 4. Manually prep and bake the data ---
# This is the most important debugging step

# "prep" learns the recipe steps (e.g., means, SDs) from the training data
prepped_recipe <- prep(recipe_obj, training = train_data)

# "bake" applies the learned steps to new data
# Use new_data = NULL to get the processed *training* data
baked_train_data <- bake(prepped_recipe, new_data = NULL)
baked_test_data  <- bake(prepped_recipe, new_data = test_data)

manual_fit <- fit(model_obj,
                  formula = formula(prepped_recipe),
                  data = baked_train_data)

# Predict on the *processed* testing data
yp <- predict(manual_fit, new_data = baked_test_data)
