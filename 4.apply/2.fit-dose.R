# Fit Hill Equations (4 parameter) to Drug Repurposing cell health model predictions
#
# The Drug Repurposing Project measured Cell Painting data for 1,570 compounds across
# six dose points. One way to demonstrate accurate model predictions, is to observe
# the goodness of fit of dose response curves.
#
# We use the CRAN package `drc` https://cran.r-project.org/web/packages/drc/index.html
# for the dose response curve analysis
#
# The following script will fit a dose response curve for every model for every compound
#
# Usage:
#   Rscript --vanilla 2.fit-dose.R
#
# Output:
#   Summary table containing goodness of fit dose response curve estimates

library(drc)
library(dplyr)
library(ggplot2)
source(file.path("scripts", "dose_utils.R"))


moa_file <- file.path(
  "repurposing_cellhealth_shiny",
  "data",
  "moa_cell_health_modz.tsv.gz"
)

moa_cols <- readr::cols(
  .default = readr::col_double(),
  Image_Metadata_Well = readr::col_character(),
  Metadata_broad_core_id = readr::col_character(),
  Metadata_broad_sample = readr::col_character(),
  Metadata_dose_recode = readr::col_integer(),
  pert_id = readr::col_character(),
  pert_iname = readr::col_character(),
  pert_type = readr::col_character(),
  moa = readr::col_character()
)

moa_df <- readr::read_tsv(moa_file, col_types = moa_cols)
moa_long_df <- moa_df %>% reshape2::melt(id.vars = c(
  "Image_Metadata_Well",
  "Metadata_broad_core_id",
  "Metadata_broad_sample",
  "Metadata_dose_recode",
  "Metadata_mmoles_per_liter",
  "umap_x",
  "umap_y",
  "pert_id",
  "pert_iname",
  "pert_type",
  "moa"),
  variable.name = "model",
  value.name = "model_score"
)

all_compounds <- unique(moa_long_df$Metadata_broad_sample)
all_models <- unique(moa_long_df$model)

result_list <- list()
for (model in all_models) {
  result_list[[model]] <- list()

  for (compound in all_compounds) {

    if (compound == "DMSO") {
      next
    }

    dose_focus_df <- moa_long_df %>%
      dplyr::filter(Metadata_broad_sample == !!compound,
                    model == !!model) %>%
      dplyr::select(
        Metadata_mmoles_per_liter,
        pert_id,
        pert_iname,
        moa,
        model,
        model_score
      ) %>%
      dplyr::mutate(model_score_transform = zero_one_normalize(model_score))

    pert_id <- dose_focus_df$pert_id[1]
    pert_iname <- dose_focus_df$pert_iname[1]
    moa <- dose_focus_df$moa[1]

    if (nrow(dose_focus_df) < 4) {
      model_result <- c(
          compound, model, pert_id, pert_iname, moa, rep(NA, 9), "not_enough_dose"
      )
      result_list[[model]][[compound]] <- model_result
      next
    }

    tryCatch(
      expr = {
        fit <- drc::drm(
            model_score_transform ~ Metadata_mmoles_per_liter,
            data = dose_focus_df,
            fct = LL2.4(
                fixed = c(NA, NA, NA, NA),
                names = c("b", "c", "d", "e")
            )
        )
        model_output_status <- "fit"
      },
      error = function(e) {
        model_output_status <- "fail_convergence"
      }
    )

    if (model_output_status == "fail_convergence") {
      model_result <- c(
          compound, model, pert_id, pert_iname, moa, rep(NA, 9), model_output_status
        )
      result_list[[model]][[compound]] <- model_result
      next
    }

    fit_summary <- summary(fit)
    residual_rse <- fit_summary$resVar

    summary_coef <- fit_summary$coefficients
    rownames(summary_coef) <- c("slope", "bottom", "top", "ic50")

    slope <- summary_coef["slope", "Estimate"]
    slope_error <- summary_coef["slope", "Std. Error"]
    slope_t <- summary_coef["slope", "t-value"]
    slope_p <- summary_coef["slope", "p-value"]

    ic_fifty <- summary_coef["ic50", "Estimate"]
    ic_fifty_error <- summary_coef["ic50", "Std. Error"]
    ic_fifty_t <- summary_coef["ic50", "t-value"]
    ic_fifty_p <- summary_coef["ic50", "p-value"]

    # Build a dataframe for each compound for each cell health target
    model_result <- c(
        compound,
        model,
        pert_id,
        pert_iname,
        moa,
        slope,
        slope_error,
        slope_t,
        slope_p,
        ic_fifty,
        ic_fifty_error,
        ic_fifty_t,
        ic_fifty_p,
        residual_rse,
        model_output_status
    )

    result_list[[model]][[compound]] <- model_result
  }

  result_list[[model]] <- do.call(rbind, result_list[[model]]) %>%
    dplyr::as_tibble(.name_repair = "minimal")

  colnames(result_list[[model]]) <- c(
      "compound",
      "model",
      "pert_id",
      "pert_iname",
      "moa",
      "slope",
      "slope_error",
      "slope_t",
      "slope_p",
      "ic_fifty",
      "ic_fifty_error",
      "ic_fifty_t",
      "ic_fifty_p",
      "residual",
      "status"
  )
}

full_results_df <- do.call(rbind, result_list)
full_results_df$slope <- as.numeric(paste(full_results_df$slope))
full_results_df$slope_error <- as.numeric(paste(full_results_df$slope_error))
full_results_df$slope_t <- as.numeric(paste(full_results_df$slope_t))
full_results_df$slope_p <- as.numeric(paste(full_results_df$slope_p))
full_results_df$ic_fifty <- as.numeric(paste(full_results_df$ic_fifty))
full_results_df$ic_fifty_error <- as.numeric(paste(full_results_df$ic_fifty_error))
full_results_df$ic_fifty_t <- as.numeric(paste(full_results_df$ic_fifty_t))
full_results_df$ic_fifty_p <- as.numeric(paste(full_results_df$ic_fifty_p))
full_results_df$residual <- as.numeric(paste(full_results_df$residual))

# Output results
output_file <- file.path("results", "dose_response_curve_fit_results.tsv")

full_results_df %>%
    dplyr::arrange(residual) %>%
    dplyr::mutate(ic_fifty_transform = 2 ** ic_fifty) %>%
    readr::write_tsv(output_file)
