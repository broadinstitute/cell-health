# Utility functions to import into dose response curve analysis

zero_one_normalize <- function(x){(x - min(x)) / (max(x) - min(x))}


get_curve_fit <- function(moa_df, dose_df, compound, model) {
  moa_focus_df <- moa_df %>%
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

  # Use drc to fit the Hill equation
  fit <- drc::drm(
    model_score_transform ~ Metadata_mmoles_per_liter,
    data = moa_focus_df,
    fct = drc::LL2.4(
      fixed = c(NA, NA, NA, NA),
      names = c("b", "c", "d", "e")
    )
  )

  dose_focus_df <- dose_df %>%
    dplyr::filter(compound == !!compound,
                  model == !!model)

  return(list("dose" = dose_focus_df, "fit" = fit, "moa" = moa_focus_df))
}

get_dose_curve <- function(moa_long_df, dose_df, model, pert_name, cell_health_model) {
    compound <- unique(
    dose_df %>%
        dplyr::filter(pert_iname == !!pert_name) %>%
        dplyr::pull(compound)
    )[1]

    example_curve <- get_curve_fit(moa_long_df, dose_df, compound, model)

    # Sample data
    newdata <- expand.grid(conc=exp(seq(log(0.04), log(10), length=1000)))
    # predictions and confidence intervals
    pm <- stats::predict(example_curve$fit, newdata=newdata, level = 0.95, interval="confidence")
    newdata$p <- pm[,1]
    newdata$pmin <- pm[,2]
    newdata$pmax <- pm[,3]

    dose_curve_gg <- ggplot(example_curve$moa,
                            aes(x = Metadata_mmoles_per_liter, y = model_score_transform)) +
        geom_point() +
        coord_trans(x="log10") +
        geom_ribbon(data=newdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
        geom_line(data=newdata, aes(x=conc, y=p)) +
        theme_bw() +
        xlab(paste0("Micromoles per Liter\n", pert_name)) +
        ylab(paste0("Cell Health Model\n", cell_health_model))

    return(dose_curve_gg)
}
