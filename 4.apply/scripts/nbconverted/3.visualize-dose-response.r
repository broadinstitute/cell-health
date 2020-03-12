suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source(file.path("scripts", "dose_utils.R"))

figure_dir <- file.path("figures", "dose_response")

moa_file <- file.path("repurposing_cellhealth_shiny", "data", "moa_cell_health_modz.tsv.gz")

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

print(dim(moa_long_df))
head(moa_long_df, 3)

dose_file <- file.path("results", "dose_response_curve_fit_results.tsv")

dose_cols <- readr::cols(
    .default = readr::col_character(),
    slope = readr::col_double(),
    slope_error = readr::col_double(),
    slope_t = readr::col_double(),
    slope_p = readr::col_double(),
    ic_fifty = readr::col_double(),
    ic_fifty_error = readr::col_double(),
    ic_fifty_t = readr::col_double(),
    ic_fifty_p = readr::col_double(),
    residual = readr::col_double()
)

dose_df <- readr::read_tsv(dose_file, col_types = dose_cols)

print(dim(dose_df))
head(dose_df, 3)

residual_gg <- ggplot(dose_df, aes(x = residual)) +
    geom_density(fill="grey") +
    theme_bw() +
    xlab("Residual of Dose Curve Fit") +
    ylab("Density")

summary_gg <- ggplot(dose_df, aes(x = slope_t, y = -log10(slope_p))) +
    geom_point(alpha = 0.5, size = 0.3, pch = 19, color = "black", fill = "black") +
    theme_bw() +
    xlab("Slope Fit Statistic") +
    ylab("-log10 Slope p Value")

ic_fifty_p_density_gg <- ggplot(dose_df, aes(x = -log10(ic_fifty_p))) +
    geom_density(fill = "grey") +
    theme_bw() +
    coord_flip() +
    xlab("") +
    ylab("Density")

output_file <- file.path(figure_dir, "dose_summary.png")

dose_summary_gg <- cowplot::plot_grid(
    residual_gg,
    cowplot::plot_grid(
        summary_gg,
        ic_fifty_p_density_gg,
        ncol = 2,
        nrow = 1,
        align = "h",
        rel_widths = c(0.8, 0.2)
    ),
    nrow = 2,
    align = "h",
    rel_heights = c(0.4, 0.6)
)

cowplot::save_plot(output_file, dose_summary_gg, base_height = 5, base_width = 5)

dose_summary_gg

dose_curve_theme <- theme(
    axis.text.x = element_text(size = 4),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 8)
)

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
        geom_point(size = 0.5) +
        coord_trans(x="log10") +
        geom_ribbon(data=newdata, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +
        geom_line(data=newdata, aes(x=conc, y=p), lwd = 0.5) +
        theme_bw() +
        xlab(paste0("Micromoles per Liter\n", pert_name)) +
        ylab(paste0("Cell Health Model\n", cell_health_model))

    return(dose_curve_gg)
}

model <- "cell_health_modz_target_vb_ros_mean"
pert_name <- "bortezomib"
bortezomib_gg <- get_dose_curve(moa_long_df, dose_df, model, pert_name, "ROS Mean") + dose_curve_theme

model <- "cell_health_modz_target_cc_g1_n_spots_mean"
cell_health_model <- "DNA Damage in G1 Cells"
pert_name <- "MLN-4924"
mln_g1_gg <- get_dose_curve(moa_long_df, dose_df, model, pert_name, cell_health_model) + dose_curve_theme

model <- "cell_health_modz_target_cc_edu_pos_n_objects"
cell_health_model <- "Number of Proliferating Cells"
pert_name <- "midostaurin"
midostaurin_edu_gg <- get_dose_curve(moa_long_df, dose_df, model, pert_name, cell_health_model) + dose_curve_theme

output_file <- file.path(figure_dir, "select_dose_examples.png")

dose_examples_gg <- cowplot::plot_grid(
    bortezomib_gg,
    mln_g1_gg,
    midostaurin_edu_gg,
    nrow = 1,
    ncol = 3,
    align = "v"
)

cowplot::save_plot(output_file, dose_examples_gg, base_height = 3, base_width = 8.5)

dose_examples_gg
