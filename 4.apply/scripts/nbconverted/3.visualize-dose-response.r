suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source(file.path("repurposing_cellhealth_shiny", "dose_utils.R"))

figure_dir <- file.path("figures", "dose_response")

dose_curve_theme <- theme(
    axis.text.x = element_text(size = 4),
    axis.text.y = element_text(size = 5),
    axis.title = element_text(size = 7)
)

moa_file <- file.path(
    "repurposing_cellhealth_shiny", "data", "moa_cell_health_modz.tsv.gz"
)

moa_cols <- readr::cols(
  .default = readr::col_double(),
  Metadata_Plate_Map_Name = readr::col_character(),
  Metadata_broad_core_id = readr::col_character(),
  Metadata_broad_sample = readr::col_character(),
  Metadata_pert_well = readr::col_character(),
  Metadata_dose_recode = readr::col_character(),
  Metadata_mmoles_per_liter = readr::col_double(),
  broad_id = readr::col_character(),
  pert_iname = readr::col_character(),
  InChIKey14 = readr::col_character(),
  moa = readr::col_character(),
  target = readr::col_character(),
  broad_date = readr::col_character(),
  clinical_phase = readr::col_character(),
  alternative_moa = readr::col_character(),
  alternative_target = readr::col_character()
)

moa_df <- readr::read_tsv(moa_file, col_types = moa_cols)
moa_long_df <- moa_df %>% reshape2::melt(id.vars = c(
  "Metadata_Plate_Map_Name",
  "Metadata_pert_well",
  "Metadata_broad_core_id",
  "InChIKey14",
  "Metadata_broad_sample",
  "Metadata_dose_recode",
  "Metadata_mmoles_per_liter",
  "umap_x",
  "umap_y",
  "broad_id",
  "broad_date",
  "clinical_phase",
  "pert_iname",
  "moa",
  "target",
  "alternative_moa",
  "alternative_target"),
  variable.name = "model",
  value.name = "model_score"
)

print(dim(moa_long_df))
head(moa_long_df, 3)

dose_file <- file.path("repurposing_cellhealth_shiny", "data", "dose_response_curve_fit_results.tsv")

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

dose_df <- readr::read_tsv(dose_file, col_types = dose_cols) %>% tidyr::drop_na()

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

model <- "cell_health_modz_target_vb_ros_mean"
pert_name <- "bortezomib"
bortezomib_gg <- suppressWarnings(
    get_dose_curve(moa_long_df, dose_df, model, pert_name, "ROS Mean") + dose_curve_theme
)

dose_df %>% dplyr::filter(pert_iname == !!pert_name, model == !!model)

model <- "cell_health_modz_target_cc_g1_n_objects"
cell_health_model <- "G1 - Number of Objects"
pert_name <- "HMN-214"
hmn214_g1_gg <- suppressWarnings(
    get_dose_curve(moa_long_df, dose_df, model, pert_name, cell_health_model) + dose_curve_theme
)

dose_df %>% dplyr::filter(pert_iname == !!pert_name, model == !!model)

output_file <- file.path(figure_dir, "select_dose_examples.png")

dose_examples_gg <- cowplot::plot_grid(
    bortezomib_gg,
    hmn214_g1_gg,
    nrow = 2,
    ncol = 1,
    align = "v"
)

cowplot::save_plot(output_file, dose_examples_gg, base_height = 3.25, base_width = 3)

dose_examples_gg
