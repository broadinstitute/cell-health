suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))

source(file.path("scripts", "assay_themes.R"))

consensus <- "modz"

results_dir <- "results"
figure_dir <- file.path("figures", "regression", consensus)
individual_fig_dir <- file.path(
    "figures",
    "individual_target_performance",
    "regression",
    consensus
)

dir.create(results_dir, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_fig_dir, recursive = TRUE, showWarnings = FALSE)

# Complete regression results
regression_file <- file.path(
    results_dir, 
    paste0("full_cell_health_regression_", consensus, ".tsv.gz")
)
all_regression_metrics_df <- readr::read_tsv(regression_file, col_types = readr::cols())

# Model coefficients
coef_file <- file.path(
    results_dir,
    paste0("full_cell_health_coefficients_", consensus, ".tsv.gz")
)
full_coef_df <- readr::read_tsv(coef_file, col_types = readr::cols()) %>%
    dplyr::filter(y_transform == "raw")

# Metadata information
metadata_file <- file.path(
    "..",
    "1.generate-profiles",
    "data",
    "profile_id_metadata_mapping.tsv"
)
metadata_df <- readr::read_tsv(metadata_file, col_types = readr::cols())

# Ground truth (y) values
y_file <- file.path(
    results_dir,
    paste0("full_cell_health_y_labels_", consensus, ".tsv.gz")
)
y_df <- readr::read_tsv(y_file, col_types = readr::cols()) %>%
    dplyr::filter(y_transform == "raw")

# Label variables with specific cell health classes
label_file <- file.path(
    "..",
    "1.generate-profiles",
    "data",
    "labels",
    "feature_mapping_annotated.csv"
)
label_df <- readr::read_csv(label_file, col_types = readr::cols())

head(label_df)

# Combine data for downstream processing
y_binary_subset_true_df <- y_df %>%
    dplyr::filter(y_type == "y_true")

y_binary_subset_pred_df <- y_df %>%
    dplyr::filter(y_type == "y_pred")

# Process data for plotting
y_plot_df <- y_binary_subset_true_df %>%
    dplyr::inner_join(
        y_binary_subset_pred_df,
        by = c("Metadata_profile_id",
               "target",
               "data_type",
               "shuffle",
               "y_transform"),
        suffix = c("_true", "_pred")) %>%
    dplyr::left_join(metadata_df, by = "Metadata_profile_id")

y_plot_df$data_type <- dplyr::recode(
    y_plot_df$data_type,
    "train" = "Train",
    "test" = "Test"
)

print(dim(y_plot_df))
head(y_plot_df, 3)

all_regression_metrics_df$data_fit <- dplyr::recode(
    all_regression_metrics_df$data_fit,
    "train" = "Train",
    "test" = "Test"
)

all_regression_metrics_df$shuffle <- dplyr::recode(
    all_regression_metrics_df$shuffle,
    "shuffle_true" = "Shuffle",
    "shuffle_false" = "Real"
)

all_regression_metrics_df <- all_regression_metrics_df %>%
    dplyr::rename(data_type = data_fit)

all_regression_metrics_df$value <- round(all_regression_metrics_df$value, 2)

regression_metrics_df <- all_regression_metrics_df %>%
    dplyr::filter(cell_line == "all")

print(dim(regression_metrics_df))
head(regression_metrics_df, 3)

mse_df <- regression_metrics_df %>%
    dplyr::filter(metric == "mse",
                  y_transform == "raw",
                  data_type == "Test") %>%
    dplyr::left_join(label_df, by = c("target" = "id"))

# Take absolute value of mean squared error
# see https://github.com/scikit-learn/scikit-learn/issues/2439
mse_df$value = abs(mse_df$value)

# Recode metadata measurent to other
mse_df$measurement <- dplyr::recode(mse_df$measurement, "metadata" = "other")

# Sort mse by minimum MSE in test set
target_order <- mse_df %>%
    dplyr::filter(data_type == "Test",
                  shuffle == "Real") %>%
    dplyr::arrange(desc(value)) %>%
    dplyr::select(readable_name)

mse_df$readable_name <- factor(
    mse_df$readable_name,
    levels=target_order$readable_name
)

print(dim(mse_df))
head(mse_df, 4)

ggplot(mse_df,
       aes(x = readable_name,
           y = value)) +
    geom_bar(aes(fill = measurement),
             stat = "identity",
             alpha = 0.7,
             position = position_dodge()) +
    facet_grid(~shuffle, scales="free_y") +
    coord_flip() +
    theme_bw() +
    ylab("Mean Squared Error") +
    xlab("Target") +
    scale_fill_manual(name = "Measurement",
                      values = measurement_colors,
                      labels = measurement_labels) +
    theme(axis.text.x = element_text(size = 8, angle = 90),
          axis.text.y = element_text(size = 5.5),
          axis.title = element_text(size = 10),
          legend.position = "right",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

file <- file.path(
    figure_dir,
    paste0("mse_test_summary_", consensus, ".png")
)

ggsave(file, dpi = 500, width = 7, height = 6)

# Split shuffle column for scatter plot
mse_spread_df <- mse_df %>% tidyr::spread(shuffle, value)

head(mse_spread_df, 2)

ggplot(mse_spread_df,
       aes(x = Real,
           y = Shuffle,
           color = measurement)) +
    geom_abline(intercept = 0,
                lwd = 0.1,
                slope = 1,
                linetype = "dotted",
                alpha = 0.7,
                color = "red") +
    geom_point(size = 0.4,
               alpha = 0.7,
               pch = 16) +
    xlim(c(0, 0.7)) +
    ylim(c(0, 1.7)) +
    xlab("Real Data (Mean Squared Error)") +
    ylab("Permuted Data (Mean Squared Error)") +
    ggtitle("Regression Performance") +
    theme_bw() +
    dye_theme +
    scale_color_manual(
        name = "Cell Health\nPhenotypes",
        values = measurement_colors,
        labels = measurement_labels
    )

output_file <- file.path(
    figure_dir,
    paste0("mse_comparison_scatter_", consensus, ".png")
)
ggsave(output_file, width = 2, height = 1.5, dpi = 500, units = "in")

r2_df <- regression_metrics_df %>%
    dplyr::filter(metric == "r_two",
                  y_transform == "raw") %>%
    dplyr::left_join(label_df, by = c("target" = "id"))

# Sort mse by r2 in test set
target_order <- r2_df %>%
    dplyr::filter(data_type == "Test",
                  shuffle == "Real") %>%
    dplyr::arrange(value) %>%
    dplyr::select(readable_name)

r2_df$readable_name <- factor(r2_df$readable_name, levels=target_order$readable_name)

# Recode metadata measurent to other
r2_df$measurement <- dplyr::recode(r2_df$measurement, "metadata" = "other")

r2_df$measurement <- factor(r2_df$measurement, levels = c(
    "cell_viability",
    "death",
    "apoptosis",
    "ros",
    "dna_damage",
    "g1_phase",
    "s_phase",
    "g2_phase",
    "early_mitosis",
    "mitosis",
    "late_mitosis",
    "cell_cycle_count",
    "shape",
    "other"
))

head(r2_df, 4)

rsquared_bar_gg <- ggplot(r2_df %>% dplyr::filter(shuffle == "Real"),
       aes(x = readable_name,
           y = value)) +
    geom_bar(aes(fill = measurement),
             stat = "identity",
             alpha = 0.7,
             position = position_dodge()) +
    facet_grid(~data_type, scales = "free_y") +
    coord_flip() +
    theme_bw() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(bquote(R^2)) +
    xlab("Cell Health Target") +
    scale_fill_manual(name = "Measurement",
                      values = measurement_colors,
                      labels = measurement_labels) +
    theme(axis.text.x = element_text(size = 8, angle = 90),
          axis.text.y = element_text(size = 5.5),
          axis.title = element_text(size = 10),
          legend.position = "right",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

file <- file.path(
    figure_dir,
    paste0("r_squared_model_summary_", consensus, ".png")
)
ggsave(file, rsquared_bar_gg, dpi = 500, width = 7, height = 6)

rsquared_bar_gg

# Process data for cell line specific plot
cellline_compare_regression_df <- all_regression_metrics_df %>%
    dplyr::filter(cell_line != "all",
                  data_type == "Test",
                  metric == "r_two",
                  shuffle == "Real") %>%
    dplyr::mutate(outlier = ifelse(value < 0, "R2: < 0", "R2: 0 to 1")) %>%
    dplyr::left_join(label_df, by = c("target" = "id"))

cellline_compare_regression_df$measurement <-
    dplyr::recode(cellline_compare_regression_df$measurement, "metadata" = "other")

cellline_compare_regression_df$readable_name <- factor(
    cellline_compare_regression_df$readable_name,
    levels=target_order$readable_name
)

head(cellline_compare_regression_df, 2)

rsquared_bar_cellline_gg <- ggplot(cellline_compare_regression_df,
       aes(x = value,
           y = readable_name,
           color = cell_line)) +
    geom_point(alpha = 0.8) +
    theme_bw() +
    xlab(bquote("Test Set "~R^2~"")) +
    ylab("Cell Health Target") +
    facet_wrap(~outlier, scales="free_x") +
    scale_color_manual(
        name = "Cell Line",
        labels = cell_line_labels,
        values = cell_line_colors
    ) +
    scale_y_discrete(position = "right") +
    theme(axis.text.x = element_text(size = 8, angle = 90),
          axis.text.y = element_text(size = 5.5),
          axis.title = element_text(size = 10),
          legend.position = c(0.15, 0.925),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.key = element_rect(size = 4),
          legend.key.size = unit(0.4, 'lines'),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

output_file <- file.path(
    figure_dir,
    paste0("cell_line_differences_rsquared_", consensus, ".png")
)
ggsave(output_file, rsquared_bar_cellline_gg, dpi = 500, width = 7, height = 6)
rsquared_bar_cellline_gg

regression_legend <- cowplot::get_legend(rsquared_bar_gg)

left_panel_margin <- margin(l = -0.8, r = 0.5, t = 0.2, b = 0.2, unit = "cm")
right_panel_margin <- margin(l = 0, r = -0.8, t = 0.2, b = 0.2, unit = "cm")

regression_performance_figure <- cowplot::plot_grid(
    rsquared_bar_gg +
        theme(legend.position = 'none',
              plot.margin = left_panel_margin) + xlab(""),
    regression_legend,
    rsquared_bar_cellline_gg +
        theme(plot.margin = right_panel_margin) + ylab(""),
    labels = c("a", "", "b"),
    ncol = 3,
    nrow = 1,
    align = "h",
    hjust = c(-0.5, -0.5, 3),
    rel_widths = c(1, 0.3, 1)
)

output_file <- file.path(
    figure_dir,
    paste0("regression_performance_figure_", consensus, ".png")
)
cowplot::save_plot(output_file, regression_performance_figure, base_width = 10, base_height = 6)

regression_performance_figure

concise_fig_df <- r2_df %>%
    dplyr::filter(shuffle == "Real", data_type == "Test") %>%
    dplyr::select(value, metric, target, data_type, shuffle, y_transform, cell_line, readable_name, measurement)

concise_fig_df$data_fit <- dplyr::recode(
    concise_fig_df$data_type,
    "Test" = "Test Set Performance\n(Data not used in training)"
)

head(concise_fig_df, 2)

cell_line_focus_df <- cellline_compare_regression_df %>%
    dplyr::select(value, metric, target, data_type, shuffle, y_transform, cell_line, readable_name, measurement) %>%
    dplyr::distinct() %>%
    dplyr::mutate(truncated_value = ifelse(value < -1, -1, value))

cell_line_focus_df$data_fit <- dplyr::recode(
    cell_line_focus_df$data_type,
    "Test" = "Test Set Performance\n(Data not used in training)"
)

head(cell_line_focus_df, 2)

measurement_legend <- cowplot::get_legend(
    ggplot(concise_fig_df,
       aes(x = readable_name)) +
    geom_bar(aes(fill = measurement, y = value),
             stat = "identity",
             alpha = 0.7,
             position = position_dodge()) +
    facet_grid(~data_fit, scales = "free_y") +
    theme_bw() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(bquote(R^2)) +
    xlab("Cell Health Target") +
    scale_fill_manual(name = "Measurement",
                      values = measurement_colors,
                      labels = measurement_labels) +
    theme(axis.text.x = element_text(size = 9, angle = 90),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 10),
          legend.position = "right",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9),
          strip.text = element_text(size = 9),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))
)

cell_line_legend <- cowplot::get_legend(
    ggplot(concise_fig_df,
       aes(x = readable_name)) +
    geom_point(data = cell_line_focus_df, aes(shape = cell_line, y = truncated_value, fill = measurement)) +
    facet_grid(~data_fit, scales = "free_y") +
    theme_bw() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(bquote(R^2)) +
    xlab("Cell Health Target") +
    scale_fill_manual(name = "Measurement",
                      values = measurement_colors,
                      labels = measurement_labels) +
    scale_shape_manual(name = "Cell Line",
                       values = c("A549" = 21, "ES2" = 24, "HCC44" = 22),
                       labels = c("A549" = "A549", "ES2" = "ES2", "HCC44" = "HCC44")) +
    theme(axis.text.x = element_text(size = 9, angle = 90),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 10),
          legend.position = "right",
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9),
          strip.text = element_text(size = 9),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    guides(fill = FALSE,
           shape = guide_legend(order = 1))
)

concise_legend_gg <- cowplot::plot_grid(
    cell_line_legend,
    cowplot::ggdraw(),
    measurement_legend,
    cowplot::ggdraw(),
    nrow = 4,
    align = "v",
    axis = "l",
    rel_heights = c(1, 0.3, 1, 0.7)
)

concise_gg <- ggplot(concise_fig_df,
       aes(x = readable_name)) +
    geom_bar(aes(fill = measurement, y = value),
             stat = "identity",
             alpha = 0.7,
             color = "black",
             lwd = 0.2,
             position = position_dodge()) +
    geom_point(data = cell_line_focus_df, aes(shape = cell_line, y = truncated_value, fill = measurement)) +
    facet_grid(~data_fit, scales = "free_y") +
    theme_bw() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(bquote(R^2)) +
    xlab("Cell Health Target") +
    scale_fill_manual(name = "Measurement",
                      values = measurement_colors,
                      labels = measurement_labels) +
    scale_shape_manual(name = "Cell Line",
                       values = c("A549" = 21, "ES2" = 24, "HCC44" = 22),
                       labels = c("A549" = "A549", "ES2" = "ES2", "HCC44" = "HCC44")) +
    theme(axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 10),
          legend.position = "none",
          strip.text = element_text(size = 9),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    scale_y_continuous(
        limits = c(-1, 1),
        breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
        labels = c("\u2264 -1.0", "-0.5", "0", "0.5", "1.0")
        )

output_file <- file.path(
    figure_dir,
    paste0("concise_regression_summary_", consensus, ".png")
)

concise_gg_full <- cowplot::plot_grid(
    concise_gg,
    concise_legend_gg,
    ncol = 2,
    rel_widths = c(1, 0.3)
)

cowplot::save_plot(output_file, concise_gg_full, dpi = 600, base_width = 7, base_height = 6)

concise_gg_full

# Split shuffle column for scatter plot
r2_spread_df <- r2_df %>% tidyr::spread(shuffle, value)
r2_spread_df$data_type <- factor(r2_spread_df$data_type, levels = c("Train", "Test"))

head(r2_spread_df, 2)

ggplot(r2_spread_df,
       aes(x = Real,
           y = Shuffle,
           color = measurement)) +
    geom_abline(intercept = 0,
                lwd = 0.5,
                slope = 1,
                linetype = "dashed",
                alpha = 0.7,
                color = "red") +
    geom_point(size = 1,
               alpha = 0.7,
               pch = 16) +
    xlab("Real Data R-Squared") +
    ylab("Permuted Data R-Squared") +
    ggtitle("Regression Performance") +
    theme_bw() +
    facet_wrap("~data_type") +
    coord_fixed() +
    scale_color_manual(
        name = "Cell Health\nPhenotypes",
        values = measurement_colors,
        labels = measurement_labels
    ) +
    dye_theme +
    theme(axis.text.x = element_text(size = 5.5),
          axis.text.y = element_text(size = 5.5),
          axis.title = element_text(size = 8),
          strip.text = element_text(size = 7),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

output_file <- file.path(
    figure_dir,
    paste0("r_squared_comparison_scatter_", consensus, ".png")
)
ggsave(output_file, width = 5.5, height = 3, dpi = 600, units = "in")

good_example_models <- c(
    "vb_percent_dead",
    "cc_s_n_objects",
    "cc_g1_n_spots_h2ax_mean",
    "vb_percent_all_apoptosis"
)

bad_example_models <- c(
    "cc_polynuclear_n_spots_h2ax_mean",
    "cc_cc_late_mitosis"
)

plot_list <- list()
for (good_model in c(good_example_models, bad_example_models)) {
    
    sup_fig_good_subset_df <- y_plot_df %>%
        dplyr::filter(
            target == !!good_model,
            y_transform == "raw"
        )
    sup_fig_good_subset_df$shuffle <- sup_fig_good_subset_df$shuffle %>%
        dplyr::recode_factor(
            "shuffle_false" = "Real",
            "shuffle_true" = "Permuted"
        )
    
    r2_df <- regression_metrics_df %>%
    dplyr::filter(
        metric == "r_two",
        target == !!good_model,
        y_transform == "raw"
    ) %>%
    dplyr::rename(r2 = value) %>%
    dplyr::mutate(x = max(sup_fig_good_subset_df$recode_target_value_true) -
                  sd(sup_fig_good_subset_df$recode_target_value_true) * 3,
                  y = max(sup_fig_good_subset_df$recode_target_value_pred) -
                  sd(sup_fig_good_subset_df$recode_target_value_true) * 2.5)

    r2_df$shuffle <- r2_df$shuffle %>%
            dplyr::recode_factor(
                "Shuffle" = "Permuted"
            )
    
    readable_title <- label_df %>%
        dplyr::filter(id == !!good_model) %>%
        dplyr::pull(readable_name)
    
    plot_list[[good_model]] <- ggplot(sup_fig_good_subset_df,
           aes(x = recode_target_value_true,
               y = recode_target_value_pred)) +
        geom_smooth(method = 'lm', formula = y~x) +
        geom_point(aes(color = Metadata_cell_line),
                   size = 0.5,
                   alpha = 0.5) +
        facet_grid(data_type~shuffle) +
        ggtitle(readable_title) +
        theme_bw() +
        xlab("True Values") +
        ylab("Predicted Values") +
        scale_color_manual(
            name = "Cell Line",
            labels = cell_line_labels,
            values = cell_line_colors
        ) +
        geom_rect(data = sup_fig_good_subset_df %>%
                    dplyr::filter(data_type == "Test", shuffle == "Real") %>%
                    dplyr::distinct(data_type, shuffle, .keep_all = TRUE), 
                  fill = NA, alpha = 1, colour = "red", linetype = "solid", size = 2,
                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
        geom_text(data = r2_df, size = 3, aes(label = paste("R2 =", r2), x = x, y = y)) +
        theme(strip.text = element_text(size = 10),
              strip.background = element_rect(colour = "black",
                                              fill = "#fdfff4"))
}

# Compile and save figure
cell_line_legend <- get_legend(
  plot_list[["vb_percent_dead"]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)

top_row_gg <- cowplot::plot_grid(
    plot_list[["vb_percent_dead"]] + theme(legend.position = "none"),
    plot_list[["cc_s_n_objects"]] + theme(legend.position = "none"),
    nrow = 1,
    align = "h"
)

middle_row_gg <- cowplot::plot_grid(
    plot_list[["cc_g1_n_spots_h2ax_mean"]] + theme(legend.position = "none"),
    plot_list[["vb_percent_all_apoptosis"]] + theme(legend.position = "none"),
    nrow = 1,
    align = "h"
)

sup_fig_part_a <- cowplot::plot_grid(
    top_row_gg,
    middle_row_gg,
    nrow = 2,
    align = "h"
)

sup_fig_part_a <- cowplot::plot_grid(
    sup_fig_part_a,
    cell_line_legend,
    rel_widths = c(3, 0.4),
    ncol = 2
)

bottom_row_gg <- cowplot::plot_grid(
    plot_list[["cc_polynuclear_n_spots_h2ax_mean"]] + theme(legend.position = "none"),
    plot_list[["cc_cc_late_mitosis"]] + theme(legend.position = "none"),
    nrow = 1,
    align = "h"
)

sup_fig_part_b <- cowplot::plot_grid(
    bottom_row_gg,
    cell_line_legend,
    rel_widths = c(3, 0.4),
    ncol = 2,
    align = "v"
)

sup_fig_full <- cowplot::plot_grid(
    sup_fig_part_a,
    sup_fig_part_b,
    rel_heights = c(4, 2),
    nrow = 2,
    labels = c("a", "b"),
    align = "v"
)

sup_fig_full

 # Save figure
supfig_file <- file.path(figure_dir, "supplementary_figure_example_distributions.png")

cowplot::save_plot(
    filename = supfig_file,
    plot = sup_fig_full,
    base_height = 8,
    base_width = 8,
    dpi = 600
)

label_thresh_value = 0.925

pdf_file <- file.path(
    figure_dir,
    paste0("all_regression_performance_metrics_", consensus, ".pdf")
)
pdf(pdf_file, width = 6, height = 8, onefile = TRUE)

for (target in unique(y_plot_df$target)) {
    # Subset all dataframes
    y_subset_df <- y_plot_df %>% dplyr::filter(target == !!target)
    
    y_subset_df$shuffle <- y_subset_df$shuffle %>%
        dplyr::recode_factor("shuffle_false" = "Real",
                             "shuffle_true" = "Permuted")
    
    coef_subset_df <- full_coef_df %>%
        dplyr::filter(target == !!target)
    metrics_subset_df <- regression_metrics_df %>%
        dplyr::filter(target == !!target)
    
    for (y_transform in unique(y_subset_df$y_transform)) {
        y_subset_transform_df <- y_subset_df %>%
            dplyr::filter(y_transform == !!y_transform)
        
        coef_subset_transform_df <- coef_subset_df %>%
            dplyr::filter(y_transform == !!y_transform,
                          shuffle == "shuffle_false") %>%
            dplyr::mutate(weight_rank = row_number(weight))
        
        metrics_subset_transform_df <- metrics_subset_df %>%
            dplyr::filter(y_transform == !!y_transform)
        
        pred_scatter_gg <-
           ggplot(y_subset_transform_df,
                  aes(x = recode_target_value_true,
                      y = recode_target_value_pred)) +
                geom_point(aes(color = Metadata_cell_line),
                           size = 0.5, alpha = 0.8) +
                facet_grid(data_type~shuffle) +
                theme_bw() +
                xlab("True Values") +
                ylab("Predicted Values") +
                scale_color_manual(
                    name = "Cell Line",
                    labels = cell_line_labels,
                    values = cell_line_colors
                ) +
                geom_smooth(method='lm', formula=y~x) +
                theme(strip.text = element_text(size = 10),
                      strip.background = element_rect(colour = "black",
                                                      fill = "#fdfff4"))
        
        # Setup labeling thresholds
        non_zero_coef <- coef_subset_transform_df$abs_weight[coef_subset_transform_df$abs_weight > 0]
        label_thresh <- quantile(non_zero_coef, label_thresh_value)
        label_logic <- coef_subset_transform_df$abs_weight > label_thresh

        coef_gg <-
            ggplot(coef_subset_transform_df,
                   aes(x = weight_rank,
                       y = weight)) +
                geom_point(size = 0.2,
                           alpha = 0.6) +
                xlab("Weight Rank") +
                ylab("Weight") +
                geom_text_repel(
                    data = subset(coef_subset_transform_df, label_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    box.padding = 0.4,
                    point.padding = 0.1,
                    segment.size = 0.5,
                    segment.alpha = 0.6,
                    size = 1.5,
                    fontface = "italic",
                    aes(label = feature,
                        x = weight_rank,
                        y = weight)) +
                theme_bw()
        
         # Build table for plotting performance metrics
        mse_df <- metrics_subset_transform_df %>%
            dplyr::filter(metric == "mse") %>%
            dplyr::select(-metric)
        mse_df$value = abs(mse_df$value)

        r2_df <- metrics_subset_transform_df %>%
            dplyr::filter(metric == "r_two") %>%
            dplyr::rename(r2 = value) %>%
            dplyr::select(-metric)

        metric_table_df <- r2_df %>%
            dplyr::inner_join(mse_df,
                              by=c("target", "data_type", "shuffle", "y_transform")) %>%
            dplyr::select(data_type, shuffle, y_transform, r2, value) %>%
            dplyr::rename(fit = data_type, transform = y_transform) %>%
            dplyr::arrange(shuffle)

        metric_table_df$shuffle <- dplyr::recode(
            metric_table_df$shuffle,
            shuffle_true = "True",
            shuffle_false = "False"
        )

        # Plot all performance metrics together with cowplot
        table_theme <- gridExtra::ttheme_default(
            core = list(fg_params=list(cex = 0.7)),
            colhead = list(fg_params=list(cex = 0.8))
        )

        table_gg <- gridExtra::tableGrob(metric_table_df,
                                         theme = table_theme,
                                         rows = NULL)
        
        bottom_row_gg <- cowplot::plot_grid(
            table_gg,
            coef_gg,
            rel_widths = c(0.8, 1),
            nrow = 1
        )
        regression_perf_gg <- cowplot::plot_grid(
            pred_scatter_gg,
            bottom_row_gg,
            rel_heights = c(1, 0.8),
            nrow = 2
        )
        
        target_title <- cowplot::ggdraw() + 
          cowplot::draw_label(
            paste("Performance:", target, "\nTransform:", y_transform),
            fontface = 'bold',
            x = 0,
            hjust = 0
          ) +
          theme(
            plot.margin = margin(0, 0, 0, 7)
          )

        regression_perf_gg <- cowplot::plot_grid(
            target_title,
            regression_perf_gg,
            ncol = 1,
            rel_heights = c(0.1, 1)
        )
        
        # Save figure
        cowplot_file <- file.path(
            individual_fig_dir,
            paste0(target, "_", y_transform, "_performance_", consensus, ".png")
        )

        cowplot::save_plot(
            filename = cowplot_file,
            plot = regression_perf_gg,
            base_height = 6,
            base_width = 6
        )
        
        print(regression_perf_gg)
    }
}

dev.off()

 # Subset all dataframes
y_subset_df <- y_plot_df %>% dplyr::filter(target == !!target)

y_subset_df$shuffle <- y_subset_df$shuffle %>%
    dplyr::recode_factor("shuffle_false" = "Real",
                         "shuffle_true" = "Permuted")

coef_subset_df <- full_coef_df %>%
    dplyr::filter(target == !!target)
metrics_subset_df <- regression_metrics_df %>%
    dplyr::filter(target == !!target)


y_subset_transform_df <- y_subset_df %>%
    dplyr::filter(y_transform == !!y_transform)

coef_subset_transform_df <- coef_subset_df %>%
    dplyr::filter(y_transform == !!y_transform,
                  shuffle == "shuffle_false") %>%
    dplyr::mutate(weight_rank = row_number(weight))

metrics_subset_transform_df <- metrics_subset_df %>%
    dplyr::filter(y_transform == !!y_transform)

pred_scatter_gg <-
   ggplot(y_subset_transform_df,
          aes(x = recode_target_value_true,
              y = recode_target_value_pred)) +
        geom_point(aes(color = Metadata_cell_line),
                   size = 0.5, alpha = 0.8) +
        facet_grid(data_type~shuffle) +
        theme_bw() +
        xlab("True Values") +
        ylab("Predicted Values") +
        scale_color_manual(
            name = "Cell Line",
            labels = cell_line_labels,
            values = cell_line_colors
        ) +
        geom_smooth(method='lm', formula=y~x) +
        theme(strip.text = element_text(size = 10),
              strip.background = element_rect(colour = "black",
                                              fill = "#fdfff4"))

# Setup labeling thresholds
non_zero_coef <- coef_subset_transform_df$abs_weight[coef_subset_transform_df$abs_weight > 0]
label_thresh <- quantile(non_zero_coef, label_thresh_value)
label_logic <- coef_subset_transform_df$abs_weight > label_thresh

coef_gg <-
    ggplot(coef_subset_transform_df,
           aes(x = weight_rank,
               y = weight)) +
        geom_point(size = 0.2,
                   alpha = 0.6) +
        xlab("Weight Rank") +
        ylab("Weight") +
        geom_text_repel(
            data = subset(coef_subset_transform_df, label_logic),
            arrow = arrow(length = unit(0.01, "npc")),
            box.padding = 0.4,
            point.padding = 0.1,
            segment.size = 0.5,
            segment.alpha = 0.6,
            size = 1.5,
            fontface = "italic",
            aes(label = feature,
                x = weight_rank,
                y = weight)) +
        theme_bw()

 # Build table for plotting performance metrics
mse_df <- metrics_subset_transform_df %>%
    dplyr::filter(metric == "mse") %>%
    dplyr::select(-metric)
mse_df$value = abs(mse_df$value)

r2_df <- metrics_subset_transform_df %>%
    dplyr::filter(metric == "r_two") %>%
    dplyr::rename(r2 = value) %>%
    dplyr::select(-metric)

metric_table_df <- r2_df %>%
    dplyr::inner_join(mse_df,
                      by=c("target", "data_type", "shuffle", "y_transform")) %>%
    dplyr::select(data_type, shuffle, y_transform, r2, value) %>%
    dplyr::rename(fit = data_type, transform = y_transform) %>%
    dplyr::arrange(shuffle)

metric_table_df$shuffle <- dplyr::recode(
    metric_table_df$shuffle,
    shuffle_true = "True",
    shuffle_false = "False"
)

# Plot all performance metrics together with cowplot
table_theme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.7)),
    colhead = list(fg_params=list(cex = 0.8))
)

table_gg <- gridExtra::tableGrob(metric_table_df,
                                 theme = table_theme,
                                 rows = NULL)

bottom_row_gg <- cowplot::plot_grid(
    table_gg,
    coef_gg,
    rel_widths = c(0.8, 1),
    nrow = 1
)
regression_perf_gg <- cowplot::plot_grid(
    pred_scatter_gg,
    bottom_row_gg,
    rel_heights = c(1, 0.8),
    nrow = 2
)

target_title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    paste("Performance:", target, "\nTransform:", y_transform),
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )

regression_perf_gg <- cowplot::plot_grid(
    target_title,
    regression_perf_gg,
    ncol = 1,
    rel_heights = c(0.1, 1)
)

# Save figure
cowplot_file <- file.path(
    individual_fig_dir,
    paste0(target, "_", y_transform, "_performance_", consensus, ".png")
)
