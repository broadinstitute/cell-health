suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source(file.path("..", "3.train", "scripts", "assay_themes.R"))

# Set seed
set.seed(123)

# Set constants
consensus <- "modz"
data_dir <- "results"
figure_dir <- "figures"
original_results_dir <- file.path("..", "3.train", "results")

# Input file names
cell_line_file <- file.path(
    data_dir, paste0("cell_line_holdout_robustness_results_", consensus, ".tsv.gz")
)
sample_titration_file <- file.path(
    data_dir, paste0("sample_titration_robustness_results_", consensus, ".tsv.gz")
)
feature_group_file <- file.path(
    data_dir, paste0("feature_group_subset_removal_robustness_results_", consensus, ".tsv.gz")
)
original_results_file <- file.path(
    original_results_dir, paste0("full_cell_health_regression_", consensus, ".tsv.gz")
)
label_file <- file.path(
    "..",
    "1.generate-profiles",
    "data",
    "labels",
    "feature_mapping_annotated.csv"
)

# Output file names
cell_line_figure_file <- file.path(
    figure_dir, paste0("cell_line_holdout_", consensus, ".png")
)

sample_titration_figure_file <- file.path(
    figure_dir, paste0("sample_titration_", consensus, ".png")
)

feature_removal_figure_file <- file.path(
    figure_dir, paste0("feature_removal_", consensus, ".png")
)

# Load data
cell_line_df <- readr::read_tsv(cell_line_file, col_types = readr::cols())
sample_df <- readr::read_tsv(sample_titration_file, col_types = readr::cols())
feature_df <- readr::read_tsv(feature_group_file, col_types = readr::cols())
original_df <- readr::read_tsv(original_results_file, col_types = readr::cols())

original_df$shuffle <- original_df$shuffle %>%
    dplyr::recode_factor(
        "shuffle_false" = "Real",
        "shuffle_true" = "Permuted"
    )

# Label variables with specific cell health classes
label_df <- readr::read_csv(label_file, col_types = readr::cols())

cell_line_df <- cell_line_df %>%
    dplyr::filter(metric == "r_two", data_fit == "test") 

cell_line_df$shuffle <- cell_line_df$shuffle %>%
    dplyr::recode_factor(
        "shuffle_false" = "Real",
        "shuffle_true" = "Permuted"
    )

cell_line_df <- cell_line_df %>%
    dplyr::bind_rows(
        original_df %>%
            dplyr::filter(metric == "r_two",
                          y_transform == "raw",
                          cell_line == "all",
                          data_fit == "test")
    ) %>%
    dplyr::left_join(label_df, by = c("target" = "id")) %>%
    dplyr::mutate(truncated_value = ifelse(value < -1, -1, value))

cell_line_df$measurement <- dplyr::recode(cell_line_df$measurement, "metadata" = "other")

cell_line_df$measurement_ordered <- dplyr::recode_factor(
    cell_line_df$measurement, !!!measurement_labels
)
cell_line_df$measurement_ordered <- dplyr::recode_factor(
    cell_line_df$measurement_ordered,
    `Reactive Oxygen Species` = "ROS",
    `Infection Efficiency` = "Inf. Efficiency"
)

assay_measurement_order <- cell_line_df %>%
    dplyr::filter(shuffle == "Real", cell_line == "A549") %>%
    dplyr::group_by(measurement_ordered) %>%
    dplyr::mutate(median_measure = median(value)) %>%
    dplyr::select(measurement_ordered, median_measure) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(median_measure)) %>%
    dplyr::pull(measurement_ordered)

cell_line_df$measurement_ordered <- factor(
    cell_line_df$measurement_ordered, levels = assay_measurement_order
)

cell_line_df$cell_line <- cell_line_df$cell_line %>%
    dplyr::recode_factor(
        "A549" = "A549 (Holdout)",
        "ES2" = "ES2 (Holdout)",
        "HCC44" = "HCC44 (Holdout)",
        "all" = "Original test set",
    )

cell_line_gg <- ggplot(cell_line_df,
       aes(y = truncated_value, x = measurement_ordered)) +
    geom_boxplot(fill = "grey",
                 alpha = 0.72,
                 size = 0.5,
                 color = "black",
                 outlier.alpha = 0) +
     geom_jitter(fill = "black",
                 alpha = 0.72,
                 size = 1,
                 shape = 21,
                 width = 0.1,
                 color = "black") +
    scale_y_continuous(
        limits = c(-1, 1),
        breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
        labels = c("\u2264 -1.0", "-0.5", "0", "0.5", "1.0")
        ) +
    ylab(bquote("Test Set Model Performance ("~R^2~")")) +
    xlab("") +
    geom_hline(yintercept = 0,
               alpha = 0.5, 
               linetype = "dashed") +
    facet_grid(cell_line~shuffle,
               scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          legend.position = "none",
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

ggsave(cell_line_figure_file, cell_line_gg, dpi = 500, width = 8, height = 7)
cell_line_gg

# Append original results to sample titration
sample_df <- sample_df %>%
    dplyr::filter(
        metric == "r_two", shuffle == "shuffle_false", cell_line == "all", data_fit == "test"
    ) %>%
    dplyr::bind_rows(
        original_df %>%
            dplyr::filter(
                metric == "r_two", shuffle == "Real", cell_line == "all", data_fit == "test"
            ) %>%
            dplyr::mutate(num_samples_dropped = 0, iteration = 0)
        ) %>%
    dplyr::left_join(label_df, by = c("target" = "id")) %>%
    dplyr::group_by(target, num_samples_dropped) %>%
    dplyr::mutate(mean_value = mean(value)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(mean_value) %>%
    dplyr::distinct(target, num_samples_dropped, mean_value, readable_name, measurement)

target_order <- sample_df %>%
    dplyr::arrange(desc(mean_value)) %>%
    dplyr::distinct(readable_name) %>%
    dplyr::pull(readable_name)

sample_df$readable_name <- factor(sample_df$readable_name, levels=target_order)

# Recode metadata measurent to other
sample_df$measurement <- dplyr::recode(sample_df$measurement, "metadata" = "other")

sample_df$measurement <- factor(sample_df$measurement, levels = c(
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

print(dim(sample_df))
head(sample_df, 3)

avg_results_full_df <- sample_df %>%
    dplyr::group_by(num_samples_dropped) %>%
    dplyr::mutate(mean_total_estimate = mean(mean_value)) %>%
    dplyr::distinct(num_samples_dropped, mean_total_estimate)

head(avg_results_full_df)

sample_gg <- (
    ggplot(sample_df, aes(x = num_samples_dropped, y = mean_value)) +
        geom_line(aes(color = measurement, group = target)) +
        ylab(bquote("Test set "~R^2~"")) +
        xlab("Number of samples\ndropped from training") +
        scale_color_manual(
            name = "Measurement",
            values = measurement_colors,
            labels = measurement_labels
        ) +
        theme_bw()
    )

avg_drop_gg <- ggplot(
    avg_results_full_df,
    aes(x = num_samples_dropped,
        y = mean_total_estimate)
    ) +
    geom_point(
        data = ,
        color = "black",
        size = 1.5)  +
    geom_line(
        data = avg_results_full_df,
        color = "black",
        size = 0.75
    ) +
    theme_bw() +
    xlab("Number of samples\ndropped from training") +
    ylab(bquote("Test set "~R^2~"\n(Mean for all Cell Health targets)"))

sample_titration_gg <- cowplot::plot_grid(
    sample_gg,
    avg_drop_gg,
    ncol = 2,
    align = "h",
    axis = "l",
    rel_widths = c(1, 0.7),
    labels = c("a", "b")
)

cowplot::save_plot(sample_titration_figure_file, sample_titration_gg, base_width = 10, base_height = 6)
sample_titration_gg

feature_df <- feature_df %>%
    dplyr::filter(metric == "r_two", data_fit == "test") %>%
    dplyr::left_join(label_df, by = c("target" = "id")) %>%
    dplyr::mutate()

feature_df$measurement <- dplyr::recode(feature_df$measurement, "metadata" = "other")

feature_df$measurement_ordered <- dplyr::recode_factor(feature_df$measurement, !!!measurement_labels)
feature_df$measurement_ordered <- dplyr::recode_factor(
    feature_df$measurement_ordered, `Reactive Oxygen Species` = "ROS", `Infection Efficiency` = "Inf. Efficiency"
)

assay_measurement_order <- feature_df %>%
    dplyr::group_by(measurement_ordered) %>%
    dplyr::mutate(median_measure = median(value)) %>%
    dplyr::select(measurement_ordered, median_measure) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(median_measure)) %>%
    dplyr::pull(measurement_ordered)

feature_df$measurement_ordered <- factor(
    feature_df$measurement_ordered, levels = assay_measurement_order
)

feature_df <- feature_df %>% dplyr::left_join(
    original_df %>%
        dplyr::filter(
            shuffle == "Real",
            cell_line == "all",
            data_fit == "test",
            metric == "r_two"
        ) %>%
        dplyr::select(value, target),
    by = "target",
    suffix = c("", "_original")
    ) %>%
    dplyr::mutate(feature_group_difference = value_original - value,
                  feature_group_name = paste0(feature_type_dropped, " (n = ", num_dropped, ")"))

feature_order <- feature_df %>%
    dplyr::group_by(readable_name) %>%
    dplyr::mutate(median_measure = mean(feature_group_difference)) %>%
    dplyr::select(readable_name, median_measure) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::arrange(median_measure) %>%
    dplyr::pull(readable_name)

feature_df$readable_name <- factor(
    feature_df$readable_name, levels = feature_order
)

head(feature_df, 2)

feature_theme <- theme(
    axis.title = element_text(size = 7),
    axis.text.y = element_text(size = 6),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
)

panel_a_gg <- ggplot(
        feature_df %>% dplyr::filter(feature_category == "feature_group"),
        aes(x = feature_group_name, y = feature_group_difference)
    ) +
    coord_flip() +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() +
    feature_theme +
    ylim(c(
        min(feature_df$feature_group_difference) - 0.1,
        max(feature_df$feature_group_difference + 0.1)
    )) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab("Feature group\ndropped") +
    ylab("Performance difference\n(Original - dropped)")

panel_b_gg <- ggplot(
        feature_df %>% dplyr::filter(feature_category == "channel"),
        aes(x = feature_group_name, y = feature_group_difference)
    ) +
    coord_flip() +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() +
    feature_theme +
    ylim(c(
        min(feature_df$feature_group_difference) - 0.1,
        max(feature_df$feature_group_difference + 0.1)
    )) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab("Channel\ndropped") +
    ylab("Performance difference\n(Original - dropped)")

panel_c_gg <- ggplot(
        feature_df %>% dplyr::filter(feature_category == "compartment"),
        aes(x = feature_group_name, y = feature_group_difference)
    ) +
    coord_flip() +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() +
    feature_theme +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    ylim(c(
        min(feature_df$feature_group_difference) - 0.1,
        max(feature_df$feature_group_difference + 0.1)
    )) +
    xlab("Compartment\ndropped") +
    ylab("Performance difference\n(Original - dropped)")

panel_d_gg <- ggplot(feature_df,
                     aes(x = readable_name,
                         y = feature_group_difference,
                         fill = measurement)) +
    geom_boxplot(outlier.size = 0.5) +
    coord_flip() +
    scale_fill_manual(
        name = "Measurement",
        values = measurement_colors,
        labels = measurement_labels
    ) +
    theme_bw() +
    feature_theme +
    geom_hline(yintercept=0, color = "red", linetype = "dashed") +
    ylab("Performance difference\n(Original - dropped)") +
    xlab("")

left_panel <- cowplot::plot_grid(
    panel_a_gg,
    panel_b_gg,
    panel_c_gg,
    nrow = 3,
    align = "hv",
    axis = "l",
    rel_heights = c(1, 0.75, 0.55),
    labels = c("a", "b", "c")
)

feature_removal_gg <- cowplot::plot_grid(
    left_panel,
    panel_d_gg,
    ncol = 2,
    align = "v",
    axis = "b",
    rel_widths = c(0.5, 1),
    labels = c("", "d")
)

cowplot::save_plot(
    feature_removal_figure_file, feature_removal_gg, base_width = 9, base_height = 7
)
feature_removal_gg
