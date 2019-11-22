suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

set.seed(123)

cell_lines <- c("A549", "ES2", "HCC44")

# Annotated Cell Health Features
feat_file <- file.path("..", "1.generate-profiles", "data", "labels", "feature_mapping_annotated.csv")
label_df <- readr::read_csv(feat_file, col_types = readr::cols())

head(label_df)

regression_file <- file.path("results", "full_cell_health_regression_results.tsv.gz")
regression_metrics_df <- readr::read_tsv(regression_file, col_types = readr::cols()) %>%
    dplyr::filter(cell_line %in% cell_lines)
    
head(regression_metrics_df)

ggplot(regression_metrics_df %>% dplyr::filter(metric == "mse"),
       aes(x = shuffle,
           y = value)) +
    geom_jitter(width = 0.2, size = 0.7, alpha = 0.8, pch = 16) +
    facet_grid(cell_line~data_fit) +
    xlab("Data Type") +
    ylab("Within Cell Line MSE\nof Cell Health Feature") +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

output_file <- file.path("figures", "cell_line_mse_differences.png")
ggsave(output_file, height = 5, width = 5, dpi = 500)

ggplot(regression_metrics_df %>% dplyr::filter(metric == "r_two"),
       aes(x = shuffle,
           y = value)) +
    geom_jitter(width = 0.2, size = 0.7, alpha = 0.8, pch = 16) +
    facet_wrap(cell_line~data_fit, ncol = 2, scales = "free_x") +
    xlab("Data Type") +
    ylab("Within Cell Line R Squared\nof Cell Health Feature") +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

output_file <- file.path("figures", "cell_line_rsquared_differences.png")
ggsave(output_file, height = 5, width = 5, dpi = 500)

# Compile Results
results_df <- regression_metrics_df %>%
    dplyr::left_join(label_df, by = c("target" = "updated_name")) %>%
    dplyr::mutate(plot_group = paste(metric, target, shuffle))

dim(results_df)
head(results_df, 2)

# Set some plotting defaults
measurement_colors <- c(
    "apoptosis" = "#a6cee3",
    "cell_cycle_arrest" = "#1f78b4",
    "cell_viability" = "#b2df8a",
    "death" = "#33a02c",
    "dna_damage" = "#fb9a99", 
    "g1_arrest" = "#fdbf6f",
    "g2_arrest" = "#ff7f00",
    "g2_m_arrest" = "#005c8c",
    "mitosis" = "green",
    "other" = "black",
    "s_arrest" = "#cab2d6",
    "toxicity" = "#6a3d9a"
)

measurement_labels <- c(
    "apoptosis" = "Apoptosis",
    "cell_cycle_arrest" = "Cell Cycle Arrest",
    "cell_viability" = "Cell Viability",
    "death" = "Death",
    "dna_damage" = "DNA Damage", 
    "g1_arrest" = "G1 Arrest",
    "g2_arrest" = "G2 Arrest",
    "g2_m_arrest" = "G2/M Arrest",
    "mitosis" = "Mitosis",
    "other" = "Other",
    "s_arrest" = "S Arrest",
    "toxicity" = "Toxicity"
)

dye_colors <- c(
    "hoechst" = "#639B94",
    "edu" = "#E45242",
    "gh2ax" = "#E2C552",
    "ph3" = "#7B9C32",
    "hoechst_gh2ax" = "#535f52",
    "hoechst_edu" = "#73414b",
    "edu_gh2ax" = "#e37a48",
    "caspase" = "#F7B1C1",
    "draq" = "#FF6699",
    "draq_caspase" = "#7f4a72",
    "many_cell_cycle" = "#E9DFC3",
    "crispr_efficiency" = "black"
)

dye_labels <- c(
    "hoechst" = "Hoechst",
    "edu" = "EdU",
    "gh2ax" = "gH2AX",
    "ph3" = "pH3",
    "hoechst_gh2ax" = "Hoechst + gH2AX",
    "hoechst_edu" = "Hoechst + EdU",
    "edu_gh2ax" = "EdU + gH2AX",
    "caspase" = "Caspase 3/7",
    "draq" = "DRAQ7",
    "draq_caspase" = "DRAQ7 + Caspase 3/7",
    "many_cell_cycle" = "Cell Cycle (Many Dyes)",
    "crispr_efficiency" = "CRISPR Efficiency"
)

head(results_df, 2)

ggplot(results_df %>%
       dplyr::filter(data_fit == "test"),
       aes(x = cell_line,
           y = value,
           group = plot_group)) +
    geom_jitter(aes(color = assay), width = 0.01) +
    geom_line(aes(color = assay),
              alpha = 0.5) +
    scale_color_manual(name = "Measurement",
                       values = dye_colors,
                       labels = dye_labels) +
    facet_wrap(metric~shuffle, scales = "free") +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

output_file <- file.path("figures", "cell_line_differences_target_linked_full.png")
ggsave(output_file, height = 5, width = 5, dpi = 500)

ggplot(results_df %>%
       dplyr::filter(data_fit == "test",
                     value > -1),
       aes(x = cell_line,
           y = value,
           group = plot_group)) +
    geom_jitter(aes(color = assay), width = 0.01) +
    geom_line(aes(color = assay),
              alpha = 0.5) +
    ylab("value (-1 set as min value)") +
    scale_color_manual(name = "Measurement",
                       values = dye_colors,
                       labels = dye_labels) +
    facet_wrap(metric~shuffle, scales = "free") +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

output_file <- file.path("figures", "cell_line_differences_target_linked_subset.png")
ggsave(output_file, height = 5, width = 5, dpi = 500)

filtered_results_df <- results_df %>%
    dplyr::filter(cell_line == "A549",
                  data_fit == "test",
                  metric == "r_two") %>%
    dplyr::select(value, metric, shuffle, target, original_name,
                  feature_type, measurement, assay, description) %>%
    tidyr::spread(shuffle, value) %>%
    dplyr::arrange(desc(shuffle_false))

filtered_results_df$target <- factor(filtered_results_df$target,
                                     levels = rev(unique(filtered_results_df$target)))
filtered_results_df$original_name <- factor(filtered_results_df$original_name,
                                            levels = rev(unique(filtered_results_df$original_name)))


# Output ranked models
output_file <- file.path("..", "4.apply", "repurposing_cellhealth_shiny", "data", "A549_ranked_models.tsv")
readr::write_tsv(filtered_results_df, output_file)

print(dim(filtered_results_df))
head(filtered_results_df, 10)

ggplot(filtered_results_df, aes(x = shuffle_true, y = shuffle_false)) +
    geom_point(aes(color = assay),
               size = 0.5,
               alpha = 0.8) +
    xlab("Random Shuffle Regression Performance (Test Set)") +
    ylab("Real Values Regression Performance (Test Set)") +
    ggtitle("A549 Cell Line") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 7, angle = 90),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"))

output_file = file.path("figures", "ranked_models_A549_with_shuffle.png")
ggsave(output_file, dpi = 300, height = 3.5, width = 4)

ggplot(filtered_results_df %>% dplyr::filter(shuffle_false > 0),
       aes(x = original_name, y = shuffle_false)) +
    geom_bar(aes(fill = assay), stat="identity") +
    ylab("Test Set Regression Performance") +
    xlab("") +
    ggtitle("A549 Cell Line") +
    coord_flip() +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 7, angle = 90),
          axis.title = element_text(size = 8),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.3, "cm"))

output_file = file.path("figures", "ranked_models_A549.png")
ggsave(output_file, dpi = 300, height = 6, width = 6)
