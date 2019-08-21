suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpmisc))

results_dir <- "results"

# Regression Results
regression_file <- file.path(results_dir, "full_cell_health_regression_results.tsv.gz")
regression_metrics_df <- readr::read_tsv(regression_file, col_types = readr::cols())

# Classification Results
roc_file <- file.path(results_dir, "full_cell_health_roc_results.tsv.gz")
full_roc_df <- readr::read_tsv(roc_file, col_types = readr::cols())

pr_file <- file.path(results_dir, "full_cell_health_pr_results.tsv.gz")
full_pr_df <- readr::read_tsv(pr_file, col_types = readr::cols())

# Model Coefficients
coef_file <- file.path(results_dir, "full_cell_health_coefficients.tsv.gz")
full_coef_df <- readr::read_tsv(coef_file, col_types = readr::cols()) %>%
    dplyr::filter(y_transform %in% c("raw", "zero-one"))

# Model Predictions
y_file <- file.path(results_dir, "full_cell_health_y_labels.tsv.gz")
y_df <- readr::read_tsv(y_file, col_types = readr::cols()) %>%
    dplyr::filter(y_transform %in% c("raw", "zero-one"))

# Annotated Cell Health Features
feat_file <- file.path("..", "0.generate-profiles", "data", "labels", "feature_mapping_annotated.csv")
label_df <- readr::read_csv(feat_file, col_types = readr::cols())

regression_subset_df <- regression_metrics_df %>%
    dplyr::filter(y_transform == "raw",
                  data_fit == "test",
                  shuffle == "shuffle_false") %>%
    tidyr::spread(key = "metric", value = "mse") %>%
    dplyr::select(-y_transform)

print(dim(regression_subset_df))
head(regression_subset_df, 3)

auroc_df <- full_roc_df %>%
    dplyr::distinct(metric, target, auc, data_fit,
                    shuffle, y_transform, min_class_count)

aupr_df <- full_pr_df %>%
    dplyr::distinct(metric, target, auc, data_fit,
                    shuffle, y_transform, min_class_count)

auc_df <- dplyr::bind_rows(auroc_df, aupr_df) %>%
    dplyr::filter(shuffle == "shuffle_false",
                  data_fit == "test")

auc_df$metric <- dplyr::recode_factor(
    auc_df$metric,
    "roc" = "AUROC",
    "aupr" = "AUPR"
)

auc_df <- auc_df %>%
    tidyr::spread(key = "metric", value = "auc") %>%
    dplyr::select(-y_transform)

print(dim(auc_df))
head(auc_df, 2)

metric_df <- regression_subset_df %>%
    dplyr::inner_join(auc_df, by=c("target", "data_fit", "shuffle")) %>%
    dplyr::left_join(label_df, by=c("target" = "updated_name"))

metric_df$mse = abs(metric_df$mse)

metric_df$maria_thumbs_up <- tidyr::replace_na(metric_df$maria_thumbs_up, 0)
metric_df$measurement <- tidyr::replace_na(metric_df$measurement, "other")

metric_df$measurement <- factor(metric_df$measurement,
                                levels = c("apoptosis", "death", "cell_viability", "toxicity", "dna_damage",
                                           "cell_cycle_arrest", "g1_arrest", "g2_arrest", "g2_m_arrest",
                                           "s_arrest", "other"))
print(dim(metric_df))
head(metric_df, 3)

ggplot(metric_df, aes(x = AUROC, y = mse)) +
    geom_point(alpha = 0.7,
               aes(color = measurement,
                   shape = as.factor(maria_thumbs_up))) +
    xlab("Classification\nAUROC (Test Set)") +
    ylab("Regression\nMean Squared Error (Test Set)") +
    geom_vline(xintercept = 0.5,
               alpha = 0.5, 
               linetype = "dashed") +
    scale_shape_manual(name = "Maria :+1:",
                       values = c(16, 3)) +
    scale_color_manual(name = "Measurement",
                       values = c("apoptosis" = "#ff0000",
                                  "cell_cycle_arrest" = "#ff5900",
                                  "cell_viability" = "#ff9900",
                                  "death" = "#ffcc00",
                                  "dna_damage" = "#E8B360", 
                                  "g1_arrest" = "#00c20d",
                                  "g2_arrest" = "#008f4c",
                                  "g2_m_arrest" = "#005c8c",
                                  "other" = "black",
                                  "s_arrest" = "#0033bf",
                                  "toxicity" = "#feab00"),
                       labels = c("apoptosis" = "Apoptosis",
                                  "cell_cycle_arrest" = "Cell Cycle Arrest",
                                  "cell_viability" = "Cell Viability",
                                  "death" = "Death",
                                  "dna_damage" = "DNA Damage", 
                                  "g1_arrest" = "G1 Arrest",
                                  "g2_arrest" = "G2 Arrest",
                                  "g2_m_arrest" = "G2/M Arrest",
                                  "other" = "Other",
                                  "s_arrest" = "S Arrest",
                                  "toxicity" = "Toxicity")) +
    theme_bw()

file <- file.path("figures", "performance_summary.png")
ggsave(file, dpi = 300, width = 6, height = 4.5)
