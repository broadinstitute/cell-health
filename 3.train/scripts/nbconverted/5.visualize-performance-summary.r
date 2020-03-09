suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))

source(file.path("scripts", "assay_themes.R"))

consensus <- "modz"

results_dir <- "results"
figure_dir <- file.path("figures", "summary", consensus)
cytominer_compare_dir <- file.path("figures", "cytominer_comparison", consensus)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cytominer_compare_dir, recursive = TRUE, showWarnings = FALSE)

# Regression Results
regression_file <- file.path(
    results_dir, 
    paste0("full_cell_health_regression_", consensus, ".tsv.gz")
)
regression_metrics_df <- readr::read_tsv(regression_file, col_types = readr::cols()) %>%
    dplyr::filter(cell_line == "all")

# Classification Results
roc_file <- file.path(
    results_dir,
    paste0("full_cell_health_roc_results_", consensus, ".tsv.gz")
)
full_roc_df <- readr::read_tsv(roc_file, col_types = readr::cols()) %>%
    dplyr::filter(cell_line == "all")

pr_file <- file.path(
    results_dir,
    paste0("full_cell_health_pr_results_", consensus, ".tsv.gz")
)
full_pr_df <- readr::read_tsv(pr_file, col_types = readr::cols()) %>%
    dplyr::filter(cell_line == "all")

# Model Coefficients
coef_file <- file.path(
    results_dir,
    paste0("full_cell_health_coefficients_", consensus, ".tsv.gz")
)
full_coef_df <- readr::read_tsv(coef_file, col_types = readr::cols())

# Model Predictions
y_file <- file.path(
    results_dir,
    paste0("full_cell_health_y_labels_", consensus, ".tsv.gz")
)
y_df <- readr::read_tsv(y_file, col_types = readr::cols()) %>%
    dplyr::filter(y_transform == "raw")

# Annotated Cell Health Features
feat_file <- file.path(
    "..",
    "1.generate-profiles",
    "data",
    "labels",
    "feature_mapping_annotated.csv"
)
label_df <- readr::read_csv(feat_file, col_types = readr::cols())

regression_subset_df <- regression_metrics_df %>%
    dplyr::filter(y_transform == "raw",
                  data_fit == "test",
                  shuffle == "shuffle_false") %>%
    tidyr::spread(key = "metric", value = "value") %>%
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
    dplyr::filter(shuffle == "shuffle_false")

auc_df$metric <- dplyr::recode_factor(
    auc_df$metric,
    "roc" = "AUROC",
    "aupr" = "AUPR"
)

auc_df <- auc_df %>%
    tidyr::spread(key = "metric", value = "auc") %>%
    dplyr::select(-y_transform)

test_auc_df <- auc_df %>%
    dplyr::filter(data_fit == "test")

train_auc_df <- auc_df %>%
    dplyr::filter(data_fit == "train") %>%
    dplyr::select(-data_fit)

auc_df <- train_auc_df %>%
    dplyr::full_join(test_auc_df,
                     by = c("target", "shuffle", "min_class_count"),
                     suffix = c("_train", "_test"))

print(dim(auc_df))
head(auc_df, 2)

metric_df <- regression_subset_df %>%
    dplyr::inner_join(auc_df, by = c("target", "data_fit", "shuffle")) %>%
    dplyr::left_join(label_df, by = c("target" = "id"))

metric_df$mse = abs(metric_df$mse)

metric_df$measurement <- tidyr::replace_na(metric_df$measurement, "other")

metric_df$measurement <- factor(
    metric_df$measurement,
    levels = measurement_levels
)

metric_df$assay <- factor(
    metric_df$assay,
    levels = dye_levels
)

print(dim(metric_df))
head(metric_df, 3)

# Note that the points that failed to plot did not contain enough positive samples in test set
ggplot(metric_df,
       aes(x = AUROC_test,
           y = AUROC_train)) +
    geom_point(alpha = 0.6,
               size = 1.7,
               aes(color = assay)) +
    ggtitle("Classification Summary") +
    xlab("Test Set - AUROC") +
    ylab("Training Set - AUROC") +
    geom_vline(xintercept = 0.5,
               alpha = 0.5, 
               linetype = "dashed") +
    geom_hline(yintercept = 0.5,
               alpha = 0.5, 
               linetype = "dashed") +
    geom_abline(slope = 1,
                intercept = 0, 
                color = "red",
                linetype = "dashed") +
    coord_fixed() +
    scale_color_manual(name = "Assay",
                       values = dye_colors,
                       labels = dye_labels) +
    theme_bw()

file <- file.path(
    figure_dir,
    paste0("performance_summary_auroc_assay_", consensus, ".png")
)
ggsave(file, dpi = 300, width = 6, height = 4.5)

ggplot(metric_df, aes(x = AUROC_test,
                      y = r_two)) +
    geom_point(alpha = 0.6,
               size = 1.7,
               aes(color = assay)) +
    xlab("Test Set - AUROC") +
    ylab("Test Set - R Squared") +
    geom_vline(xintercept = 0.5,
               alpha = 0.5, 
               linetype = "dashed") +
    geom_hline(yintercept = 0,
               alpha = 0.5, 
               linetype = "dashed") +
    coord_fixed() +
    scale_color_manual(name = "Assay",
                       values = dye_colors,
                       labels = dye_labels) +
    theme_bw()

file <- file.path(
    figure_dir,
    paste0(
        "performance_summary_assay_classification_vs_regression_",
        consensus,
        ".png"
    )
)
ggsave(file, dpi = 300, width = 6, height = 4.5)

# Get data ready for plotting
auc_test_full_df <- dplyr::bind_rows(auroc_df, aupr_df) %>%
    dplyr::left_join(label_df, by = c("target" = "id")) %>%
    dplyr::filter(data_fit == "test")

auc_test_full_real_df <- auc_test_full_df %>%
    dplyr::filter(shuffle == "shuffle_false") %>%
    dplyr::arrange(target, metric)
auc_test_full_shuffle_df <- auc_test_full_df %>%
    dplyr::filter(shuffle == "shuffle_true") %>%
    dplyr::arrange(target, metric)

real_minus_shuffle_df <- auc_test_full_real_df %>%
    dplyr::mutate(
        auc = auc_test_full_real_df$auc - auc_test_full_shuffle_df$auc,
        shuffle = "real_minus_shuffle"
    )

target_order <- real_minus_shuffle_df %>%
    dplyr::filter(metric == "roc") %>%
    dplyr::arrange(desc(auc)) %>%
    dplyr::pull(target)

auc_test_full_df <- auc_test_full_df %>%
    dplyr::bind_rows(real_minus_shuffle_df)

auc_test_full_df$metric <- auc_test_full_df$metric %>%
    dplyr::recode("aupr" = "AUPR", "roc" = "AUROC")

auc_test_full_df$target <- factor(auc_test_full_df$target, levels = rev(target_order))
auc_test_full_df$metric <- factor(auc_test_full_df$metric, levels = c("AUROC", "AUPR"))

print(dim(auc_test_full_df))
head(auc_test_full_df, 2)

# Get side by side plot
metric_df$target <- factor(metric_df$target, levels = rev(target_order))
metric_df <- metric_df %>%
    dplyr::mutate(r_two_rank = as.numeric(paste(rank(metric_df$r_two, ties.method = "first"))))

head(metric_df, 2)

r_two_df <- regression_metrics_df %>%
    dplyr::filter(metric == "r_two",
                  shuffle == "shuffle_false",
                  y_transform == "raw") %>%
    tidyr::spread(data_fit, value) %>%
    dplyr::left_join(label_df, by=c("target" = "id"))

r_two_df$measurement <- factor(
    r_two_df$measurement,
    levels = measurement_levels
)
r_two_df$measurement <- tidyr::replace_na(r_two_df$measurement, "other")

r_two_df$assay <- factor(
    r_two_df$assay,
    levels = dye_levels
)

print(dim(r_two_df))
head(r_two_df, 2)

ggplot(r_two_df,
       aes(y = train, x = test)) +
    geom_point(alpha = 0.6,
               size = 2,
               aes(color = measurement)) +
    geom_vline(xintercept = 0,
               alpha = 0.5, 
               linetype = "dashed") +
    geom_hline(yintercept = 0,
               alpha = 0.5,
               linetype = "dashed") +
    coord_fixed() +
    ylab("Training R squared") +
    xlab("Testing R squared") +
    geom_abline(intercept = 0,
            slope = 1,
            linetype = "dashed",
            color = "red",
            alpha = 0.7) +
    scale_color_manual(name = "Measurement",
                       values = measurement_colors,
                       labels = measurement_labels) +
    theme_bw()

file <- file.path(
    figure_dir,
    paste0("performance_summary_rsquared_", consensus, ".png")
)
ggsave(file, dpi = 300, width = 6, height = 4.25)

ggplot(r_two_df,
       aes(y = train, x = test)) +
    geom_point(alpha = 0.6,
               size = 1.7,
               aes(color = assay)) +
    ylab("Training -  R squared") +
    xlab("Test Set - R squared") +
    geom_vline(xintercept = 0,
               alpha = 0.5, 
               linetype = "dashed") +
    geom_hline(yintercept = 0,
               alpha = 0.5, 
               linetype = "dashed") +
    geom_abline(slope = 1,
                intercept = 0, 
                color = "red",
                linetype = "dashed") +
    coord_fixed() +
    scale_color_manual(name = "Assay",
                       values = dye_colors,
                       labels = dye_labels) +
    theme_bw()

file <- file.path(
    figure_dir,
    paste0("performance_summary_rsquared_assay_", consensus, ".png")
)
ggsave(file, dpi = 300, width = 6, height = 4.25)

# Cytominer results are archived on github
hash = "26d1095c209d402102494c0c28e978476643e572"

cyto_file = paste0(
    "https://github.com/broadinstitute/cell-health/raw/",
    hash,
    "/3.train/results/full_cell_health_regression_results.tsv.gz"
)

cyto_regression_df = readr::read_tsv(cyto_file, col_types=readr::cols()) %>%
    dplyr::mutate(package = "cytominer")

head(cyto_regression_df, 2)

regression_metrics_df <- regression_metrics_df %>%
    dplyr::mutate(package = "pycytominer")

head(regression_metrics_df, 2)

# Note that mse_diff is coded as "value" in pycytominer and "mse" in cytominer
all_regression_df <- regression_metrics_df %>%
    dplyr::inner_join(
        cyto_regression_df,
        by = c("metric", "target", "data_fit", "shuffle", "y_transform"),
        suffix = c("_pycytominer", "_cytominer")
    ) %>%
    dplyr::mutate(mse_diff = value - mse)

head(all_regression_df, 2)

for (metric in unique(all_regression_df$metric)) {
    all_regression_subset_df <- all_regression_df %>%
        dplyr::filter(metric == metric)
    
    if (metric == "mse") {
        all_regression_subset_df$value <- abs(all_regression_subset_df$value)
        all_regression_subset_df$mse <- abs(all_regression_subset_df$mse)
    }
    
    mse_gg <- ggplot(all_regression_subset_df,
       aes(x=value, 
           y=mse)) +
    geom_point(aes(color = data_fit),
               alpha = 0.7,
               size = 0.8) +
    facet_grid(y_transform~shuffle) +
    xlab("pycytominer") +
    ylab("cytominer") +
    ggtitle(metric) +
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dashed",
                color = "black",
                alpha = 0.7) +
    theme_bw() +
    theme(strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))
    
    print(mse_gg)
    outfile <- file.path(
        cytominer_compare_dir,
        paste0("compare_pycytominer_cytominer_", metric, "_", consensus, ".png")
    )
    ggsave(outfile, height = 5, width = 6, dpi = 300)
}



all_regression_df %>%
    dplyr::group_by(metric, data_fit, shuffle, y_transform) %>%
    dplyr::mutate(percent_pycytominer_better = sum(mse_diff > 0) / n()) %>%
    dplyr::distinct(metric, data_fit, shuffle, y_transform, percent_pycytominer_better)

# Cytominer results are archived on github
cyto_file = paste0(
    "https://github.com/broadinstitute/cell-health/raw/",
    hash,
    "/3.train/results/full_cell_health_roc_results.tsv.gz"
)

cyto_roc_df = readr::read_tsv(cyto_file, col_types=readr::cols())

head(cyto_roc_df, 2)

cyto_file = paste0(
    "https://github.com/broadinstitute/cell-health/raw/",
    hash,
    "/3.train/results/full_cell_health_pr_results.tsv.gz"
)

cyto_pr_df = readr::read_tsv(cyto_file, col_types=readr::cols())

head(cyto_pr_df, 2)

cyto_auroc_df <- cyto_roc_df %>%
    dplyr::distinct(metric, target, auc, data_fit,
                    shuffle, y_transform, min_class_count)

cyto_aupr_df <- cyto_pr_df %>%
    dplyr::distinct(metric, target, auc, data_fit,
                    shuffle, y_transform, min_class_count)

cyto_auc_df <- dplyr::bind_rows(cyto_auroc_df, cyto_aupr_df) %>%
    dplyr::filter(shuffle == "shuffle_false",
                  data_fit == "test")

cyto_auc_df$metric <- dplyr::recode_factor(
    cyto_auc_df$metric,
    "roc" = "AUROC",
    "aupr" = "AUPR"
)

cyto_auc_df <- cyto_auc_df %>%
    tidyr::spread(key = "metric", value = "auc") %>%
    dplyr::select(-y_transform) %>%
    dplyr::mutate(package = "cytominer")

print(dim(cyto_auc_df))
head(cyto_auc_df, 2)

auc_df <- auc_df %>%
    dplyr::mutate(package = "pycytominer")

all_classification_df <- auc_df %>%
    dplyr::inner_join(cyto_auc_df,
                      by = c("target", "data_fit", "shuffle"),
                      suffix = c("_pycytominer", "_cytominer")) %>%
    dplyr::mutate(auroc_diff = AUROC_test - AUROC,
                  aupr_diff = AUPR_test - AUPR) %>%
    dplyr::filter(data_fit == "test",
                  shuffle == "shuffle_false") 

head(all_classification_df, 2)

metric <- "AUROC"

label_logic <- abs(all_classification_df$auroc_diff) > 0.08

ggplot(all_classification_df,
       aes(x = AUROC_test, 
           y = AUROC)) +
    geom_point(alpha = 0.7,
               size = 0.8) +
    xlab("pycytominer") +
    ylab("cytominer") +
    geom_hline(yintercept = 0.5,
               linetype = "dashed",
               color = "red",
               alpha = 0.7) +
    geom_vline(xintercept = 0.5,
               linetype = "dashed",
               color = "red",
               alpha = 0.7) +
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dashed",
                color = "black",
                alpha = 0.7) +
    ggtitle(metric) +
    geom_text_repel(data = subset(all_classification_df, label_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    box.padding = 0.6,
                    point.padding = 0.3,
                    segment.size = 0.2,
                    segment.alpha = 0.6,
                    size = 1.2,
                    fontface = "italic",
                    aes(label = target,
                        x = AUROC_test,
                        y = AUROC)) +
    theme_bw()

outfile <- file.path(
    cytominer_compare_dir,
    paste0("compare_pycytominer_cytominer_", metric, "_", consensus, ".png")
)
ggsave(outfile, height = 3, width = 3, dpi = 300)

metric <- "AUPR"

label_logic <- abs(all_classification_df$aupr_diff) > 0.08

ggplot(all_classification_df,
       aes(x = AUPR_test, 
           y = AUPR)) +
    geom_point(alpha = 0.7,
               size = 0.8) +
    xlab("pycytominer") +
    ylab("cytominer") +
    geom_abline(intercept = 0,
                slope = 1,
                linetype = "dashed",
                color = "black",
                alpha = 0.7) +
    ggtitle(metric) +
    geom_text_repel(data = subset(all_classification_df, label_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    box.padding = 0.6,
                    point.padding = 0.3,
                    segment.size = 0.2,
                    segment.alpha = 0.6,
                    size = 1.2,
                    fontface = "italic",
                    aes(label = target,
                        x = AUPR_test,
                        y = AUPR)) +
    theme_bw()

outfile <- file.path(
    cytominer_compare_dir,
    paste0("compare_pycytominer_cytominer_", metric, "_", consensus, ".png")
)
ggsave(outfile, height = 3, width = 3, dpi = 300)

all_classification_df %>%
    tidyr::drop_na() %>%
    dplyr::mutate(percent_auroc_pycytominer_better = sum(auroc_diff > 0) / n(),
                  percent_aupr_pycytominer_better = sum(aupr_diff > 0) / n()) %>%
    dplyr::distinct(data_fit, shuffle, percent_auroc_pycytominer_better, percent_aupr_pycytominer_better)
