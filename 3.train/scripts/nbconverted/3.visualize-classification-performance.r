suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))

consensus <- "modz"

results_dir <- "results"
figure_dir <- file.path("figures", "classification")
individual_fig_dir <- file.path(
    "figures",
    "individual_target_performance",
    "classification",
    consensus
)

dir.create(results_dir, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_fig_dir, recursive = TRUE, showWarnings = FALSE)

y_file <- file.path(
    results_dir,
    paste0("full_cell_health_y_labels_", consensus, ".tsv.gz")
)

y_df <- readr::read_tsv(y_file, col_types = readr::cols())
print(dim(y_df))
head(y_df)

y_binary_df <- y_df %>%
    dplyr::filter(y_transform == "binarize",
                  shuffle == "shuffle_false",
                  y_type == "y_true") %>%
    dplyr::group_by(target, data_type) %>%
    dplyr::mutate(pos_count = sum(recode_target_value),
                  pos_prop = (pos_count / dplyr::n())) %>%
    dplyr::distinct(target, data_type, pos_count, pos_prop) %>%
    dplyr::select(target, data_type, pos_count, pos_prop) %>%
    dplyr::arrange(pos_count, target)

target_order <- y_binary_df %>%
    dplyr::filter(data_type == "test") %>%
    dplyr::arrange(pos_count) %>%
    dplyr::pull(target)

y_binary_df$target <- factor(y_binary_df$target,
                             levels = target_order)

y_binary_df$data_type <- y_binary_df$data_type %>%
    dplyr::recode("test" = "Test", "train" = "Train")

head(y_binary_df)

min_pos_prop <- 0.02

cutoff_gg <- ggplot(y_binary_df, aes(x = target,
                        y = pos_prop,
                        fill = pos_count)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_continuous(name = "Positive\nClass\nCount") +
    xlab("Target") +
    ylab("Positive Proportion") +
    facet_grid("~data_type") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 90),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +  
    geom_hline(yintercept = min_pos_prop, linetype = "dashed", color = "red")

file <- file.path(figure_dir,
                  paste0("classification_model_cutoff_", consensus, ".png"))
ggsave(file, dpi = 300, width = 6, height = 6)
cutoff_gg

# Define which targets to use downstream based on cutoff
use_targets <- paste(
    y_binary_df %>%
        dplyr::filter(data_type == "Test",
                      pos_prop >= min_pos_prop) %>%
        dplyr::distinct(target) %>%
        dplyr::pull(target)
    )

print(paste("After filtering models by minumum number of samples, we have remaining:", length(use_targets)))

roc_file <- file.path(results_dir,
                      paste0("full_cell_health_roc_results_", consensus, ".tsv.gz"))
full_roc_df <- readr::read_tsv(roc_file, col_types = readr::cols()) %>%
    dplyr::filter(cell_line == "all", target %in% use_targets) 

pr_file <- file.path(results_dir,
                     paste0("full_cell_health_pr_results_", consensus, ".tsv.gz"))
full_pr_df <- readr::read_tsv(pr_file, col_types = readr::cols()) %>%
    dplyr::filter(cell_line == "all", target %in% use_targets)

coef_file <- file.path(results_dir,
                       paste0("full_cell_health_coefficients_", consensus, ".tsv.gz"))
full_coef_df <- readr::read_tsv(coef_file, col_types = readr::cols())

auroc_df <- full_roc_df %>%
    dplyr::distinct(metric, target, auc, data_fit, shuffle, y_transform, min_class_count)

aupr_df <- full_pr_df %>%
    dplyr::distinct(metric, target, auc, data_fit, shuffle, y_transform, min_class_count)

auc_df <- dplyr::bind_rows(auroc_df, aupr_df)

auc_df$metric <- dplyr::recode_factor(
    auc_df$metric,
    "roc" = "AUROC",
    "aupr" = "AUPR"
)

target_order <- auc_df %>%
    dplyr::filter(metric == "AUROC",
                  data_fit == "test",
                  shuffle == "shuffle_false") %>%
    dplyr::arrange(desc(auc))

target_order <- paste(target_order$target)

auc_df$target <- factor(auc_df$target, levels = rev(target_order))

head(auc_df, 3)

summary_gg <- ggplot(auc_df,
       aes(x = target,
           y = auc)) +
    geom_point(aes(size = min_class_count,
                   fill = data_fit,
                   alpha = shuffle),
               pch = 21) +
    geom_hline(data = subset(auc_df,metric == "AUROC"),
               aes(yintercept = 0.5),
               linetype = "dashed",
               color="black",
               alpha = 0.9) +
    ylab("AUC") +
    xlab("") +
    coord_flip() +
    scale_fill_manual(name = "Data Split",
                      values = c("train" = "#0FA3AD",
                                 "test" = "#F76916"),
                      labels = c("train" = "Train",
                                 "test" = "Test")) +
    scale_alpha_manual(name = "Shuffled",
                       values = c("shuffle_false" = 1,
                                  "shuffle_true" = 0.2),
                       labels = c("shuffle_false" = "False",
                                  "shuffle_true" = "True")) +
    scale_size_continuous(name = "Num Positive\nTraining n =",
                          range = c(0.5, 2)) +
    facet_wrap(~metric, nrow = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 90),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 9),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    guides(
        fill = guide_legend(order = 1),
        alpha = guide_legend(order = 2),
        size = guide_legend(order = 3)
    )

file <- file.path(figure_dir,
                  paste0("classification_summary_", consensus, ".png"))
ggsave(file, dpi = 300, width = 6, height = 6)
summary_gg

full_roc_df$target <- factor(full_roc_df$target, levels = target_order)

ggplot(full_roc_df,
       aes(x = fpr, y = tpr)) +
    coord_fixed() +
    facet_wrap(~target) +
    geom_step(aes(color = data_fit,
                  linetype = shuffle),
              alpha = 0.6) +
    geom_abline(intercept = 0,
                slope = 1,
                alpha = 0.5,
                linetype = "dashed") +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    scale_color_manual(name = "",
                       values = c("train" = "#0FA3AD",
                                  "test" = "#F76916"),
                       labels = c("train" = "Train",
                                  "test" = "Test")) +
    scale_linetype_manual(name = "Shuffled Data",
                          values = c("shuffle_true" = "dotted",
                                     "shuffle_false" = "solid"),
                          labels = c("shuffle_true" = "Shuffled",
                                     "shuffle_false" = "Real")) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90),
        strip.text = element_text(size = 3.5),
        strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")
    )

file <- file.path(figure_dir, paste0("roc_curves_", consensus, ".png"))
ggsave(file, dpi = 300, width = 9, height = 9)

full_pr_df$target <- factor(full_pr_df$target, levels = target_order)

ggplot(full_pr_df, aes(x = recall, y = precision)) +
    coord_fixed() +
    facet_wrap(~target) +
    geom_step(aes(color = data_fit,
                  linetype = shuffle),
              alpha = 0.7) +
    xlab("Recall") +
    ylab("Precision") +
    scale_color_manual(name = "",
                       values = c("train" = "#0FA3AD",
                                  "test" = "#F76916"),
                       labels = c("train" = "Train",
                                  "test" = "Test")) +
    scale_linetype_manual(name = "Shuffled Data",
                          values = c("shuffle_true" = "dotted",
                                     "shuffle_false" = "solid"),
                          labels = c("shuffle_true" = "Shuffled",
                                     "shuffle_false" = "Real")) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90),
        strip.text = element_text(size = 3.5),
        strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")
    )

file <- file.path(figure_dir, paste0("pr_curves_", consensus, ".png"))
ggsave(file, dpi = 300, width = 9, height = 9)

y_binary_subset_true_df <- y_df %>%
    dplyr::filter(y_transform == "binarize",
                  y_type == "y_true")

y_binary_subset_pred_df <- y_df %>%
    dplyr::filter(y_transform == "binarize",
                  y_type == "y_pred")

y_plot_df <- y_binary_subset_true_df %>%
    dplyr::inner_join(y_binary_subset_pred_df,
                      by = c("Metadata_profile_id",
                             "target",
                             "data_type",
                             "shuffle",
                             "y_transform"),
                      suffix = c("_true", "_pred"))

y_plot_df$data_type <- dplyr::recode(y_plot_df$data_type,
                                     "train" = "Train",
                                     "test" = "Test")

head(y_plot_df, 3)

label_thresh_value = 0.925

pdf_file <- file.path(
    figure_dir,
    paste0("all_classification_performance_metrics_", consensus, ".pdf")
)
pdf(pdf_file, width = 7, height = 6, onefile = TRUE)

for (target in use_targets) {
    subset_roc_df <- full_roc_df %>%
        dplyr::filter(target == !!target)

    subset_pr_df <- full_pr_df %>%
        dplyr::filter(target == !!target)
    
    # Plot ROC Curves
    roc_gg <-
        ggplot(subset_roc_df,
               aes(x = fpr,
                   y = tpr)) +
        geom_step(aes(color = data_fit,
                      linetype = shuffle),
                  alpha = 0.6) +
        geom_abline(intercept = 0,
                    slope = 1,
                    alpha = 0.5,
                    linetype = "dashed") +
        xlab("False Positive Rate") +
        ylab("True Positive Rate") +
        scale_color_manual(name = "Fit:",
                           values = c("train" = "#0FA3AD",
                                      "test" = "#F76916"),
                           labels = c("train" = "Train",
                                      "test" = "Test"),
                           guide = FALSE) +
        scale_linetype_manual(name = "Data:",
                              values = c("shuffle_true" = "dotted",
                                         "shuffle_false" = "solid"),
                              labels = c("shuffle_true" = "Shuffled",
                                         "shuffle_false" = "Real")) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 9),
              legend.box = "vertical",
              axis.text.x = element_text(angle = 90),
              strip.text = element_text(size = 3.5),
              strip.background = element_rect(colour = "black",
                                              fill = "#fdfff4"))
    # Plot PR Curves
    pr_gg <- ggplot(subset_pr_df,
                    aes(x = recall,
                        y = precision)) +
        geom_step(aes(color = data_fit,
                      linetype = shuffle),
                  alpha = 0.7) +
        xlab("Recall") +
        ylab("Precision") +
        ylim(c(0, 1)) +
        xlim(c(0, 1)) +
        scale_color_manual(name = "Fit:",
                           values = c("train" = "#0FA3AD",
                                      "test" = "#F76916"),
                           labels = c("train" = "Train",
                                      "test" = "Test")) +
        scale_linetype_manual(name = "Data",
                              values = c("shuffle_true" = "dotted",
                                         "shuffle_false" = "solid"),
                              labels = c("shuffle_true" = "Shuffled",
                                         "shuffle_false" = "Real"),
                              guide = FALSE) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 9),
              legend.box = "vertical",
              axis.text.x = element_text(angle = 90),
              strip.text = element_text(size = 3.5),
              strip.background = element_rect(colour = "black",
                                              fill = "#fdfff4"))
    
    # Plot Machine Learning Coefficients
    subset_coef_df <- full_coef_df %>%
        dplyr::filter(target == !!target,
                      shuffle == "shuffle_false",
                      y_transform == "binarize") %>%
        dplyr::mutate(weight_rank = row_number(weight))
    
    # Setup labeling thresholds
    non_zero_coef <- subset_coef_df$abs_weight[subset_coef_df$abs_weight > 0]
    label_thresh <- quantile(non_zero_coef, label_thresh_value)
    label_logic <- subset_coef_df$abs_weight > label_thresh

    coef_gg <-
        ggplot(subset_coef_df,
               aes(x = weight_rank,
                   y = weight)) +
        geom_point(size = 0.2,
                   alpha = 0.6) +
        xlab("Weight Rank") +
        ylab("Weight") +
        geom_text_repel(data = subset(subset_coef_df, label_logic),
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
    
    # Isolate feature
    y_binary_subset_df <- y_plot_df %>%
        dplyr::filter(target == !!target)

    y_binary_subset_df$plot_group <- paste0(y_binary_subset_df$recode_target_value_true,
                                            y_binary_subset_df$shuffle)
    y_binary_subset_df$shuffle <- y_binary_subset_df$shuffle %>%
        dplyr::recode("shuffle_true" = "Shuffled", "shuffle_false" = "Real")
    feature_distrib_gg <-
        ggplot(y_binary_subset_df,
               aes(x = factor(recode_target_value_true),
                   y = recode_target_value_pred,
                   color = shuffle,
                   group = plot_group)) +
        scale_color_manual(name = "Shuffled",
                       labels = c("Shuffled" = "True",
                                  "Real" = "False"),
                       values = c("Shuffled" = "#614650",
                                  "Real" = "#62936D")) +
        facet_grid(shuffle~data_type, scales = "free") +
        geom_boxplot(outlier.alpha = 0) +
        geom_point(alpha = 0.5,
                   size = 0.3,
                   position = position_jitterdodge(jitter.width = 0.2)) +
        theme_bw() +
        theme(strip.text = element_text(size = 8),
              strip.background = element_rect(colour = "black",
                                              fill = "#fdfff4")) +
        xlab("Class") +
        ylab("Prediction")
    
    # Build table for plotting AUC
    auroc_auc_df <- subset_roc_df %>%
        dplyr::distinct(auc, data_fit, shuffle, min_class_count) %>%
        dplyr::rename(auroc = auc) %>%
        dplyr::left_join(
            subset_pr_df %>%
                dplyr::distinct(auc, data_fit, shuffle, min_class_count) %>%
                dplyr::rename(aupr = auc),
            by = c("data_fit", "shuffle", "min_class_count")
        ) %>%
        dplyr::select(auroc, aupr, data_fit, shuffle, min_class_count) %>%
        dplyr::mutate(auroc = round(auroc, 2), aupr = round(aupr, 2))

    auroc_auc_df$shuffle <- dplyr::recode(auroc_auc_df$shuffle,
                                          shuffle_true = "True",
                                          shuffle_false = "False")
    auroc_auc_df$data_fit <- dplyr::recode(auroc_auc_df$data_fit,
                                           train = "Train",
                                           test = "Test")
    auroc_auc_df <- auroc_auc_df %>%
        dplyr::arrange(shuffle) %>%
        dplyr::rename(pos_n = min_class_count,
                      fit = data_fit,
                      AUROC = auroc,
                      AUPR = aupr)
    
    # Plot all performance metrics together with cowplot
    table_theme <- gridExtra::ttheme_default(
        core = list(fg_params=list(cex = 0.4)),
        colhead = list(fg_params=list(cex = 0.5))
    )

    table_gg <- gridExtra::tableGrob(auroc_auc_df,
                                     theme = table_theme,
                                     rows = NULL)
    color_legend_gg <- cowplot::get_legend(
        roc_gg + theme(legend.background = element_rect(fill = "transparent"))
    )
    line_legend_gg <- cowplot::get_legend(
        pr_gg + theme(legend.background = element_rect(fill = "transparent"))
    )
    top_right_panel_gg <- cowplot::plot_grid(
        color_legend_gg,
        line_legend_gg,
        table_gg,
        NULL,
        nrow = 4,
        rel_heights = c(0.4, 0.4, 1, 0.3)
    )
    
    top_row_perf_gg <- cowplot::plot_grid(
        roc_gg + theme(legend.position = "none"),
        pr_gg + theme(legend.position = "none"),
        top_right_panel_gg,
        ncol = 3,
        labels = c("a", "b", ""),
        rel_widths = c(1, 1, 0.8)
    )

    bottom_row_perf_gg <- cowplot::plot_grid(
        coef_gg,
        feature_distrib_gg,
        ncol = 2,
        labels = c("c", "d")
    )

    performance_gg <- cowplot::plot_grid(
        top_row_perf_gg,
        bottom_row_perf_gg,
        nrow = 2,
        ncol = 1,
        labels = c("", "")
    )

    target_title <- ggdraw() + 
      draw_label(
        paste("Performance:", target),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        plot.margin = margin(0, 0, 0, 7)
      )
    
    performance_gg <- cowplot::plot_grid(
        target_title,
        performance_gg,
        ncol = 1,
        rel_heights = c(0.1, 1)
    )
    
    # Save figure
    cowplot_file <- file.path(
        individual_fig_dir,
        paste0(target, "_performance_", consensus, ".png")
    )
    
    cowplot::save_plot(filename = cowplot_file,
                       plot = performance_gg,
                       base_height = 6,
                       base_width = 8)
    
    print(performance_gg)
}

dev.off()
