suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))

consensus <- "modz"

coef_dir <- file.path("figures", "coefficients")
figure_dir <- file.path(coef_dir, consensus)
results_dir <- "results"
empty_grob <- cowplot::ggdraw()

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

coef_plot <- function(
    df, target_name, compartment_features, coef_theme, top_plot_num = 15
) {
    # Compile a series of plots that describe model coefficients
    #
    # Arguments:
    # df - coefficient dataframe with feature metadata and weights
    # target_name - a string of a specific cell health model
    # compartment_features - a vector of which features to plot
    # coef_theme - a ggplot theme object to apply to all plots
    # top_plot_num - the number of top individual features to visualize
    
    # Subset the input dataframe to focus on the specific target
    subset_coef_df <- coef_df %>% dplyr::filter(target == !!target_name)

    # Extract and process specific feature sets
    area_df <- subset_coef_df %>%
        dplyr::filter(feature_group %in% c("AreaShape", "Neighbors")) %>%
        dplyr::group_by(shuffle, compartment, feature_group) %>%
        dplyr::top_n(n = 1, wt = abs_weight) %>%
        dplyr::group_by(compartment, feature_group, shuffle) %>%
        dplyr::mutate(abs_max_weight = max(abs_weight)) %>%
        dplyr::select(compartment, feature_group, shuffle, abs_max_weight) %>%
        dplyr::distinct()

    compartment_df <- subset_coef_df %>%
        dplyr::filter(feature_group %in% !!compartment_features) %>%
        dplyr::group_by(shuffle, compartment, feature_group, channel) %>%
        dplyr::top_n(n = 1, wt = abs_weight) %>%
        dplyr::group_by(compartment, feature_group, channel, shuffle) %>%
        dplyr::mutate(abs_max_weight = max(abs_weight)) %>%
        dplyr::select(compartment, feature_group, channel, shuffle, abs_max_weight) %>%
        dplyr::distinct()

    correlation_df <- subset_coef_df %>%
        dplyr::filter(feature_group == "Correlation") %>%
        dplyr::group_by(channel, parameter1, compartment, shuffle) %>%
        dplyr::mutate(abs_max_weight = max(abs_weight)) %>%
        dplyr::select(channel, parameter1, compartment, shuffle, abs_max_weight) %>%
        dplyr::distinct()

    # Process individual feature name info
    total_features <- length(unique(subset_coef_df$feature))
    total_non_zero_features <- nrow(
        subset_coef_df %>%
            dplyr::filter(shuffle == "Real") %>%
            dplyr::filter(abs_weight > 0)
        )
    
    top_plot_num <- ifelse(total_non_zero_features < top_plot_num,
                           total_non_zero_features,
                           top_plot_num)
    
    feature_title <- paste0(
        round((total_non_zero_features / total_features) * 100, 2),
        "% Non-Zero"
    )

    top_n_features <- subset_coef_df %>%
        dplyr::filter(shuffle == "Real") %>%
        dplyr::top_n(n = top_plot_num, wt = abs_weight) %>%
        dplyr::pull(feature)

    feature_order <- subset_coef_df %>%
        dplyr::filter(shuffle == "Real") %>%
        dplyr::arrange(weight) %>%
        dplyr::pull(feature)

    subset_coef_features_df <- subset_coef_df %>%
        dplyr::filter(shuffle == "Real") %>%
        dplyr::filter(feature %in% !!top_n_features)

    subset_coef_features_df$feature <- factor(
        subset_coef_features_df$feature, levels = feature_order
    )

    min_gradient <- min(subset_coef_df$abs_weight)
    max_gradient <- max(subset_coef_df$abs_weight)

    # Draw an empty plot for later
    empty_grob <- ggdraw()

    # 1st Plot - Area Features
    # First, create an area background to fill missing elements
    area_comparments <- unique(area_df$compartment)
    area_feature_groups <- unique(area_df$feature_group)
    area_shuffle <- unique(area_df$shuffle)
    area_background <- tidyr::crossing(area_comparments, area_feature_groups, area_shuffle)
    colnames(area_background) <- c("compartment", "feature_group", "shuffle")
    
    # Now plot
    area_gg <- ggplot(area_df, aes(x = compartment, y = feature_group)) +
        geom_point(data = area_background, fill="grey", size = 3, pch = 22) +
        geom_point(aes(fill = abs_max_weight), size = 4, pch = 22) +
        facet_wrap(~shuffle) +
        scale_fill_gradient2(
            name = "Max\nAbs. Weight",
            low = "white",
            high = "blue",
            limits = c(min_gradient, max_gradient)
        ) +
        ylab("Feature Group") +
        xlab("Compartment") +
        theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
        theme_bw() +
        coord_fixed() +
        geom_text(aes(label = round(abs_max_weight, 2)), size = 1.5) +
        coef_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # 2nd Plot - Other Compartment Features
    # First, create a background to fill missing elements
    comp_comparments <- unique(compartment_df$compartment)
    comp_channels <- unique(compartment_df$channel)
    comp_feature_groups <- unique(compartment_df$feature_group)
    comp_shuffle <- unique(compartment_df$shuffle)
    comp_background <- tidyr::crossing(comp_comparments, comp_channels, comp_feature_groups, comp_shuffle)
    colnames(comp_background) <- c("compartment", "channel", "feature_group", "shuffle")
    
    na_legend <- cowplot::get_legend(
        ggplot(comp_background %>% dplyr::filter(shuffle == "Real"),
               aes(x = channel, y = compartment)) +
        geom_point(aes(fill = shuffle), size = 4, pch = 22) +
        scale_fill_manual(name = "", values = "grey", labels = "N/A") +
        coef_theme
        ) 
    
    # Now, plot
    compartment_gg <- ggplot(compartment_df,
                             aes(x = channel, y = feature_group)) +
        geom_point(data = comp_background, fill="grey", size = 3, pch = 22) +
        geom_point(aes(fill = abs_max_weight), size = 4, pch = 22) +
        facet_grid(compartment~shuffle) +
        scale_fill_gradient2(
            name = "Max\nAbs. Weight",
            low = "white",
            high = "blue",
            limits = c(min_gradient, max_gradient)
        ) +
        ylab("Feature Group") +
        xlab("Channel") +
        coord_fixed() +
        theme_bw() +
        geom_text(aes(label = round(abs_max_weight, 2)), size = 1.5) +
        coef_theme +
        theme(axis.text.x = element_text(angle = 90))

    # 3rd Plot - Correlation Features
    # First, create a background to fill missing elements
    corr_comparments <- unique(correlation_df$compartment)
    corr_channels <- unique(correlation_df$channel)
    corr_feature_groups <- unique(correlation_df$parameter1)
    corr_shuffle <- unique(correlation_df$shuffle)
    corr_background <- tidyr::crossing(corr_comparments, corr_channels, corr_feature_groups, corr_shuffle)
    colnames(corr_background) <- c("compartment", "channel", "parameter1", "shuffle")
    
    # Now, plot
    correlation_gg <- ggplot(correlation_df,
                             aes(x = channel, y = parameter1)) +
        geom_point(data = corr_background, fill="grey", size = 3, pch = 22) +
        geom_point(aes(fill = abs_max_weight), size = 4, pch = 22) +
        facet_wrap(~shuffle) +
        scale_fill_gradient2(
            name = "Max\nAbs. Weight",
            low = "white",
            high = "blue",
            limits = c(min_gradient, max_gradient)
        ) +
        ylab("Channel Colocalization") +
        xlab("Channel Colocalization") +
        facet_grid(~compartment) + 
        coord_fixed() +
        geom_text(aes(label = round(abs_max_weight, 2)), size = 1.5) +
        theme_bw() +
        coef_theme +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
        facet_grid(compartment~shuffle) 

    # 4th Plot - Individual Feature Names
    feature_name_gg <- ggplot(subset_coef_features_df,
                              aes(x = feature, y = weight, fill = abs(weight))) +
        geom_bar(color = "black", stat = "identity") +
        scale_fill_gradient2(
                name = "Abs. Weight",
                low = "white",
                high = "blue",
                limits = c(min_gradient, max_gradient)
            ) +
        ggtitle(feature_title) +
        coord_flip() +
        xlab("") +
        ylab("Model Coefficient") +
        theme_bw() +
        coef_theme +
        theme(legend.position = "none")

    # Get cowplot title
    use_title <- label_df %>%
        dplyr::filter(id == !!target) %>%
        dplyr::pull(readable_name)

    full_title <- ggdraw() + 
      draw_label(
          use_title,
          fontface = 'bold',
          x = 0,
          hjust = -0.1
      )
    coef_legend <- cowplot::get_legend(correlation_gg)
    full_legend <- cowplot::plot_grid(
        empty_grob,
        coef_legend,
        na_legend,
        empty_grob,
        nrow = 4,
        rel_heights = c(1, 0.4, 0.4, 1)
    )
    
    right_panel <- cowplot::plot_grid(
        area_gg + theme(legend.position = "none",
                        plot.margin = margin(1, 1, 1, 1)),
        cowplot::plot_grid(
            correlation_gg + theme(legend.position = "none",
                                   plot.margin = margin(1, 1, 1, 1)),
            compartment_gg + theme(legend.position = "none",
                                   plot.margin = margin(1.3, 1, 1, 1)),
            nrow = 2,
            align = "vh",
            axis = "1"
        ),
        rel_heights = c(0.2, 1),
        nrow = 2,
        align = "none"
    )

    full_panel <- cowplot::plot_grid(
        feature_name_gg + theme(plot.margin = margin(1, 1, 1, 1)),
        right_panel,
        align = "none"
    )

    coef_full_gg <- cowplot::plot_grid(
        full_panel,
        full_legend,
        rel_widths = c(1, 0.1)
    )
    
    coef_full_gg <- cowplot::plot_grid(
        full_title,
        coef_full_gg,
        nrow = 2,
        rel_heights = c(0.1, 1)
        )
    
    return(coef_full_gg)
}

# Annotated Cell Health Features
feat_file <- file.path(
    "..",
    "1.generate-profiles",
    "data",
    "labels",
    "feature_mapping_annotated.csv"
)

label_df <- readr::read_csv(feat_file, col_types = readr::cols())

# Load Regression Results
regression_file <- file.path(
    results_dir, 
    paste0("full_cell_health_regression_", consensus, ".tsv.gz")
)

regression_metrics_df <- readr::read_tsv(regression_file, col_types = readr::cols()) %>%
    dplyr::filter(cell_line == "all", metric == "r_two", data_fit == "test")

regression_metrics_df$shuffle <- dplyr::recode(
    regression_metrics_df$shuffle,
    "shuffle_true" = "Permuted",
    "shuffle_false" = "Real"
)
regression_metrics_df$shuffle <- factor(
    regression_metrics_df$shuffle,
    levels = c("Real", "Permuted")
)

print(dim(regression_metrics_df))
head(regression_metrics_df)

coef_file <- file.path(
    "results", paste0("full_cell_health_coefficients_", consensus, ".tsv.gz")
)

coef_df <- readr::read_tsv(
    coef_file,
    col_types = readr::cols(
        feature = readr::col_character(),
        weight = readr::col_double(),
        abs_weight = readr::col_double(),
        target = readr::col_character(),
        y_transform = readr::col_character(),
        shuffle = readr::col_character()
    )
) %>%
    dplyr::filter(y_transform == "raw") %>%
    tidyr::separate(
        feature,
        into = c(
            "compartment",
            "feature_group",
            "measurement",
            "channel", 
            "parameter1", 
            "parameter2"
        ),
        sep = "_",
        remove = FALSE
    )

coef_df$shuffle <- dplyr::recode(
    coef_df$shuffle,
    "shuffle_true" = "Permuted",
    "shuffle_false" = "Real"
)
coef_df$shuffle <- factor(
    coef_df$shuffle,
    levels = c("Real", "Permuted")
)

print(dim(coef_df))
head(coef_df, 5)

compartment_features <- c(
    "Texture",
    "Intensity",
    "RadialDistribution",
    "Granularity"
)

top_plot_num <- 15

coef_theme <- theme(
    strip.text = element_text(size = 6,
                              color = "black",
                              margin = margin(1, 1, 1, 1)),
    strip.background = element_rect(colour = "black",
                                    fill = "#fdfff4"),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 6),
    plot.title = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm"),
    legend.key.size = unit(0.3, "cm")
)

point_size <- 5.5
text_label_size <- 1.5

pdf_file <- file.path(
    coef_dir,
    paste0("all_model_coefficients_", consensus, ".pdf")
)

pdf(pdf_file, width = 7, height = 8, onefile = TRUE)
for (target in unique(coef_df$target)) {
    coef_gg <- coef_plot(
        df = coef_df,
        target_name = target,
        compartment_features = compartment_features,
        coef_theme = coef_theme,
        top_plot_num = top_plot_num
    )
    output_file <- file.path(
        figure_dir,
        paste0("model_", consensus, "_", target, ".png")
    )
    cowplot::save_plot(output_file, coef_gg, base_height = 7.5, base_width = 6.5)
    
    print(coef_gg)
}

dev.off()

coef_summary_df <- coef_df %>%
    dplyr::left_join(regression_metrics_df, by = c("target", "y_transform", "shuffle")) %>%
    dplyr::mutate(weighted_coef = weight * value, abs_weighted_coef = abs(weight * value)) %>%
    dplyr::group_by(feature, shuffle) %>%
    dplyr::mutate(
        abs_max_weight = max(abs_weight),
        abs_mean_weight = mean(abs_weight),
        abs_95percentile = quantile(abs_weight, 0.95)
    ) %>%
    dplyr::select(
        feature,
        abs_max_weight,
        abs_mean_weight,
        abs_95percentile,
        compartment,
        feature_group,
        measurement,
        channel,
        parameter1,
        shuffle
    ) %>%
    dplyr::distinct()

head(coef_summary_df)

min_gradient <- 0
point_size <- 8
summary_metrics <- c("abs_mean", "abs_max", "abs_95th")

for (shuffle_option in c("Permuted", "Real")) {
    summary_subset_df <- coef_summary_df %>% dplyr::filter(shuffle == !!shuffle_option)
    
    for (summary_metric in summary_metrics) {

        if (summary_metric == "abs_mean") {
            column_choice <- "abs_mean_weight"
            stat_string <- "aggregated_mean"
            legend_name <- "Abs. Mean\nWeighted Coef\n (x100)"
            max_gradient <- max(summary_subset_df[, column_choice]) * 75
        } else if (summary_metric == "abs_max") {
            column_choice <- "abs_max_weight"
            stat_string <- "aggregated_max"
            legend_name <- "Abs. Max\nWeighted Coef"
            max_gradient <- max(summary_subset_df[, column_choice])
        } else if (summary_metric == "abs_95th") {
            column_choice <- "abs_95percentile"
            stat_string <- "aggregated_95th"
            legend_name <- "Abs. 95th Percent\nWeighted Coef"
            max_gradient <- max(summary_subset_df[, column_choice])
        }

        # Process correlation different from other features
        correlation_df <- summary_subset_df %>%
            dplyr::filter(feature_group == "Correlation") %>%
            dplyr::group_by(channel, parameter1, compartment, shuffle) %>%
            dplyr::mutate(
                aggregated_max = max(abs_max_weight),
                aggregated_mean = mean(abs_mean_weight) * 100,
                aggregated_95th = quantile(abs_95percentile, 0.95)
            ) %>%
            dplyr::select(
                channel, parameter1, compartment, shuffle, aggregated_max, aggregated_mean, aggregated_95th
            ) %>%
            dplyr::distinct() %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
                aggregated_max_round = round(aggregated_max, 2),
                aggregated_mean_round = round(aggregated_mean, 2),
                aggregated_95th_round = round(aggregated_95th, 2)
            )

        corr_comparments <- unique(correlation_df$compartment)
        corr_channels <- unique(correlation_df$channel)
        corr_feature_groups <- unique(correlation_df$parameter1)
        corr_shuffle <- unique(correlation_df$shuffle)
        corr_background <- tidyr::crossing(corr_comparments, corr_channels, corr_feature_groups, corr_shuffle)
        colnames(corr_background) <- c("compartment", "channel", "parameter1", "shuffle")

        correlation_gg <- ggplot(correlation_df, aes(x = channel, y = parameter1)) +
            geom_point(data = corr_background, fill="grey", size = 4, pch = 22) +
            geom_point(aes_string(fill = stat_string), size = point_size, pch = 22) +
            facet_wrap(~shuffle) +
            scale_fill_gradient2(
                name = legend_name,
                low = "white",
                high = "blue",
                limits = c(min_gradient, max_gradient)
            ) +
            ylab("Channel Colocalization") +
            xlab("Channel Colocalization") +
            facet_grid(~compartment) + 
            coord_fixed() +
            geom_text(aes_string(label = paste0(stat_string, "_round")), size = text_label_size) +
            theme_bw() +
            coef_theme +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5))

        # Process area different from other features
        area_df <- summary_subset_df %>%
            dplyr::filter(feature_group %in% c("AreaShape", "Neighbors")) %>%
            dplyr::group_by(compartment, feature_group, shuffle) %>%
            dplyr::mutate(
                aggregated_max = max(abs_max_weight),
                aggregated_mean = mean(abs_mean_weight) * 100,
                aggregated_95th = quantile(abs_95percentile, 0.95)
            ) %>%
            dplyr::select(
                compartment, feature_group, shuffle, aggregated_max, aggregated_mean, aggregated_95th
            ) %>%
            dplyr::distinct() %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
                aggregated_max_round = round(aggregated_max, 2),
                aggregated_mean_round = round(aggregated_mean, 2),
                aggregated_95th_round = round(aggregated_95th, 2)
            )

        area_df$stat_output <- round(area_df[, stat_string], 2)

        area_comparments <- unique(area_df$compartment)
        area_feature_groups <- unique(area_df$feature_group)
        area_shuffle <- unique(area_df$shuffle)
        area_background <- tidyr::crossing(area_comparments, area_feature_groups, area_shuffle)
        colnames(area_background) <- c("compartment", "feature_group", "shuffle")

        area_gg <- ggplot(area_df, aes(x = compartment, y = feature_group)) +
            geom_point(data = area_background, fill = "grey", size = 4, pch = 22) +
            geom_point(aes_string(fill = stat_string), size = point_size, pch = 22) +
            scale_fill_gradient2(
                name = legend_name,
                low = "white",
                high = "blue",
                limits = c(min_gradient, max_gradient)
            ) +
            ylab("Feature Group") +
            xlab("Compartment") +
            geom_text(aes_string(label = paste0(stat_string, "_round")), size = text_label_size) +
            theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
            theme_bw() +
            coord_fixed() +
            coef_theme

        # Process compartment features different from other features
        compartment_df <- summary_subset_df %>%
            dplyr::filter(feature_group %in% !!compartment_features) %>%
            dplyr::group_by(compartment, channel, feature_group, shuffle) %>%
            dplyr::mutate(
                aggregated_max = max(abs_max_weight),
                aggregated_mean = mean(abs_mean_weight) * 100,
                aggregated_95th = quantile(abs_95percentile, 0.95)
            ) %>%
            dplyr::select(
                compartment, channel, feature_group, shuffle, aggregated_max, aggregated_mean, aggregated_95th
            ) %>%
            dplyr::distinct() %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
                aggregated_max_round = round(aggregated_max, 2),
                aggregated_mean_round = round(aggregated_mean, 2),
                aggregated_95th_round = round(aggregated_95th, 2)
            )

        compartment_df$stat_output <- round(compartment_df[, stat_string], 2)

        comp_comparments <- unique(compartment_df$compartment)
        comp_channels <- unique(compartment_df$channel)
        comp_feature_groups <- unique(compartment_df$feature_group)
        comp_shuffle <- unique(compartment_df$shuffle)
        comp_background <- tidyr::crossing(comp_comparments, comp_channels, comp_feature_groups, comp_shuffle)
        colnames(comp_background) <- c("compartment", "channel", "feature_group", "shuffle")

        na_legend <- cowplot::get_legend(
            ggplot(comp_background %>% dplyr::filter(shuffle == "Real"),
                   aes(x = channel, y = compartment)) +
            geom_point(aes(fill = shuffle), size = 4, pch = 22) +
            scale_fill_manual(name = "", values = "grey", labels = "N/A") +
            coef_theme
            ) 

        compartment_gg <- ggplot(compartment_df,
                                 aes(x = channel, y = feature_group)) +
            geom_point(data = comp_background, fill = "grey", size = 4, pch = 22) +
            geom_point(aes_string(fill = stat_string), size = point_size, pch = 22) +
            geom_text(aes_string(label = paste0(stat_string, "_round")), size = text_label_size) +
            facet_grid(~compartment) +
            scale_fill_gradient2(
                name = legend_name,
                low = "white",
                high = "blue",
                limits = c(min_gradient, max_gradient)
            ) +
            ylab("Feature Group") +
            xlab("Channel") +
            coord_fixed() +
            theme_bw() +
            coef_theme +
            theme(axis.text.x = element_text(angle = 90))

        # Compile full panel plot
        bottom_panel <- cowplot::plot_grid(
            correlation_gg + theme(legend.position = "none",
                                   plot.margin = margin(1.3, 1.3, 1.3, 1.3)),
            compartment_gg + theme(legend.position = "none",
                                   plot.margin = margin(1.3, 1.3, 1.3, 1.3)),
            nrow = 2,
            align = "hv"
            )

        main_legend_gg <- cowplot::get_legend(correlation_gg)
        full_legend <- cowplot::plot_grid(
            empty_grob,
            main_legend_gg,
            na_legend,
            empty_grob,
            nrow = 4,
            rel_heights = c(1, 0.4, 0.4, 1)
        )

        big_fig <- cowplot::plot_grid(
            cowplot::plot_grid(
                area_gg + theme(legend.position = "none"),
                bottom_panel,
                nrow = 2,
                rel_heights = c(0.3, 1)
            ),
            full_legend,
            ncol = 2,
            rel_widths = c(1, 0.15)
        )

        print(big_fig)

        output_file <- file.path(
            coef_dir,
            paste0(
                "coefficient_summary_", shuffle_option, "_", consensus, "_", summary_metric, ".png"
            )
        )
        cowplot::save_plot(output_file, big_fig, base_height = 5, base_width = 6, dpi = 500)
    }
}
