suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))

consensus <- "modz"

coef_dir <- file.path("figures", "coefficients")
figure_dir <- file.path(coef_dir, consensus)

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
        dplyr::filter(feature_group  == "AreaShape") %>%
        dplyr::group_by(shuffle, compartment, feature_group) %>%
        dplyr::top_n(n = 1, wt = abs_weight)

    compartment_df <- subset_coef_df %>%
        dplyr::filter(feature_group %in% !!compartment_features) %>%
        dplyr::group_by(shuffle, compartment, feature_group, channel) %>%
        dplyr::top_n(n = 1, wt = abs_weight)

    # Process individual feature name info
    total_features <- length(unique(subset_coef_df$feature))
    total_non_zero_features <- nrow(
        subset_coef_df %>%
            dplyr::filter(shuffle == "Real") %>%
            dplyr::filter(abs_weight > 0)
        )
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

    # 1st Plot - Area Features
    area_gg <- ggplot(area_df,
                      aes(x = compartment, y = feature_group)) +
        geom_point(aes(fill = abs_weight), size = 3, pch = 21) +
        facet_wrap(~shuffle) +
        scale_fill_viridis_c(name = "Abs Weight") +
        ylab("") +
        xlab("Compartment") +
        theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
        theme_bw() +
        coef_theme

    # 2nd Plot - Other Compartment Features
    compartment_gg <- ggplot(compartment_df,
                             aes(x = channel, y = feature_group)) +
        geom_point(aes(fill = abs_weight), size = 3, pch = 21) +
        facet_grid(compartment~shuffle) +
        scale_fill_viridis_c(name = "Abs Weight") +
        ylab("Feature Group") +
        xlab("Channel") +
        theme_bw() +
        coef_theme +
        theme(axis.text.x = element_text(angle = 90))

    # 3rd Plot - Individual Feature Names
    feature_name_gg <- ggplot(subset_coef_features_df,
                              aes(x = feature, y = weight)) +
        geom_bar(fill = "#8DB495", color = "black", stat = "identity") +
        ggtitle(feature_title) +
        coord_flip() +
        xlab("") +
        ylab("Model Coefficient") +
        theme_bw() +
        coef_theme

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
    
    # Compile together into single cowplot
    coef_full_gg <- cowplot::plot_grid(
        full_title,
        cowplot::plot_grid(
            feature_name_gg,
            cowplot::plot_grid(
                area_gg,
                compartment_gg,
                nrow = 2,
                rel_heights = c(0.25, 1),
                align = "v",
                axis = "l"
            ),
            ncol = 2
        ),
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
    "Correlation",
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

pdf_file <- file.path(
    coef_dir,
    paste0("all_model_coefficients_", consensus, ".pdf")
)

pdf(pdf_file, width = 7, height = 5, onefile = TRUE)
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
    cowplot::save_plot(output_file, coef_gg, base_height = 5, base_width = 7)
    
    print(coef_gg)
}

dev.off()
