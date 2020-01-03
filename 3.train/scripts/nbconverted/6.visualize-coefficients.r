suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))

consensus <- "median"

process_coefficients <- function(coef_file) {
    coef_df <- readr::read_tsv(
        coef_file,
        col_types=readr::cols(
            .default = readr::col_double(),
            feature = readr::col_character(),
            target = readr::col_character(),
            y_transform = readr::col_character(),
            shuffle = readr::col_character())
    )
    
    number_models_per_feature <- rowSums(coef_df > 0) %>% tibble::as_tibble()
    number_models_per_feature <- cbind(coef_df$feature, number_models_per_feature)
    colnames(number_models_per_feature) <- c("feature", "num_models_above_zero")
    
    coef_sums_df <- abs(coef_df %>% dplyr::select(-feature, -target, -y_transform, -shuffle)) %>%
        base::rowSums() %>% tibble::enframe() %>%
        dplyr::mutate(feature = coef_df$feature) %>%
        dplyr::select(-name) %>%
        dplyr::rename(value_sum = value) %>%
        dplyr::arrange(desc(value_sum)) %>%
        dplyr::left_join(number_models_per_feature, by = "feature")
    
    split_coef_df <- coef_sums_df %>%
        tidyr::separate(feature,
                        into=c("compartment",
                               "feature_group",
                               "measurement",
                               "channel", 
                               "parameter1", 
                               "parameter2"), sep="_") %>%
        dplyr::mutate(feature_id = coef_sums_df$feature)

    split_coef_df$feature_id <- factor(split_coef_df$feature_id, levels=split_coef_df$feature_id)
    
    return(split_coef_df)
}

coef_file <- file.path("results",
                       paste0("full_cell_health_coefficients_", consensus, ".tsv.gz"))
coef_df <- process_coefficients(coef_file)

head(coef_df, 10)

coef_file <- file.path("results", "all_model_coefficients_shuffled.tsv")
coef_shuffle_df <- process_coefficients(coef_file)

head(coef_shuffle_df, 2)

real_hist_gg <- ggplot(coef_df,
                       aes(x = value_sum)) +
    geom_histogram(bins = 100) +
    theme_bw() +
    xlab("Coefficient Sum (Real)") +
    ylab("Count")

shuffled_hist_gg <- ggplot(coef_shuffle_df, aes(x=value_sum)) +
    geom_histogram(bins = 100) +
    theme_bw() +
    xlab("Coefficient Sum (Shuffled)") +
    ylab("Count")

hist_gg <- cowplot::plot_grid(
    real_hist_gg,
    shuffled_hist_gg,
    ncol = 1,
    labels = c("a", "b")
)

hist_gg

coef_merge_df <- coef_df %>%
    dplyr::inner_join(coef_shuffle_df,
                      by = c("compartment", "feature_group", "measurement",
                             "channel", "parameter1", "parameter2", "feature_id"),
                      suffix = c("", "_shuffle"))

channels <- c("DNA", "AGP", "ER", "Mito", "RNA")
coef_merge_df[!coef_merge_df$channel %in% channels, "channel"] <- "NA"

head(coef_merge_df, 2)

label_logic <- coef_merge_df$value_sum > 1

feature_distrib_gg <- ggplot(coef_merge_df,
                             aes(x = value_sum, y = value_sum_shuffle)) +
    geom_point(aes(color = compartment),
               alpha = 0.6) +
    xlab("Coefficient Sum") +
    ylab("Shuffled Coefficient Sum") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_text_repel(data = subset(coef_merge_df, label_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    box.padding = 0.4,
                    point.padding = 0.1,
                    segment.size = 0.5,
                    segment.alpha = 0.6,
                    size = 1.2,
                    fontface = "italic",
                    aes(label = feature_id,
                        x = value_sum,
                        y = value_sum_shuffle)) +
    scale_color_manual(name = "Compartment",
                         values = c("Cells" = "#7B1336",
                                    "Cytoplasm" = "#313053",
                                    "Nuclei" = "#65B999")) +

    theme_bw()

feature_distrib_gg

coef_channel_mean <- coef_merge_df %>%
    dplyr::group_by(channel) %>%
    dplyr::mutate(channel_mean= mean(abs(value_sum)), 
                  channel_mean_shuffle = mean(abs(value_sum_shuffle))) %>%
    dplyr::distinct(channel, channel_mean, channel_mean_shuffle) %>%
    reshape2::melt(id.vars = 'channel',
                   variable.name = "channel_type",
                   value.name = "channel_mean")

bar_coef_gg <- ggplot(coef_channel_mean, aes(x = channel, y = channel_mean, fill = channel_type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    xlab("Channel") +
    ylab("Abs. Feature Mean") +
    scale_fill_manual(name = "",
                      labels = c("channel_mean" = "Real Data",
                                 "channel_mean_shuffle" = "Shuffled Data"), 
                      values = c("channel_mean" = "#2B76A4",
                                 "channel_mean_shuffle" = "#FCAA2A")) +
    theme_bw()

bar_coef_gg

feature_type_gg <- cowplot::plot_grid(
    bar_coef_gg,
    feature_distrib_gg,
    ncol = 1,
    labels = c("c", "d")
)

all_coef_gg <- cowplot::plot_grid(
    hist_gg,
    feature_type_gg,
    nrow = 1,
    rel_widths = c(0.6, 1)
)

all_coef_gg

cowplot_file <- file.path("figures", "model_coefficient_summary.png")

cowplot::save_plot(filename = cowplot_file,
                           plot = all_coef_gg,
                           base_height = 6,
                           base_width = 8)

coef_all_df <- coef_df %>%
    dplyr::mutate(shuffle = "shuffle_false") %>%
    dplyr::bind_rows(
        coef_shuffle_df %>%
        dplyr::mutate(shuffle = "shuffle_true")
    )

print(dim(coef_all_df))
head(coef_all_df, 3)

all_feature_gg <- ggplot(coef_all_df,
       aes(x = feature_id,
           y = value_sum)) +
    geom_bar(aes(fill = num_models_above_zero), stat="identity") +
    scale_fill_continuous(name = "Number of Models") +
    theme_bw() + 
    facet_wrap(~shuffle, ncol = 1) +
    theme(axis.text.x = element_text(size = 5, angle = 90),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

output_file <- file.path("figures", "coefficient_sum_full.png")
ggsave(output_file, all_feature_gg, height = 10, width = 18)

all_feature_gg

subset_feature_gg <- ggplot(coef_all_df %>%
                            dplyr::top_n(50, value_sum),
       aes(x = feature_id,
           y = value_sum)) +
    geom_bar(aes(fill = num_models_above_zero), stat="identity") +
    scale_fill_continuous(name = "Number of Models") +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 5, angle = 90),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    coord_flip()

output_file <- file.path("figures", "coefficient_sum_subset.png")
ggsave(output_file, subset_feature_gg, height = 8, width = 8)

subset_feature_gg
