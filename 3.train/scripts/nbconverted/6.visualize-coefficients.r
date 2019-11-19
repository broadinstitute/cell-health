suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))

process_coefficients <- function(coef_file) {
    coef_df <- readr::read_tsv(coef_file,
                               col_types=readr::cols(
                                   .default = readr::col_double(),
                                   features = readr::col_character())
                              )
    
    coef_sums_df <- abs(coef_df %>% dplyr::select(-features)) %>%
        base::rowSums() %>% tibble::enframe() %>%
        dplyr::mutate(features = coef_df$features) %>%
        dplyr::select(-name) %>%
        dplyr::rename(value_sum = value) %>%
        dplyr::arrange(desc(value_sum))
    
    split_coef_df <- coef_sums_df %>%
        tidyr::separate(features,
                        into=c("compartment",
                               "feature_group",
                               "measurement",
                               "channel", 
                               "parameter1", 
                               "parameter2"), sep="_") %>%
        dplyr::mutate(feature_id = coef_sums_df$features)

    split_coef_df$feature_id <- factor(split_coef_df$feature_id, levels=split_coef_df$feature_id)
    
    return(split_coef_df)
}

coef_file <- file.path("results", "all_model_coefficients.tsv")
coef_df <- process_coefficients(coef_file)

head(coef_df, 10)

coef_file <- file.path("results", "all_model_coefficients_shuffled.tsv")
coef_shuffle_df <- process_coefficients(coef_file)

head(coef_shuffle_df, 2)

real_hist_gg <- ggplot(coef_df,
                       aes(x=value_sum)) +
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
