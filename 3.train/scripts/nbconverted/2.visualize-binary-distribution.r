suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

consensus <- "median"

# Load Data
y_cols <- readr::cols(
    Metadata_profile_id = readr::col_character(),
    recode_target_value = readr::col_double(),
    target = readr::col_character(),
    data_type = readr::col_character(),
    shuffle = readr::col_character(),
    y_transform = readr::col_character(),
    y_type = readr::col_character()
)

y_file <- file.path("results",
                    paste0("full_cell_health_y_labels_", consensus, ".tsv.gz"))
y_df <- readr::read_tsv(y_file,
                        col_types = y_cols)

y_binary_df <- y_df %>%
    dplyr::filter(shuffle == "shuffle_false",
                  y_transform == "binarize",
                  y_type == "y_true")

y_raw_scores_df <- y_df %>%
    dplyr::filter(shuffle == "shuffle_false",
                  y_transform == "raw",
                  y_type == "y_true")

# Process data for plotting
y_plot_df <- y_raw_scores_df %>%
    dplyr::inner_join(y_binary_df,
                      by = c("Metadata_profile_id",
                             "target",
                             "data_type",
                             "shuffle",
                             "y_type"),
                      suffix = c("_raw", "_binary"))

y_plot_df$data_type <- dplyr::recode(y_plot_df$data_type,
                                     "train" = "Train",
                                     "test" = "Test")

head(y_plot_df, 3)

# Generate and save figures
pdf_file <- file.path("figures",
                      paste0("all_binary_distributions_", consensus, ".pdf"))
pdf(pdf_file, width = 5, height = 3.5, onefile = TRUE)

for (target in unique(y_plot_df$target)) {
    y_plot_subset_df = y_plot_df %>%
        dplyr::filter(target == !!target)

    target_gg <- 
        ggplot(y_plot_subset_df,
               aes(x = recode_target_value_raw,
                   fill = as.factor(recode_target_value_binary))) +
            geom_histogram(bins = 50, alpha = 0.6) +
            facet_grid(~ data_type,
                       scales = "free_y") +
            scale_fill_manual(name = "Binary\nRecoding",
                              labels = c("0" = "0", "1" = "1"),
                              values = c("0" = "#AEA367", "1" = "#403019")) +
            xlab(target) +
            ylab("Count") +
            theme_bw() +
            ggtitle(target) +
            theme(axis.text = element_text(size = 8),
                  axis.title = element_text(size = 9),
                  strip.text = element_text(size = 7),
                  legend.title = element_text(size = 8),
                  title = element_text(size = 12),
                  strip.background = element_rect(colour = "black",
                                                  fill = "#fdfff4"))

    output_file <- file.path("figures",
                             "feature_distribution",
                             paste0(target, "_dist_", consensus, ".png"))
    
    ggsave(filename = output_file,
           plot = target_gg,
           width = 5,
           height = 3.5,
           dpi = 400)
    
    print(target_gg)
}

dev.off()
