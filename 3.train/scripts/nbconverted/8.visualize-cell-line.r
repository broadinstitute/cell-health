suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

set.seed(123)

get_f_stat <- function(df) {
    feature <- c()
    all_f_stat <- c()
    for (health_feature in colnames(df)) {
        if (health_feature != "Metadata_cell_line") {
            result <- eval(
                parse(
                    text = paste0("aov(", health_feature, " ~ Metadata_cell_line, data = df)")
                )
            )
            result <- summary(result)
            f_stat <- result[[1]]$`Pr(>F)`[1]
            all_f_stat <- c(all_f_stat, f_stat)
            feature <- c(feature, health_feature)
        }
    }
    
    f_stat_df <- cbind(feature, all_f_stat) %>% dplyr::as_tibble()
    return(f_stat_df)
}

file <- file.path("results", "all_model_sample_squared_error.tsv")
all_score_error <- readr::read_tsv(file, col_types = readr::cols()) %>%
    dplyr::mutate(data_type = "real")

print(dim(all_score_error))
head(all_score_error, 2)

anova_ready_df <- all_score_error %>%
    dplyr::select(-Metadata_profile_id, -Metadata_gene_name, -Metadata_pert_name, -data_type)

f_stat_df <- get_f_stat(anova_ready_df) %>%
    dplyr::mutate(data_type = "real")
head(f_stat_df, 2)

file <- file.path("results", "all_model_sample_squared_error_shuffled.tsv")
all_shuffle_score_error <- readr::read_tsv(file, col_types = readr::cols()) %>%
    dplyr::mutate(data_type = "shuffled")

print(dim(all_shuffle_score_error))
head(all_shuffle_score_error, 2)

anova_ready_shuffle_df <- all_shuffle_score_error %>%
    dplyr::select(-Metadata_profile_id, -Metadata_gene_name, -Metadata_pert_name, -data_type)

f_stat_shuffle_df <- get_f_stat(anova_ready_shuffle_df) %>%
    dplyr::mutate(data_type = "shuffled")
head(f_stat_shuffle_df, 2)

# Merge Data
full_score_error <- all_score_error %>% 
    dplyr::bind_rows(all_shuffle_score_error)

full_fstat <- f_stat_df %>% 
    dplyr::bind_rows(f_stat_shuffle_df)

mse_mean_cell_line_df <- full_score_error %>%
    dplyr::group_by(Metadata_cell_line, data_type) %>%
    dplyr::summarise_at(vars(starts_with("cc"), starts_with("vb")), mean, na.rm = TRUE) %>%
    reshape2::melt(id.vars = c("Metadata_cell_line", "data_type"),
                   variable.name = "feature",
                   value.name = "mse_mean") %>%
    dplyr::inner_join(full_fstat, by = c("feature", "data_type")) %>%
    dplyr::arrange(desc(mse_mean))

mse_mean_cell_line_df$all_f_stat <- as.numeric(paste(mse_mean_cell_line_df$all_f_stat))
mse_mean_cell_line_df$data_type <- factor(mse_mean_cell_line_df$data_type, levels = c("shuffled", "real"))

head(mse_mean_cell_line_df, 10)

ggplot(mse_mean_cell_line_df,
       aes(x = data_type,
           y = mse_mean,
           color = all_f_stat)) +
    geom_jitter(width = 0.2, size = 0.7, alpha = 0.8, pch = 16) +
    facet_wrap(~Metadata_cell_line, nrow = 3) +
    scale_color_viridis_c(name = "F Stat") +
    xlab("Data Type") +
    ylab("Within Cell Line Mean MSE\nof Cell Health Feature") +
    coord_flip() +
    theme_bw() +
    theme(
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
          strip.text = element_text(size = 8),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

output_file <- file.path("figures", "cell_line_mse_differences.png")
ggsave(output_file, height = 5, width = 5, dpi = 500)
