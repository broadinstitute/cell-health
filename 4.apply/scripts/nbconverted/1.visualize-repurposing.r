suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source(file.path("scripts", "visualize_utils.R"))

consensus <- "modz"

output_dir <- file.path("figures", "umap", consensus)
dir.create(output_dir)

# Load Data
data_dir <- file.path("repurposing_cellhealth_shiny", "data")
real_file <- file.path(
    data_dir,
    paste0("moa_cell_health_", consensus, ".tsv.gz")
)

cp_embedding_df <- readr::read_tsv(real_file, col_types = readr::cols())

cp_embedding_df <- cp_embedding_df %>%
    dplyr::mutate(Metadata_Treatment = "Compound")

cp_embedding_df$Metadata_Treatment[cp_embedding_df$Metadata_broad_core_id == "DMSO"] = "DMSO"

print(dim(cp_embedding_df))
head(cp_embedding_df, 3)

ggplot(cp_embedding_df,
       aes(x = umap_x, y = umap_y)) +
    geom_point(aes(color = Metadata_dose_recode,
                   size = paste(Metadata_Treatment)),
               pch = 16,
               alpha = 0.6) +
    theme_bw() +
    scale_color_viridis_c(name = "Dose\nRecoded") +
    scale_size_discrete("Treatment",
                        range = c(0.5, 2)) +
    xlab("UMAP 1") +
    ylab("UMAP 2")

output_file <- file.path(
    output_dir,
    paste0("umap_repurposing_cell_painting_dose_consensus_", consensus, ".png")
)
ggsave(output_file, height = 5, width = 6, dpi = 500)

ggplot(cp_embedding_df,
       aes(x = umap_x, y = umap_y)) +
    geom_point(aes(color = Metadata_dose_recode,
                   size = paste(Metadata_Treatment)),
               pch = 16,
               alpha = 0.6) +
    theme_bw() +
    scale_color_viridis_c(name = "Dose\nRecoded") +
    scale_size_discrete("Treatment",
                        range = c(0.5, 2)) +
    xlab("UMAP 1") +
    ylab("UMAP 2")

output_file <- file.path(
    output_dir,
    paste0("umap_repurposing_cell_painting_dose_consensus_", consensus, ".png")
)
ggsave(output_file, height = 5, width = 6, dpi = 500)

ggplot(cp_embedding_df %>% dplyr::filter(Metadata_Treatment == "DMSO"),
       aes(x = umap_x, y = umap_y)) +
    geom_point(aes(color = Metadata_pert_well),
               pch = 16,
               size = 2,
               alpha = 0.6) +
    geom_point(data = cp_embedding_df,
               color = "grey",
               alpha = 0.1,
               size = 0.5) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab("UMAP 1") +
    ylab("UMAP 2")

output_file <- file.path(
    output_dir,
    paste0("umap_repurposing_cell_painting_dose_consensus_dmso_", consensus, ".png")
)

ggsave(output_file, height = 5, width = 6, dpi = 500)

# Load feature mapping
mapping_dir <- file.path("..", "1.generate-profiles", "data", "labels")
mapping_file <- file.path(mapping_dir, "feature_mapping_annotated.csv")
map_df <- readr::read_csv(
    mapping_file,
    col_types = readr::cols(.default = readr::col_character())
)

print(dim(map_df))
tail(map_df, 3)

visualize_umap(
    cp_embedding_df,
    target_variable = "cell_health_modz_target_vb_num_live_cells",
    legend_title = "Num Live Cells",
    output_dir = "none",
    save = FALSE
)

visualize_umap(
    cp_embedding_df,
    target_variable = "cell_health_modz_target_vb_live_cell_width_length",
    legend_title = "Live Cell\n(Width:Length)",
    output_dir = "none",
    save = FALSE
)

visualize_umap(
    cp_embedding_df,
    target_variable = "cell_health_modz_target_vb_live_cell_roundness",
    legend_title = "Live Cell Roundness",
    output_dir = "none",
    save = FALSE
)

visualize_umap(
    cp_embedding_df,
    target_variable = "cell_health_modz_target_cc_all_n_objects",
    legend_title = "Number of Objects",
    output_dir = "none",
    save = FALSE
)

visualize_umap(
    cp_embedding_df,
    target_variable = "cell_health_modz_target_vb_live_cell_area",
    legend_title = "Live Cell Area",
    output_dir = "none",
    save = FALSE
)

visualize_umap(
    cp_embedding_df,
    target_variable = "cell_health_modz_target_cc_cc_n_objects",
    legend_title = "Num Cell\nCycle Objects",
    output_dir = "none",
    save = FALSE
)

visualize_umap(
    cp_embedding_df,
    target_variable = "cell_health_modz_target_cc_g1_n_objects",
    legend_title = "G1 Objects",
    output_dir = "none",
    save = FALSE
)

visualize_umap(
    cp_embedding_df,
    target_variable = "cell_health_modz_target_cc_s_intensity_nucleus_area_mean",
    legend_title = "Sum S phase",
    output_dir = "none",
    save = FALSE
)

visualize_umap(
    cp_embedding_df,
    target_variable = "cell_health_modz_target_vb_ros_back_mean",
    legend_title = "ROS Background",
    output_dir = "none",
    save = FALSE
)

cell_health_variables <- colnames(
    cp_embedding_df %>%
        dplyr::select(starts_with("cell_health_modz_target_"))
    )

length(cell_health_variables)

pdf_file <- file.path(
    output_dir,
    paste0("repurposing_hub_umaps_consensus_", consensus, ".pdf")
)
pdf(pdf_file, width = 5, height = 5, onefile = TRUE)

for (cell_health_variable in cell_health_variables) {
    umap_gg <- visualize_umap(
        df = cp_embedding_df,
        target_variable = cell_health_variable,
        legend_title = "Prediction:",
        title = cell_health_variable,
        dpi = 200,
        save_figure = FALSE
    )
}

dev.off()

assay_theme_file <- file.path("..", "3.train", "scripts", "assay_themes.R")
source(assay_theme_file)

col_types <- readr::cols(
    .default = readr::col_character(),
    shuffle_false = readr::col_double(),
    shuffle_true = readr::col_double()
)

rank_file <- file.path(
    "repurposing_cellhealth_shiny",
    "data",
    paste0("A549_ranked_models_regression_", consensus, ".tsv")
)
model_rank_df <- readr::read_tsv(rank_file, col_types = col_types)

# Recode the target variable
model_rank_df$target <- paste0("cell_health_", consensus, "_target_", model_rank_df$target)

head(model_rank_df, 3)

dmso_embeddings_df <- cp_embedding_df %>%
    dplyr::filter(Metadata_Treatment == "DMSO")

non_dmso_embeddings_df <- cp_embedding_df %>%
    dplyr::filter(Metadata_Treatment != "DMSO")

std_dev_dmso_features <- apply(
    dmso_embeddings_df %>% 
        dplyr::select(matches("cell_health_modz_target")),
    2,
    sd
)
std_dev_compound_features <- apply(
    non_dmso_embeddings_df %>%
        dplyr::select(matches("cell_health_modz_target")),
    2, 
    sd
)

std_dev_all_df <- dplyr::bind_cols(
    as.data.frame(std_dev_dmso_features),
    as.data.frame(std_dev_compound_features)
) %>%
    dplyr::mutate(
        features = colnames(dmso_embeddings_df %>%
                                dplyr::select(matches("cell_health_modz_target")))
    ) %>%
    dplyr::left_join(model_rank_df, by = c("features" = "target")) 

good_performing <- std_dev_all_df %>%
    dplyr::filter(shuffle_false > 0)

bad_performing <- std_dev_all_df %>%
    dplyr::filter(shuffle_false <= 0)

std_dev_good_df <- good_performing %>%
    dplyr::mutate(performance_scaled = (
        good_performing$shuffle_false - min(good_performing$shuffle_false)
    ) / (
        max(good_performing$shuffle_false) - min(good_performing$shuffle_false)
    )
                  )

print(dim(std_dev_good_df))
head(std_dev_good_df, 2)

ggplot(std_dev_good_df,
       aes(x = std_dev_dmso_features, y = std_dev_compound_features)) +
    geom_point(aes(color = assay, size = performance_scaled),
               alpha = 0.8) +
    geom_point(data = bad_performing,
               aes(color = assay),
               size = 1,
               alpha = 0.5) +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    scale_color_manual(name = "Assay",
                       values = dye_colors,
                       labels = dye_labels) +
    scale_size_continuous(name = "Performance Scaled") + 
    xlab("Standard Deviation across DMSO Well Consensus") +
    ylab("Standard Deviation across all Consensus Compound Treatments") +
    theme(axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
         legend.key.size = unit(0.4, "cm"))

output_file <- file.path(output_dir, "dmso_vs_compound_standard_deviation.png")
ggsave(output_file, height = 5, width = 6, dpi = 500)
