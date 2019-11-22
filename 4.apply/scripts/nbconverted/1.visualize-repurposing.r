suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# Load Data
data_dir <- file.path("data", "repurposing_transformed")
real_file <- file.path(data_dir, "repurposing_umap_transformed_cell_painting.tsv.gz")

cp_embedding_df <- readr::read_tsv(real_file, col_types = readr::cols())

cell_health_file <- file.path(data_dir, "repurposing_transformed_real_models.tsv.gz")
cell_health_df <- readr::read_tsv(cell_health_file, col_types = readr::cols()) 

# Merge data together
cp_embedding_df <- cp_embedding_df %>%
    dplyr::left_join(cell_health_df,
                     by = c("Metadata_broad_sample", "Metadata_dose_recode", "Image_Metadata_Well")) %>%
    dplyr::mutate(Metadata_Treatment = cp_embedding_df$Image_Metadata_Well)

cp_embedding_df$Metadata_Treatment[cp_embedding_df$Image_Metadata_Well == "collapsed"] = "Compound"
cp_embedding_df$Metadata_Treatment[cp_embedding_df$Image_Metadata_Well != "collapsed"] = "DMSO"

print(dim(cp_embedding_df))
head(cp_embedding_df, 3)

table(cp_embedding_df$Metadata_dose_recode, cp_embedding_df$Metadata_Treatment)

ggplot(cp_embedding_df,
       aes(x = umap_x, y = umap_y)) +
    geom_point(aes(color = Metadata_dose_recode,
                   size = paste(Metadata_Treatment)),
               pch = 16,
               alpha = 0.6) +
    theme_bw() +
    scale_color_viridis_c(name = "log(dose)") +
    xlim(c(-7.75, 6)) +
    ylim(c(-6.75, 5)) +
    scale_size_discrete("Treatment",
                        range = c(0.5, 3)) +
    xlab("UMAP 1") +
    ylab("UMAP 2")

output_file <- file.path("figures", "umap_repurposing_cell_painting_dose_consensus.png")
ggsave(output_file, height = 5, width = 6, dpi = 500)

ggplot(cp_embedding_df %>% dplyr::filter(Metadata_Treatment == "DMSO"),
       aes(x = umap_x, y = umap_y)) +
    geom_point(aes(color = Image_Metadata_Well),
               pch = 16,
               size = 3,
               alpha = 0.6) +
    geom_point(data = cp_embedding_df, color = "grey", alpha = 0.1, size = 0.5) +
    xlim(c(-7.75, 6)) +
    ylim(c(-6.75, 5)) +
    theme_bw() +
    xlab("UMAP 1") +
    ylab("UMAP 2")

output_file <- file.path("figures", "umap_repurposing_cell_painting_dose_consensus_dmso.png")
ggsave(output_file, height = 5, width = 6, dpi = 500)

visualize_model <- function(target_variable, legend_title, title = "none", dpi = 500, save_figure = TRUE) {
    plot_gg <- ggplot(cp_embedding_df, aes(x = umap_x, y = umap_y)) +
        geom_point(aes_string(color = target_variable),
                   size = 0.5,
                   pch = 16,
                   alpha = 0.6) +
        theme_bw() +
        scale_color_viridis_c(name = legend_title) +
        xlab("UMAP 1") +
        ylab("UMAP 2")
    
    if (title != "none") {
        plot_gg <- plot_gg + ggtitle(title)
    }
    if (save_figure) {
        output_file <- file.path("figures",
                                 paste0("umap_repurposing_cell_painting_",
                                        target_variable,
                                        "_consensus.png"))
        ggsave(output_file, height = 5, width = 6, dpi = dpi)
    }
    
    print(plot_gg)
}

# Load feature mapping
mapping_dir <- file.path("..", "1.generate-profiles", "data", "labels")
mapping_file <- file.path(mapping_dir, "feature_mapping_annotated.csv")
map_df <- readr::read_csv(mapping_file,
                          col_types = readr::cols(.default = readr::col_character()))

print(dim(map_df))
head(map_df, 3)

visualize_model(target_variable = "Metadata_dose_recode",
                legend_title = "Compound\nDose (log)") 

map_df %>% filter(original_name == "# Live Cells")

visualize_model(target_variable = "vb_num_live_cells",
                legend_title = "Num Live Cells\n(DRAQ7)")

map_df %>% filter(original_name == "Live Width:Length")

visualize_model(target_variable = "vb_live_cell_width_length",
                legend_title = "Live Cell\n(Width:Length)\n(DRAQ7)")

map_df %>% filter(original_name == "Live Cell Roundness")

visualize_model(target_variable = "vb_live_cell_roundness",
                legend_title = "Live Cell Roundness\n(DRAQ7)")

map_df %>% filter(original_name == "ALL - Number of Objects")

visualize_model(target_variable = "cc_all_n_objects",
                legend_title = "Number of Objects\n(Hoechst)")

map_df %>% filter(original_name == "Live Cell Area")

visualize_model(target_variable = "vb_live_cell_area",
                legend_title = "Live Cell Area\n(DRAQ7)")

map_df %>% filter(original_name == "CC - Number of Objects")

visualize_model(target_variable = "cc_cc_n_objects",
                legend_title = "Num Cell\nCycle Objects\n(Hoechst)")

map_df %>% filter(original_name == "G1 - Number of Objects")

visualize_model(target_variable = "cc_g1_n_objects",
                legend_title = "G1 Objects\n(Many Dyes)")

map_df %>% filter(original_name == "edu positive - Intensity Nucleus Alexa 647 Sum - Sum per Well")

visualize_model(target_variable = "cc_edu_pos_alexa647_intensity_nucleus_area_sum",
                legend_title = "Sum S phase\n(EdU)")

map_df %>% filter(original_name == "edu positive - Number of Objects")

visualize_model(target_variable = "cc_edu_pos_n_objects",
                legend_title = "Number of\nS-phase Cells\n(EdU)")

map_df %>% filter(original_name == "ROS-back Mean")

visualize_model(target_variable = "vb_ros_back_mean",
                legend_title = "ROS Background\n(Caspase)")

cell_health_variables <- colnames(
    cp_embedding_df %>%
        dplyr::select(starts_with("cc"), starts_with("vb"))
    )

length(cell_health_variables)

pdf_file <- file.path("figures", "repurposing_hub_umaps_consensus.pdf")
pdf(pdf_file, width = 5, height = 5, onefile = TRUE)

for (cell_health_variable in cell_health_variables) {
    umap_gg <- visualize_model(target_variable = cell_health_variable,
                               legend_title = "Prediction:",
                               title = cell_health_variable,
                               dpi = 200,
                               save_figure = FALSE)
}

dev.off()

# Set some plotting defaults
measurement_colors <- c(
    "apoptosis" = "#a6cee3",
    "cell_cycle_arrest" = "#1f78b4",
    "cell_viability" = "#b2df8a",
    "death" = "#33a02c",
    "dna_damage" = "#fb9a99", 
    "g1_arrest" = "#fdbf6f",
    "g2_arrest" = "#ff7f00",
    "g2_m_arrest" = "#005c8c",
    "mitosis" = "green",
    "other" = "black",
    "s_arrest" = "#cab2d6",
    "toxicity" = "#6a3d9a"
)

measurement_labels <- c(
    "apoptosis" = "Apoptosis",
    "cell_cycle_arrest" = "Cell Cycle Arrest",
    "cell_viability" = "Cell Viability",
    "death" = "Death",
    "dna_damage" = "DNA Damage", 
    "g1_arrest" = "G1 Arrest",
    "g2_arrest" = "G2 Arrest",
    "g2_m_arrest" = "G2/M Arrest",
    "mitosis" = "Mitosis",
    "other" = "Other",
    "s_arrest" = "S Arrest",
    "toxicity" = "Toxicity"
)

dye_colors <- c(
    "hoechst" = "#639B94",
    "edu" = "#E45242",
    "gh2ax" = "#E2C552",
    "ph3" = "#7B9C32",
    "hoechst_gh2ax" = "#535f52",
    "hoechst_edu" = "#73414b",
    "edu_gh2ax" = "#e37a48",
    "caspase" = "#F7B1C1",
    "draq" = "#FF6699",
    "draq_caspase" = "#7f4a72",
    "many_cell_cycle" = "#E9DFC3",
    "crispr_efficiency" = "black"
)

dye_labels <- c(
    "hoechst" = "Hoechst",
    "edu" = "EdU",
    "gh2ax" = "gH2AX",
    "ph3" = "pH3",
    "hoechst_gh2ax" = "Hoechst + gH2AX",
    "hoechst_edu" = "Hoechst + EdU",
    "edu_gh2ax" = "EdU + gH2AX",
    "caspase" = "Caspase 3/7",
    "draq" = "DRAQ7",
    "draq_caspase" = "DRAQ7 + Caspase 3/7",
    "many_cell_cycle" = "Cell Cycle (Many Dyes)",
    "crispr_efficiency" = "CRISPR Efficiency"
)

col_types <- readr::cols(
    .default = readr::col_character(),
    shuffle_false = readr::col_double(),
    shuffle_true = readr::col_double()
)

rank_file <- file.path("..", "3.train", "results", "A549_ranked_models.tsv")
model_rank_df <- readr::read_tsv(rank_file, col_types = col_types)

head(model_rank_df, 3)

dmso_embeddings_df <- cp_embedding_df %>%
    dplyr::filter(Metadata_Treatment == "DMSO")

non_dmso_embeddings_df <- cp_embedding_df %>%
    dplyr::filter(Metadata_Treatment != "DMSO")

std_dev_dmso_features <- apply(dmso_embeddings_df %>% dplyr::select(matches("cc_|vb_")), 2, sd)
std_dev_compound_features <- apply(non_dmso_embeddings_df %>% dplyr::select(matches("cc_|vb_")), 2, sd)

std_dev_all_df <- dplyr::bind_cols(as.data.frame(std_dev_dmso_features),
                                   as.data.frame(std_dev_compound_features)) %>%
    dplyr::mutate(features = colnames(dmso_embeddings_df %>% dplyr::select(matches("cc_|vb_")))) %>%
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

ggplot(std_dev_good_df, aes(x = std_dev_dmso_features, y = std_dev_compound_features)) +
    geom_point(aes(color = assay, size = performance_scaled),
               alpha = 0.8) +
    geom_point(data = bad_performing,
               aes(color = assay),
               size = 1,
               alpha = 0.2) +
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

output_file <- file.path("figures", "dmso_vs_compound_standard_deviation.png")
ggsave(output_file, height = 5, width = 6, dpi = 500)
