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
    dplyr::bind_cols(cell_health_df)

print(dim(cp_embedding_df))
head(cp_embedding_df, 3)

ggplot(cp_embedding_df,
       aes(x = umap_x, y = umap_y)) +
    geom_point(aes(color = log10(Metadata_mmoles_per_liter)),
               size = 0.5,
               pch = 16,
               alpha = 0.6) +
    theme_bw() +
    scale_color_viridis_c(name = "log(dose)") +
    xlab("UMAP 1") +
    ylab("UMAP 2")

output_file <- file.path("figures", "umap_repurposing_cell_painting_dose_consensus.png")
ggsave(output_file, height = 5, width = 6, dpi = 500)

visualize_model <- function(target_variable, title, save_figure = TRUE) {
    plot_gg <- ggplot(cp_embedding_df, aes(x = umap_x, y = umap_y)) +
        geom_point(aes_string(color = target_variable),
                   size = 0.5,
                   pch = 16,
                   alpha = 0.6) +
        theme_bw() +
        scale_color_viridis_c(name = title) +
        xlab("UMAP 1") +
        ylab("UMAP 2")
    
    if (save_figure) {
        output_file <- file.path("figures",
                                 paste0("umap_repurposing_cell_painting_",
                                        target_variable,
                                        "_consensus.png"))
        ggsave(output_file, height = 5, width = 6, dpi = 500)
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

visualize_model(target_variable = "log(Metadata_mmoles_per_liter)",
                title = "Compound\nDose (log)") 

map_df %>% filter(original_name == "# Live Cells")

visualize_model(target_variable = "vb_num_live_cells",
                title = "Num Live Cells\n(DRAQ7)")

map_df %>% filter(original_name == "Live Width:Length")

visualize_model(target_variable = "vb_live_cell_width_length",
                title = "Live Cell\n(Width:Length)\n(DRAQ7)")

map_df %>% filter(original_name == "Live Cell Roundness")

visualize_model(target_variable = "vb_live_cell_roundness",
                title = "Live Cell Roundness\n(DRAQ7)")

map_df %>% filter(original_name == "ALL - Number of Objects")

visualize_model(target_variable = "cc_all_n_objects",
                title = "Number of Objects\n(Hoechst)")

map_df %>% filter(original_name == "Live Cell Area")

visualize_model(target_variable = "vb_live_cell_area",
                title = "Live Cell Area\n(DRAQ7)")

map_df %>% filter(original_name == "CC - Number of Objects")

visualize_model(target_variable = "cc_cc_n_objects",
                title = "Num Cell\nCycle Objects\n(Hoechst)")

map_df %>% filter(original_name == "G1 - Number of Objects")

visualize_model(target_variable = "cc_g1_n_objects",
                title = "G1 Objects\n(Many Dyes)")

map_df %>% filter(original_name == "edu positive - Intensity Nucleus Alexa 647 Sum - Sum per Well")

visualize_model(target_variable = "cc_edu_pos_alexa647_intensity_nucleus_area_sum",
                title = "Sum S phase\n(EdU)")

map_df %>% filter(original_name == "edu positive - Number of Objects")

visualize_model(target_variable = "cc_edu_pos_n_objects",
                title = "Number of\nS-phase Cells\n(EdU)")

map_df %>% filter(original_name == "ROS-back Mean")

visualize_model(target_variable = "vb_ros_back_mean",
                title = "ROS Background\n(Caspase)")

cell_health_variables <- colnames(
    cp_embedding_df %>%
        dplyr::select(starts_with("cc"), starts_with("vb"))
    )

length(cell_health_variables)

pdf_file <- file.path("figures", "repurposing_hub_umaps_consensus.pdf")
pdf(pdf_file, width = 5, height = 5, onefile = TRUE)

for (cell_health_variable in cell_health_variables) {
    umap_gg <- visualize_model(target_variable = cell_health_variable,
                               title = "Prediction",
                               save_figure = FALSE)
    
    print(umap_gg)
}

dev.off()
