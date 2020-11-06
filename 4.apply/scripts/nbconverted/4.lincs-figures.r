suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source(file.path("scripts", "visualize_utils.R"))
source(file.path("repurposing_cellhealth_shiny", "util.R"))
source(file.path("repurposing_cellhealth_shiny", "dose_utils.R"))

# Define a function to test enrichment of a particular moa in a model
get_enrichment <- function(df, moa, model, hypothesis = "greater") {
    midway_point <- quantile(df %>% dplyr::pull(!!model), 0.5)
    midway_point <- as.numeric(paste(midway_point))
    
    cp_greater_df <- df %>%
        dplyr::filter((df %>% dplyr::pull(!!model)) > midway_point)
    cp_less_df <- df %>%
        dplyr::filter((df %>% dplyr::pull(!!model)) <= midway_point)
    
    num_moa_high <- dim(cp_greater_df %>% dplyr::filter(moa == !!moa))[1]
    num_moa_low <- dim(cp_less_df %>% dplyr::filter(moa == !!moa))[1]

    num_notmoa_high <- dim(cp_greater_df %>% dplyr::filter(moa != !!moa))[1]
    num_notmoa_low <- dim(cp_less_df %>% dplyr::filter(moa != !!moa))[1]

    contingency_mat <- matrix(
        c(
            c(num_moa_high, num_moa_low),
            c(num_notmoa_high, num_notmoa_low)
        ),
        nrow = 2
    )
    
    result <- fisher.test(contingency_mat, alternative = hypothesis)
    return(result)
}

consensus = "modz"

# Load data
data <- load_data(path="repurposing_cellhealth_shiny/")
cp_embedding_df <- data[["moa"]]

cp_embedding_df <- cp_embedding_df %>%
    dplyr::mutate(Metadata_Treatment = "Compound")

cp_embedding_df$Metadata_Treatment[cp_embedding_df$Metadata_broad_core_id == "DMSO"] = "DMSO"

print(dim(cp_embedding_df))
head(cp_embedding_df, 3)

# Load Dose and Rank information
dose_df <- data[["dose"]]
rank_df <- data[["rank"]]

model_dict_df <- rank_df %>%
    dplyr::select(target, readable_name)
model_dict_df$target <- paste(model_dict_df$target)
model_dict_df$readable_name <- paste(model_dict_df$readable_name)

# Reshape the moa dataframe for different variable plotting
melt_id_vars <- c(
  "Metadata_Plate_Map_Name",
  "Metadata_pert_well",
  "Metadata_broad_core_id",
  "Metadata_broad_sample",
  "Metadata_dose_recode",
  "Metadata_mmoles_per_liter",
  "umap_x",
  "umap_y",
  "broad_id",
  "pert_iname",
  "InChIKey14",
  "moa",
  "target",
  "clinical_phase",
  "alternative_moa",
  "alternative_target",
  "broad_date"
)

moa_long_df <- reshape2::melt(
  cp_embedding_df,
  id.vars = melt_id_vars,
  value.name = "model_score",
  variable.name = "model"
)

moa_long_df$model <- paste(moa_long_df$model)
moa_long_df <- moa_long_df %>%
  dplyr::left_join(model_dict_df, by = c("model" = "target"))

moa_long_df$model_score <- as.numeric(paste(moa_long_df$model_score))

head(moa_long_df)

# Load True Cell Counts
batch <- "2016_04_01_a549_48hr_batch1"
lincs_commit <- "e6852b49992b4fa2f0d36cc07011856ebee6326a"
raw_url <- paste0("https://raw.githubusercontent.com/broadinstitute/lincs-cell-painting/", lincs_commit)
barcode_platemap_file <- paste0(raw_url, "/metadata/platemaps/", batch, "/barcode_platemap.csv")
barcode_platemap_df <- readr::read_csv(barcode_platemap_file, col_types = readr::cols())

plates <- unique(barcode_platemap_df$Assay_Plate_Barcode)
platemaps <- unique(barcode_platemap_df$Plate_Map_Name)
print(length(plates))
print(length(platemaps))

# Load platemaps
platemap_data <- list()
for (platemap in platemaps) {
    platemap_file <- paste0(
        raw_url, "/metadata/platemaps/", batch, "/", "platemap", "/", platemap, ".txt"
    )
    
    platemap_df <- readr::read_tsv(platemap_file, col_types = readr::cols())
    platemap_df$broad_sample <- tidyr::replace_na(platemap_df$broad_sample, "DMSO")
    platemap_df$mg_per_ml <- tidyr::replace_na(platemap_df$mg_per_ml, 0)
    platemap_df$mmoles_per_liter <- tidyr::replace_na(platemap_df$mmoles_per_liter, 0)
    
    platemap_data[[platemap]] <- platemap_df
}

# Load Cell Count Files
cell_count_data = list()
for (plate in plates) {
    platemap_name <- barcode_platemap_df %>%
        dplyr::filter(Assay_Plate_Barcode == !!plate) %>% dplyr::pull(Plate_Map_Name)
    
    platemap_df <- platemap_data[[platemap_name]]

    cell_count_file <- paste0(
        raw_url, "/profiles/cell_count/", batch, "/", plate, "/", plate, "_cell_count.csv"
    )
    
    cell_count_df <- tryCatch({
        readr::read_csv(cell_count_file, col_types = readr::cols()) %>%
            dplyr::left_join(platemap_df, by = c("Image_Metadata_Well" = "well_position"))
        }, error = function(w) {
            return("not found")
    })
    
    if (cell_count_df == "not found") {
        next
    } else {
        cell_count_data[[plate]] <- cell_count_df
    }
}

# Combine and process cell count files
all_count_df <- do.call(rbind, cell_count_data) %>%
    dplyr::group_by(plate_map_name, Image_Metadata_Well, broad_sample, mmoles_per_liter) %>%
    dplyr::mutate(mean_cell_count = mean(cell_count)) %>%
    dplyr::select(plate_map_name, Image_Metadata_Well, broad_sample, mmoles_per_liter, mean_cell_count) %>%
    dplyr::distinct() %>%
    dplyr::mutate(dose_merge = round(mmoles_per_liter, 4))

cell_count_info <- all_count_df %>%
    dplyr::right_join(cp_embedding_df %>%
        dplyr::mutate(dose_merge = round(Metadata_mmoles_per_liter, 4)),
                 by = c("plate_map_name" = "Metadata_Plate_Map_Name",
                        "Image_Metadata_Well" = "Metadata_pert_well",
                        "dose_merge" = "dose_merge"))

head(cell_count_info, 3)

# Load Viability Validation
validation_file = file.path(
    "..", "5.validate-repurposing", "results", "depmap_viability_validation.tsv.gz"
)

validation_df <- readr::read_tsv(validation_file, col_types = readr::cols())

print(dim(validation_df))
head(validation_df, 3)

# Subset controls and compounds
cp_embedding_df$Metadata_Controls <- "Compounds"
cp_embedding_df$Metadata_Controls[cp_embedding_df$Metadata_broad_sample == "DMSO"] <- "DMSO"
cp_embedding_df$Metadata_Controls[cp_embedding_df$pert_iname == "bortezomib"] <- "bortezomib"
cp_embedding_df$Metadata_Controls[cp_embedding_df$pert_iname == "MG-132"] <- "MG-132"

table(cp_embedding_df$Metadata_Controls)

control_df <- cp_embedding_df %>% dplyr::filter(Metadata_Controls != "Compounds")
control_df$Metadata_Controls <- factor(
    control_df$Metadata_Controls, levels = c("DMSO", "bortezomib", "MG-132")
)
compound_df <- cp_embedding_df %>%
    dplyr::filter(Metadata_Controls == "Compounds")

figure_theme <- theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
)

# Panel A - Dose Diferences
panel_a_gg = ggplot(
    validation_df, aes(y = dose, x = Metadata_mmoles_per_liter)
    ) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.6, lwd = 0.4, linetype = "dashed", color = "red") +
    geom_point(size = 0.3, alpha = 1) +
    theme_bw() +
    ylab("PRISM Dose\n(micromoles per liter)") +
    xlab("Drug Repurposing Dose\n(micromoles per liter)") +
    coord_fixed()

# Generate text for panel B
# Calculate correlation between orthogonal tests
spearman_result <- cor.test(
    validation_df$cell_health_viability,
    validation_df$depmap_viability,
    method = "spearman")

pval <- spearman_result$p.value
if (pval < 0.001) {
    print(paste0("Calculated p value: ", pval))
    # For presentation purposes:
    pval <- "p < 10e-3"
} else {
    pval <- paste0("p = ", pval)
}
stat <- round(as.numeric(paste(spearman_result$estimate)), 3)

result_text = paste0("Spearman = ", stat, "\n", pval)

result_text

# Panel B - Dose Correlation
panel_b_gg = ggplot(validation_df,
                    aes(x = cell_health_viability, y = depmap_viability)) +
    geom_point(size = 0.5, alpha = 0.3) +
    theme_bw() +
    geom_smooth(method = "lm", formula = y~x) +
    annotate("text", label = result_text, x = -2, y = -8.5, size = 2.5) +
    xlab("Cell Health Model Predictions\n(Number of Live Cells)") +
    ylab("PRISM Assay\n(A549 Phycoerythrin Intensity)") +
    coord_fixed() +
    theme()

# Subset data for panel c
focus_moa <- "PLK inhibitor"

moa_df <- compound_df %>% dplyr::filter(moa == !!focus_moa)
print(dim(moa_df))

# What are the existing PLK inhibitors
table(moa_df$pert_iname)

# Define enrichment of select MOAs in select models
moa <- "proteasome inhibitor"
model <- "vb_ros_mean"

get_enrichment(cp_embedding_df, moa, model)

# Define enrichment of select MOAs in select models
moa <- "PLK inhibitor"
model <- "cc_g1_n_objects"

get_enrichment(cp_embedding_df, moa, model, hypothesis = "less")

# Panel C
panel_c_gg <- ggplot(compound_df,
       aes(
           x = cc_g1_n_objects,
           y = vb_ros_mean,
           color = Metadata_dose_recode
       )) +
    geom_point(
        aes(shape = Metadata_Controls),
        data = compound_df %>% dplyr::filter(moa != "PLK inhibitor"),
        size = 0.5,
        alpha = 0.6
    ) +
    geom_point(
        data = control_df,
        aes(shape = Metadata_Controls),
        fill = "red",
        color = "black",
        size = 1,
        alpha = 0.4,
    ) +
    geom_point(
        data = moa_df,
        aes(fill = Metadata_dose_recode,
            shape = moa),
        color = "black",
        size = 2,
        alpha = 0.7
    ) +
    theme_bw() +
    ylab("Predicted ROS") +
    xlab("Predicted G1 Cell Count") +
    scale_shape_manual(
        name = "Point Labels",
        values = c("DMSO" = 21,
                   "bortezomib" = 23,
                   "MG-132" = 25,
                   "PLK inhibitor" = 24,
                   "Compounds" = 16),
        labels = c("DMSO" = "DMSO",
                   "bortezomib" = "Bortezomib",
                   "MG-132" = "MG-132",
                   "PLK inhibitor" = "PLK inhibitor",
                   "Compounds" = "Other")
      ) +
    scale_color_continuous(name = "Dose Rank") +
    scale_fill_continuous(name = "Dose Rank") +
    guides(shape = guide_legend(order = 1,
                                override.aes = list(size = c(2, 1, 2, 2, 3.5))))

# Generate Figure 4C Figure Legend
dose_rank_fill_legend <- cowplot::get_legend(
    ggplot(compound_df, aes(x = cc_g1_n_objects, y = vb_ros_mean, color = Metadata_dose_recode)) +
    geom_point(data = compound_df,
               size = 0.5,
               alpha = 0.6) +
    scale_color_continuous(name = "Dose Rank") +
    figure_theme +
    theme(legend.key.size = unit(0.5, "lines"))
)

control_point_label_legend <- cowplot::get_legend(
    ggplot(compound_df,
           aes(x = cc_g1_n_objects,
               y = vb_ros_mean,
               color = Metadata_dose_recode)) +
    geom_point(data = control_df,
               aes(
                   shape = Metadata_Controls),
                   fill = "red",
                   color = "black",
                   size = 1,
                   alpha = 0.4
              ) +
    scale_shape_manual(name = "Controls",
                       values = c("DMSO" = 21,
                                  "bortezomib" = 23,
                                  "MG-132" = 25),
                       labels = c("DMSO" = "DMSO",
                                  "bortezomib" = "Bortezomib",
                                  "MG-132" = "MG-132")) +
    theme_bw() +
    figure_theme +
    theme(legend.key.size = unit(0.4, "lines"))
    )

compound_point_label_legend <- cowplot::get_legend(
    ggplot(compound_df,
           aes(x = cc_g1_n_objects,
               y = vb_ros_mean)) +
    geom_point(data = compound_df,
               aes(shape = Metadata_Controls), size = 0.5) +
    geom_point(data = moa_df,
               aes(shape = moa),
               color = "black",
               size = 2,
               alpha = 0.7) +
    scale_shape_manual(name = "Compounds",
                       values = c("PLK inhibitor" = 24,
                                  "Compounds" = 16),
                       labels = c("PLK inhibitor" = "PLK inhibitor",
                                  "Compounds" = "Other")) +
    theme_bw() +
    guides(shape = guide_legend(
        order = 1,
        override.aes = list(size = c(0.5, 2))
    ))  +
    figure_theme +
    theme(legend.key.size = unit(0.4, "lines"))
)

# Put everything together into a single figure
panel_c_with_legend <- cowplot::ggdraw(
    panel_c_gg +
    theme(legend.position = "none") +
    figure_theme +
    cowplot::draw_plot(
        control_point_label_legend, x = 0.8, y = 50
        ) +
    cowplot::draw_plot(
        compound_point_label_legend, x = 0.8, y = 37
        ) +
    cowplot::draw_plot(
        dose_rank_fill_legend, x = 0.65, y = 22
        ) 
)

# Panel D - Dose Response
panel_d_gg <- suppressWarnings(
    get_dose_curve(
        moa_long_df,
        dose_df,
        model = "cc_g1_n_objects",
        pert_name = "HMN-214",
        cell_health_model = "G1 Cell Count"
    )
) + ylab("Predicted\nG1 Cell Count")

# Prep data for plotting panel e
dna_damage_g1_density_df <- cp_embedding_df %>%
    dplyr::mutate(density_class = "Other Compound")

aurorakinase_inhibit <- grepl("Aurora", cp_embedding_df$moa)
tubulin_inhibit <- grepl("tubulin", cp_embedding_df$moa)

dna_damage_g1_density_df[aurorakinase_inhibit, "density_class"] <- "Aurora Kinase Inhibitor"
dna_damage_g1_density_df[tubulin_inhibit, "density_class"] <- "Tubulin Inhibitor"
dna_damage_g1_density_df[tubulin_inhibit, "density_class"] <- "Tubulin Inhibitor"
dna_damage_g1_density_df$density_class[dna_damage_g1_density_df$Metadata_broad_sample == "DMSO"] <- "DMSO"
dna_damage_g1_density_df$density_class[dna_damage_g1_density_df$pert_iname == "bortezomib"] <- "Proteasome Inhibitor"
dna_damage_g1_density_df$density_class[dna_damage_g1_density_df$pert_iname == "MG-132"] <- "Proteasome Inhibitor"

dna_damage_g1_density_df$density_class <- factor(
    dna_damage_g1_density_df$density_class,
    levels = c("DMSO", "Proteasome Inhibitor", "Other Compound", "Tubulin Inhibitor", "Aurora Kinase Inhibitor")
)

dna_damage_labels <- c(
    "DMSO" = "DMSO",
    "Proteasome Inhibitor"= "Proteasome Inhibitor",
    "Other Compound"= "Other Compound",
    "Tubulin Inhibitor"= "Tubulin Inhibitor",
    "Aurora Kinase Inhibitor"= "Aurora Kinase Inhibitor"
)

dna_damage_colors <- c(
    "DMSO" = "orange",
    "Proteasome Inhibitor"= "red",
    "Other Compound"= "yellow",
    "Tubulin Inhibitor"= "#8C72C0",
    "Aurora Kinase Inhibitor"= "#1C2D40"
)

# Define enrichment of tubulin and aurora kinase inhibitors
# Cannot easily use function since we're testing if two MOAs are enriched
model <- "cc_g1_n_spots_h2ax_mean"

midway_point <- quantile(dna_damage_g1_density_df %>% dplyr::pull(!!model), 0.5)
midway_point <- as.numeric(paste(midway_point))

cp_greater_df <- dna_damage_g1_density_df %>%
    dplyr::filter((dna_damage_g1_density_df %>% dplyr::pull(!!model)) > midway_point)
cp_less_df <- dna_damage_g1_density_df %>%
    dplyr::filter((dna_damage_g1_density_df %>% dplyr::pull(!!model)) <= midway_point)

# Build contingency matrix
test_moa <- c("Tubulin Inhibitor", "Aurora Kinase Inhibitor")
num_moa_high <- dim(cp_greater_df %>% dplyr::filter(density_class %in% test_moa))[1]
num_moa_low <- dim(cp_less_df %>% dplyr::filter(density_class %in% test_moa))[1]

num_notmoa_high <- dim(cp_greater_df %>% dplyr::filter(!(density_class %in% test_moa)))[1]
num_notmoa_low <- dim(cp_less_df %>% dplyr::filter(!(density_class %in% test_moa)))[1]

contingency_mat <- matrix(
    c(
        c(num_moa_high, num_moa_low),
        c(num_notmoa_high, num_notmoa_low)
    ),
    nrow = 2
)

result <- fisher.test(contingency_mat, alternative = "greater")
result

# Panel E - DNA Damage in G1 Cell Count
panel_e_gg <- ggplot(dna_damage_g1_density_df, aes_string(x = model)) +
    geom_density(aes(fill = density_class), alpha = 0.5) +
    theme_bw() +
    figure_theme +
    theme(legend.position = c(0.75, 0.8), legend.key.size = unit(0.7, "lines")) +
    geom_vline(xintercept = midway_point, color = "red", linetype = "dashed") +
    scale_fill_manual(name = "Compound class",
                      labels = dna_damage_labels,
                      values = dna_damage_colors) +
    xlab("Predicted G1 - # of gH2AX Spots") +
    ylab("Density")

# Panel F - Dose Response for barasertib
panel_f_gg <- suppressWarnings(
    get_dose_curve(
        moa_long_df,
        dose_df,
        model = "cc_g1_n_spots_h2ax_mean",
        pert_name = "barasertib",
        cell_health_model = "G1 - # of gH2AX Spots"
    )
) + ylab("Predicted\nG1 - # of gH2AX Spots")

# Compile figure panels by row and column
panel_a_b <- cowplot::plot_grid(
    cowplot::ggdraw(),
    panel_a_gg + figure_theme + theme(plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm")),
    cowplot::ggdraw(),
    panel_b_gg + figure_theme + theme(plot.margin = unit(c(1, -0.5, 0.5, 0.5), "cm")),
    labels = c("", "a", "", "b"),
    align = "hv",
    ncol = 4,
    axis = "l",
    rel_widths = c(0.01, 1, 0.01, 1.1),
    hjust = c(0, -8, 0, -8)
)

panel_c_d <- cowplot::plot_grid(
    panel_c_with_legend,
    panel_d_gg + theme(axis.text.x = element_text(angle = 90)) + figure_theme,
    nrow = 2,
    align = "v",
    axis = "l",
    rel_heights = c(1, 0.7),
    labels = c("c", "d")
)

panel_e_f <- (
    cowplot::plot_grid(
        panel_e_gg + figure_theme,
        panel_f_gg + figure_theme + theme(axis.text.x = element_text(angle = 90)),
        nrow = 2,
        align = "v",
        axis = "l",
        rel_heights = c(1, 0.7),
        labels = c("e", "f")
    )
)

figure_4_gg <- cowplot::plot_grid(
    panel_a_b,
    cowplot::plot_grid(
        panel_c_d,
        panel_e_f,
        labels = c("", ""),
        axis = "l",
        align = "hv",
        rel_widths = c(1, 1),
        ncol = 2
    ),
    labels = c("", ""),
    align = "hv",
    axis = "b",
    nrow = 2,
    rel_heights = c(0.4, 1)
)

figure_file <- file.path("figures", "lincs_main_figure_4.png")

cowplot::save_plot(
    filename = figure_file,
    plot = figure_4_gg,
    dpi = 500,
    base_width = 8,
    base_height = 8
)

figure_4_gg

# Panel A
legend_title <- "Predicted\nG1 Cell\nCount"
sup_panel_a_gg <- visualize_umap(
    df = cp_embedding_df,
    target_variable = "cc_g1_n_objects",
    legend_title = legend_title,
    print_figure = FALSE
)

# Tweak the color scale to emphasize gradient
sup_panel_a_gg <- sup_panel_a_gg +
    scale_color_viridis_c(name = legend_title,
        values = scales::rescale(c(1, 0.8, 0.2)))

# Panel B
legend_title <- "Predicted\nROS"
sup_panel_b_gg <- visualize_umap(
    df = cp_embedding_df,
    target_variable = "vb_ros_mean",
    legend_title = legend_title,
    print_figure = FALSE
)

# Tweak the color scale to emphasize gradient
sup_panel_b_gg <- sup_panel_b_gg +
    scale_color_viridis_c(name = legend_title,
        values = scales::rescale(c(1, 0.96, 0.95)))

# Panel C
sup_panel_c_gg <- ggplot(cell_count_info, aes(x = umap_x, y = umap_y)) +
    xlab("UMAP X") +
    ylab("UMAP Y") +
    geom_point(aes(color = mean_cell_count),
               size = 0.65,
               pch = 16,
               alpha = 0.6) +
    geom_point(data = control_df,
              aes(shape = Metadata_Controls),
              fill = "grey",
              color = "black",
              size = 1.5,
              alpha = 0.4,
              show.legend = TRUE) +
    scale_shape_manual(
      name = "Controls",
      values = c("DMSO" = 21, "bortezomib" = 23, "MG-132" = 25),
      labels = c("DMSO" = "DMSO", "bortezomib" = "Bortezomib", "MG-132" = "MG-132")
    ) +
    theme_bw() +
    scale_color_viridis_c(name = "Mean\nGround Truth\nCell Count") +
    guides(shape = guide_legend(order = 1))

figure_file_ext <- file.path("figures", "lincs_umap_before_annotation")

for (extension in c(".png", ".pdf")) {
    figure_file <- paste0(figure_file_ext, extension)
    ggsave(figure_file, sup_panel_c_gg, dpi = 500, width = 10, height = 9)
}

# Create multiplot
panel_ab <- cowplot::plot_grid(
    sup_panel_a_gg + figure_theme,
    sup_panel_b_gg + figure_theme,
    labels = c("a", "b"),
    align = "hv",
    vjust = c(1, 1),
    nrow = 2
)

sup_fig <- cowplot::plot_grid(
    panel_ab,
    sup_panel_c_gg + figure_theme,
    labels = c("", "c"),
    align = "hv",
    rel_widths = c(0.6, 1)
)

figure_file <- file.path("figures", "lincs_supplementary_figure.png")

cowplot::save_plot(
    filename = figure_file,
    plot = sup_fig,
    dpi = 500,
    base_width = 8,
    base_height = 5
)

sup_fig
