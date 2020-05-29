suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

source(file.path("scripts", "visualize_utils.R"))
source(file.path("repurposing_cellhealth_shiny", "util.R"))
source(file.path("repurposing_cellhealth_shiny", "dose_utils.R"))

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

# Panel A
legend_title <- "Predicted\nG1 Cell\nCount"
panel_a_gg <- visualize_umap(
    df = cp_embedding_df,
    target_variable = "cc_g1_n_objects",
    legend_title = legend_title,
    print_figure = FALSE
)

# Tweak the color scale to emphasize gradient
panel_a_gg <- panel_a_gg +
    scale_color_viridis_c(name = legend_title,
        values = scales::rescale(c(1, 0.8, 0.2)))

# Panel B
legend_title <- "Predicted\nROS"
panel_b_gg <- visualize_umap(
    df = cp_embedding_df,
    target_variable = "vb_ros_mean",
    legend_title = legend_title,
    print_figure = FALSE
)

# Tweak the color scale to emphasize gradient
panel_b_gg <- panel_b_gg +
    scale_color_viridis_c(name = legend_title,
        values = scales::rescale(c(1, 0.96, 0.95)))

cp_embedding_df$Metadata_Controls <- "Compounds"
cp_embedding_df$Metadata_Controls[cp_embedding_df$Metadata_broad_sample == "DMSO"] <- "DMSO"
cp_embedding_df$Metadata_Controls[cp_embedding_df$pert_iname == "bortezomib"] <- "bortezomib"
cp_embedding_df$Metadata_Controls[cp_embedding_df$pert_iname == "MG-132"] <- "MG-132"

table(cp_embedding_df$Metadata_Controls)

focus_moa <- "PLK inhibitor"

control_df <- cp_embedding_df %>% dplyr::filter(Metadata_Controls != "Compounds")
compound_df <- cp_embedding_df %>%
    dplyr::filter(Metadata_Controls == "Compounds")

moa_df <- compound_df %>% dplyr::filter(moa == !!focus_moa)

print(dim(moa_df))

# What are the existing PLK inhibitors
table(moa_df$pert_iname)

# Panel C
panel_c_gg <- ggplot(compound_df,
       aes(
           x = cc_g1_n_objects,
           y = vb_ros_mean,
           color = Metadata_dose_recode
       )) +
    geom_point(
        aes(shape = Metadata_Controls),
        data = compound_df,
        size = 0.5,
        alpha = 0.6
    ) +
    geom_point(
        data = control_df,
        aes(shape = Metadata_Controls),
        fill = "grey",
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

figure_theme = theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7)
)

# Create multiplot
panel_a_b <- cowplot::plot_grid(
    panel_a_gg + figure_theme,
    panel_b_gg + figure_theme,
    labels = c("a", "b"),
    align = "hv",
    nrow = 2
)

panel_c_d <- cowplot::plot_grid(
    panel_c_gg + figure_theme +
        theme(legend.position = "none", axis.title.x = element_text(vjust = 5)),
    panel_d_gg + figure_theme + theme(axis.text.x = element_text(angle = 90)),
    labels = c("c", "d"),
    vjust = c(1.5, 0.7),
    align = "hv",
    nrow = 2,
    rel_heights = c(1, 0.5)
)

panel_c_d_legend <- cowplot::plot_grid(
    cowplot::get_legend(panel_c_gg + figure_theme),
    cowplot::ggdraw(),
    labels = c("", ""),
    align = "hv",
    nrow = 2,
    rel_heights = c(1, 0.85)
)

figure_4_gg <- cowplot::plot_grid(
    panel_a_b,
    panel_c_d,
    panel_c_d_legend,
    labels = c("", "", ""),
    align = "v",
    ncol = 3,
    rel_widths = c(1, 0.8, 0.3)
)

figure_file <- file.path("figures", "lincs_main_figure_4.png")

cowplot::save_plot(
    filename = figure_file,
    plot = figure_4_gg,
    dpi = 500,
    base_width = 8,
    base_height = 6
)

figure_4_gg

# Load Viability Validation
validation_file = file.path(
    "..", "5.validate-repurposing", "results", "depmap_viability_validation.tsv.gz"
)

validation_df <- readr::read_tsv(validation_file, col_types = readr::cols())

print(dim(validation_df))
head(validation_df, 3)

model <- "vb_num_live_cells"
sup_panel_a_gg <- ggplot(cp_embedding_df, aes(x = umap_x, y = umap_y)) +
    xlab("UMAP X") +
    ylab("UMAP Y") +
    geom_point(aes_string(color = model),
               size = 1.25,
               pch = 16,
               alpha = 0.6) +
    geom_point(data = control_df,
              aes(shape = Metadata_Controls),
              fill = "grey",
              color = "black",
              size = 2,
              alpha = 0.4,
              show.legend = TRUE) +
    scale_shape_manual(
      name = "Controls",
      values = c("DMSO" = 21, "bortezomib" = 23, "MG-132" = 25),
      labels = c("DMSO" = "DMSO", "bortezomib" = "Bortezomib", "MG-132" = "MG-132")
    ) +
    theme_bw() +
    scale_color_viridis_c(name = "Predicted\n# Live Cells",
                          values = scales::rescale(c(1, 0.8, 0.1)))

# Panel B - Dose Diferences
sup_panel_b_gg = ggplot(validation_df,
                        aes(x = dose, y = Metadata_mmoles_per_liter)) +
    geom_point(size = 1, alpha = 0.6) +
    theme_bw() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    xlab("PRISM Dose\n(mmoles per liter)") +
    ylab("Dependency Map Dose\n(mmoles per liter)") +
    coord_fixed()

# Calculate correlation between orthogonal tests
spearman_result <- cor.test(
    validation_df$cell_health_viability,
    validation_df$depmap_viability,
    method = "spearman")

pval <- format(spearman_result$p.value, digits = 3)
stat <- round(as.numeric(paste(spearman_result$estimate)), 3)

result_text = paste0("Spearman = ", stat, "\np = ", pval)

result_text

panel_c_gg = ggplot(validation_df,
                    aes(x = cell_health_viability, y = depmap_viability)) +
    geom_point(size = 0.5, alpha = 0.3) +
    theme_bw() +
    geom_smooth(method = "lm", formula = y~x) +
    annotate("text", label = result_text, x = -2, y = -8.5, size = 2.5) +
    xlab("Cell Health Model Predictions\n(Number of Live Cells)") +
    ylab("Dependency Map Viability") +
    coord_fixed() +
    theme()

# Create multiplot
panel_b_c <- cowplot::plot_grid(
    sup_panel_b_gg + figure_theme,
    panel_c_gg + figure_theme,
    labels = c("b", "c"),
    align = "hv",
    vjust = c(1, 1),
    nrow = 2,
    rel_heights = c(0.6, 1)
)

sup_fig <- cowplot::plot_grid(
    sup_panel_a_gg + figure_theme,
    panel_b_c,
    labels = c("a", ""),
    align = "hv",
    rel_widths = c(1, 0.6)
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
