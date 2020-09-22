suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

source(file.path("scripts", "assay_themes.R"))

consensus <- "modz"
results_dir <- file.path("results")

output_figure_file <- file.path(
    "figures", "summary", "crispr_umap_supplementary_figure.png"
)

point_alpha <- 0.7
point_size <- 0.8

umap_file <- file.path(
    results_dir, paste0("profile_umap_with_cell_health_", consensus, ".tsv")
)
umap_df <- readr::read_tsv(umap_file, col_types = readr::cols())

print(dim(umap_df))
head(umap_df, 3)

cell_line_gg <- ggplot(umap_df, aes(x = umap_x, y = umap_y)) +
    geom_point(aes(color = Metadata_cell_line), alpha = point_alpha, size = point_size) +
    xlab("UMAP X") +
    ylab("UMAP Y") +
    scale_color_manual(
        name = "Cell Line",
        values = cell_line_colors,
        labels = cell_line_labels
    ) +
    theme_bw()

cell_line_gg

g1_count_gg <- ggplot(umap_df, aes(x = umap_x, y = umap_y)) +
    geom_point(aes(color = cc_g1_n_objects), alpha = point_alpha, size = point_size) +
    xlab("UMAP X") +
    ylab("UMAP Y") +
    theme_bw() +
    scale_color_viridis_c(
        name = "G1 Cell\nCount",
        values = scales::rescale(c(1, 0.8, 0.2))
    )

g1_count_gg

ros_gg <- ggplot(umap_df, aes(x = umap_x, y = umap_y)) +
    geom_point(aes(color = vb_ros_mean), alpha = point_alpha, size = point_size) +
    xlab("UMAP X") +
    ylab("UMAP Y") +
    theme_bw() +
    scale_color_viridis_c(
        name = "Reactive\nOxygen\nSpecies",
        values = scales::rescale(c(10, 3.5, 2))
    )

ros_gg

scatter_gg <- ggplot(umap_df, aes(x = cc_g1_n_objects, y = vb_ros_mean)) +
    geom_point(aes(color = Metadata_cell_line), alpha = point_alpha, size = point_size) +
    xlab("G1 Cell Count") +
    ylab("Reactive Oxygen Species") +
    theme_bw() +
    scale_color_manual(
        name = "Cell Line",
        values = cell_line_colors,
        labels = cell_line_labels
    )

scatter_gg

# Create multiplot
full_figure <- cowplot::plot_grid(
    cell_line_gg,
    g1_count_gg,
    ros_gg,
    scatter_gg,
    labels = c("a", "b", "c", "d"),
    align = "hv",
    nrow = 2
)

cowplot::save_plot(
    filename = output_figure_file,
    plot = full_figure,
    dpi = 500,
    base_width = 7,
    base_height = 5
)

full_figure
