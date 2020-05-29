suppressPackageStartupMessages(library(ggplot2))

visualize_umap <- function(
    df,
    target_variable,
    legend_title,
    output_dir = "none",
    title = "none",
    dpi = 500,
    save_figure = FALSE,
    print_figure = TRUE
) {

    plot_gg <- ggplot(df, aes(x = umap_x, y = umap_y)) +
        geom_point(aes_string(color = target_variable),
                   size = 0.5,
                   pch = 16,
                   alpha = 0.6) +
        theme_bw() +
        theme(legend.text = element_text(size = 8)) +
        scale_color_viridis_c(name = legend_title) +
        xlab("UMAP 1") +
        ylab("UMAP 2")

    if (title != "none") {
        plot_gg <- plot_gg + ggtitle(title)
    }
    if (save_figure) {
        output_file <- file.path(
            output_dir,
            paste0(
                "umap_repurposing_cell_painting_",
                target_variable,
                "_consensus.png"
            )
        )
        ggsave(output_file, height = 5, width = 6, dpi = dpi)
    }

    if (print_figure) {
      print(plot_gg)
    }
    return(plot_gg)
}
