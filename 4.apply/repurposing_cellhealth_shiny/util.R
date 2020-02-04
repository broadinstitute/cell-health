# Functions to enable shiny plot functionality

library(dplyr)
library(ggplot2)
library(ggrepel)

get_target <- function(rank, model) {
  target <- rank %>%
    dplyr::filter(target == !!model) %>%
    dplyr::pull(original_name)
  
  return(target)
}

build_cell_health_scatter <- function(
  moa_full_df, compound_df, dmso_df, model, target
  ) {
  
  moa_scatter_gg <- ggplot(moa_full_df,
       aes_string(x = "vb_num_live_cells",
                  y = model)) +
  xlab("# Live Cells") +
  ylab(target) +
  geom_point(aes(color = Metadata_dose_recode),
             size = 1.25,
             pch = 16,
             alpha = 0.6) +
  geom_point(data = compound_df,
             aes(color = Metadata_dose_recode),
             size = 5,
             pch = 17,
             alpha = 0.7,
             show.legend = FALSE) +
  geom_point(data = dmso_df,
             aes(shape = Metadata_broad_sample),
             fill = "grey",
             color = "black",
             size = 3,
             alpha = 0.7,
             show.legend = TRUE) +
  scale_color_viridis_c(name = "Dose\nLevel") +
  scale_shape_manual(name = "", values = 21, labels = "DMSO") +
  theme_bw() +
  guides(
    shape = guide_legend(order = 2,
                         keywidth = 0.1,
                         keyheight = 0.1,
                         title = "",
                         override.aes = list(
                           size = 3,
                           pch = 21,
                           fill = "grey",
                           color = "black",
                           alpha = 0.7)
                         )
    )
  
  return(moa_scatter_gg)
}

build_umap_scatter <- function(
  moa_full_df, compound_df, dmso_df, model, target
  ){

  umap_scatter_gg <- ggplot(moa_full_df, aes(x = umap_x, y = umap_y)) +
    xlab("UMAP X") +
    ylab("UMAP Y") +
    geom_point(aes(color = Metadata_dose_recode),
               size = 1.25,
               pch = 16,
               alpha = 0.6) +
    geom_point(data = compound_df,
               aes(color = Metadata_dose_recode),
               size = 5,
               pch = 17,
               alpha = 0.7,
               show.legend = FALSE) +
    geom_text_repel(data = compound_df,
                    arrow = arrow(length = unit(0.01, "npc")),
                    size = 4,
                    segment.size = 0.5,
                    segment.alpha = 0.9,
                    force = 10,
                    aes(label = paste0("Dose:", Metadata_dose_recode),
                        x = umap_x,
                        y = umap_y)) +
    geom_point(data = dmso_df,
               aes(shape = Metadata_broad_sample),
               fill = "grey",
               color = "black",
               size = 3,
               alpha = 0.7,
               show.legend = TRUE) +
    ggtitle(target) +
    scale_color_viridis_c(name = "Dose\nLevel") +
    scale_shape_manual(name = "", values = 21, labels = "DMSO") +
    theme_bw() +
    guides(
      shape = guide_legend(order = 2,
                           keywidth = 0.1,
                           keyheight = 0.1,
                           title = "",
                           override.aes = list(
                             size = 3,
                             pch = 21,
                             fill = "grey",
                             color = "black",
                             alpha = 0.7))
    )
}


build_rank_plot <- function(rank_df) {
  ggplot(rank_df,
         aes(x = original_name,
             y = shuffle_false)) +
    geom_bar(aes(fill = assay, color = to_highlight),
             lwd = 1.25,
             stat="identity") +
    ylab("Test Set Regression Performance") +
    xlab("") +
    ggtitle("A549 Cell Line") +
    scale_color_manual(name = "",
                       values = c("black" = "black",
                                  "white" = "white"),
                       labels = c("black" = "black",
                                  "white" = "white"),
                       guide = "none") +
    coord_flip() +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 8, angle = 90),
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.5, "cm"))
}


build_compound_explorer_plot <- function(moa_long_df, rank_df, compound, models) {
  # Create a variable to plot different results
  moa_long_subset_df <- moa_long_df %>%
    dplyr::mutate(
      compound_type = ifelse(
        moa_long_df$pert_iname == compound, compound, "other"
        )
      ) %>%
    tidyr::drop_na(pert_id)
  
  moa_long_subset_df$compound_type[moa_long_subset_df$pert_iname == "bortezomib"] <- "bortezomib"
  moa_long_subset_df$compound_type[moa_long_subset_df$Metadata_broad_sample == "DMSO"] <- "DMSO"
  
  moa_long_subset_df <- moa_long_subset_df %>% dplyr::filter(model %in% models) %>%
      dplyr::as_tibble()
  
  moa_long_subset_df$original_name <- factor(moa_long_subset_df$original_name,
                                             levels = rank_df$original_name)
  
  compound_levels <- c("DMSO", "bortezomib", "other")
  if (compound != "bortezomib") {
    compound_levels <- c(compound_levels, compound)
  }

  moa_long_subset_df$compound_type <- factor(moa_long_subset_df$compound_type,
                                             levels = compound_levels)
  
  full_scatter_gg <- ggplot(moa_long_subset_df,
         aes(x = compound_type, y = model_score, fill = compound_type)) +
    geom_boxplot(outlier.size = 0.1) +
    facet_wrap("~original_name", nrow = length(models)) +
    theme_bw() +
    coord_flip() +
    xlab("") +
    ylab("Model Score") +
    theme(legend.position = "none",
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))
  
  compound_details <- moa_long_subset_df %>%
      dplyr::filter(compound_type == !!compound)
  compound_dose_gg <- ggplot(compound_details,
                             aes(x = Metadata_dose_recode, y = model_score)) +
    geom_bar(stat="identity") +
    facet_wrap("~original_name", nrow = length(models), scales = "free") +
    xlab(paste(compound, "Dose")) +
    ylab("Model Score") +
    theme_bw() +
    theme(legend.position = "none",
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))
  
  main_plot <- (
    cowplot::plot_grid(
      full_scatter_gg,
      compound_dose_gg,
      labels = c("", ""),
      ncol = 2,
      nrow = 1,
      rel_widths = c(1, 0.8),
      align = "b"
    )
  )
  main_plot
}
