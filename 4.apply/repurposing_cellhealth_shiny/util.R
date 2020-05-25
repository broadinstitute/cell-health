# Functions to enable shiny plot functionality

library(dplyr)
library(ggplot2)
library(ggrepel)

dye_colors <- c(
    "hoechst" = "#639B94",
    "edu" = "#E45242",
    "hoechst_edu" = "#73414b",
    "hoechst_edu_ph3" = "#7B9C32",
    "hoechst_gh2ax" = "#535f52",
    "hoechst_edu_gh2ax" = "#e37a48",
    "hoechst_edu_ph3_gh2ax" = "#E2C552",
    "draq" = "#FF6699",
    "draq_caspase" = "#7f4a72",
    'cell_rox' = "#E9DFC3",
    "dpc" = "#CA662D",
    "qc" = "black"
)

dye_labels <- c(
  "hoechst" = "Hoechst",
  "edu" = "EdU",
  "hoechst_edu" = "Hoechst + EdU",
  "hoechst_edu_ph3" = "Hoechst + EdU + PH3",
  "hoechst_gh2ax" = "Hoechst + gH2AX",
  "hoechst_edu_gh2ax" = "Hoechst + EdU + gH2AX",
  "hoechst_edu_ph3_gh2ax" = "Hoechst + EdU + PH3 + gH2AX",
  "draq" = "DRAQ7",
  "draq_caspase" = "DRAQ7 + Caspase",
  'cell_rox' = "CellROX",
  "dpc" = "DPC (Shape)",
  "qc" = "CRISPR Efficiency"
)

get_target <- function(rank, model) {
  target <- rank %>%
    dplyr::filter(target == !!model) %>%
    dplyr::pull(original_name)

  return(target)
}

build_cell_health_scatter <- function(
  moa_full_df, compound_df, control_df, model_y, model_x, target_y, target_x
  ) {

  moa_scatter_gg <- ggplot(moa_full_df,
       aes_string(x = model_x,
                  y = model_y)) +
  xlab(target_x) +
  ylab(target_y) +
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
  geom_point(data = control_df,
             aes(shape = pert_iname),
             fill = "grey",
             color = "black",
             size = 2,
             alpha = 0.4,
             show.legend = TRUE) +
  scale_color_viridis_c(name = "Dose\nLevel") +
  scale_shape_manual(
    name = "Controls",
    values = c("DMSO" = 21, "bortezomib" = 23, "MG-132" = 25),
    labels = c("DMSO" = "DMSO", "bortezomib" = "Bortezomib", "MG-132" = "MG-132")
  ) +
  theme_bw() +
  guides(
    shape = guide_legend(order = 2,
                         keywidth = 0.1,
                         keyheight = 0.1)
    )

  return(moa_scatter_gg)
}

build_umap_scatter <- function(moa_full_df, compound_df, control_df, model, target){

  umap_scatter_gg <- ggplot(moa_full_df, aes(x = umap_x, y = umap_y)) +
    xlab("UMAP X") +
    ylab("UMAP Y") +
    geom_point(aes_string(color = model),
               size = 1.25,
               pch = 16,
               alpha = 0.6) +
    geom_point(data = compound_df,
               aes_string(color = model),
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
    geom_point(data = control_df,
               aes(shape = pert_iname),
               fill = "grey",
               color = "black",
               size = 2,
               alpha = 0.4,
               show.legend = TRUE) +
    ggtitle(target) +
    scale_color_viridis_c(name = "") +
    scale_shape_manual(
      name = "Controls",
      values = c("DMSO" = 21, "bortezomib" = 23, "MG-132" = 25),
      labels = c("DMSO" = "DMSO", "bortezomib" = "Bortezomib", "MG-132" = "MG-132")
    ) +
    theme_bw() +
    guides(
      shape = guide_legend(order = 2,
                           keywidth = 0.1,
                           keyheight = 0.1)
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
    scale_fill_manual(
        name = "Assay",
        values = dye_colors,
        labels = dye_labels
    ) +
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
    tidyr::drop_na(broad_id)

  moa_long_subset_df$compound_type[is.na(moa_long_subset_df$compound_type)] <- "other"
  moa_long_subset_df$compound_type[moa_long_subset_df$pert_iname == "bortezomib"] <- "bortezomib"
  moa_long_subset_df$compound_type[moa_long_subset_df$pert_iname == "MG-132"] <- "MG-132"
  moa_long_subset_df$compound_type[moa_long_subset_df$Metadata_broad_sample == "DMSO"] <- "DMSO"

  moa_long_subset_df <- moa_long_subset_df %>% dplyr::filter(model %in% models) %>%
      dplyr::as_tibble()

  moa_long_subset_df$original_name <- factor(moa_long_subset_df$original_name,
                                             levels = rank_df$original_name)

  compound_levels <- c("DMSO", "bortezomib", "MG-132", "other")
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


load_data <- function(pos_controls=c("bortezomib", "MG-132")) {
  # Load profiles
  moa_file <- file.path("data", "moa_cell_health_modz.tsv.gz")

  # Set column dtypes for loading with readr
  moa_cols <- readr::cols(
    .default = readr::col_double(),
    Metadata_Plate_Map_Name = readr::col_character(),
    Metadata_broad_core_id = readr::col_character(),
    Metadata_broad_sample = readr::col_character(),
    Metadata_pert_well = readr::col_character(),
    broad_id = readr::col_character(),
    pert_iname = readr::col_character(),
    InChIKey14 = readr::col_character(),
    moa = readr::col_character(),
    target = readr::col_character(),
    broad_date = readr::col_character(),
    clinical_phase = readr::col_character(),
    alternative_moa = readr::col_character(),
    alternative_target = readr::col_character()
  )

  moa_df <- readr::read_tsv(moa_file, col_types = moa_cols)
  colnames(moa_df) <- gsub("cell_health_modz_target_", "", colnames(moa_df))

  # Subset DMSO from the MOA data - we always want to highlight where DMSO samples fall
  dmso_df <- moa_df %>% dplyr::filter(Metadata_broad_sample == 'DMSO')
  dmso_df$pert_iname <- "DMSO"
  pos_controls_df <- moa_df %>% dplyr::filter(pert_iname %in% pos_controls)

  rank_file <- file.path("data", "A549_ranked_models_regression_modz.tsv")
  rank_df <- readr::read_tsv(rank_file, col_types = readr::cols()) %>%
    dplyr::filter(shuffle_false > 0)
  rank_df$target <- factor(rank_df$target, levels = rev(unique(rank_df$target)))
  rank_df$original_name <- factor(rank_df$original_name,
                                  levels = rev(unique(rank_df$original_name)))

  return(
    list(
      "rank" = rank_df,
      "moa" = moa_df,
      "pos_control" = pos_controls_df,
      "dmso" = dmso_df
    )
  )
}
