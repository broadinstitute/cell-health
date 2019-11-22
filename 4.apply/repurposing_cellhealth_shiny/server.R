suppressMessages(library(shiny))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(cowplot))

set.seed(1234)

# Load profiles
moa_file <- file.path("data", "moa_cell_health.tsv.gz")
rank_file <- file.path("data", "A549_ranked_models.tsv")

moa_cols <- readr::cols(
  .default = readr::col_double(),
  Metadata_broad_sample = readr::col_character(),
  Image_Metadata_Well = readr::col_character(),
  pert_iname = readr::col_character(),
  vendor = readr::col_character(),
  catalog_no = readr::col_character(),
  vendor_name = readr::col_character(),
  smiles = readr::col_character(),
  InChIKey = readr::col_character(),
  deprecated_broad_id = readr::col_character()
)

moa_df <- readr::read_tsv(moa_file, col_types = moa_cols)

rank_df <- readr::read_tsv(rank_file, col_types = readr::cols()) %>%
  dplyr::filter(shuffle_false > 0)
rank_df$target <- factor(rank_df$target, levels = rev(unique(rank_df$target)))
rank_df$original_name <- factor(rank_df$original_name,
                                levels = rev(unique(rank_df$original_name)))
# Subset moa_df
moa_cols <- c("Metadata_broad_sample", "Metadata_dose_recode", "Image_Metadata_Well",
               "umap_x", "umap_y", "pert_iname", paste(rank_df$target))

dmso_df <- moa_df %>%
    dplyr::filter(Metadata_broad_sample == "DMSO")

shinyServer(function(input, output) {
  
  compound <- reactive({
    paste(input$compound)
  })
  
  cell_health_model <- reactive({
    target_name <- input$cell_health_model
    target <- rank_df %>%
        dplyr::filter(original_name == !!target_name) %>%
        dplyr::pull(target)
    paste(target)
  })

  output$bar_chart <- renderPlot(
    {
      compound_select <- compound() 
      cell_health_model_select <- cell_health_model()
      target <- rank_df %>%
        dplyr::filter(target == !!cell_health_model_select) %>%
          dplyr::pull(original_name)
      
      compound_df <- moa_df %>%
          dplyr::filter(pert_iname == !!compound_select)
    
      bar_gg <- ggplot(compound_df,
                       aes_string(x = "as.factor(Metadata_dose_recode)",
                                  y = cell_health_model_select)) +
        geom_bar(stat = "identity", aes(fill = vb_num_live_cells)) +
        theme_bw() +
        theme(plot.title = element_text(size = 14)) +
        scale_fill_continuous(name = "# Live Cells") +
        ggtitle(compound_select) +
        xlab("Dose Recoded (low to high)") +
        ylab(target)
      
      dmso_gg <- ggplot(dmso_df,
                        aes_string(y = cell_health_model_select,
                                   x = "Metadata_broad_sample")) +
        geom_jitter(stat = "identity", width = 0.2) +
        theme_bw() +
        ggtitle("") +
        xlab("") +
        ylab(target)
      
      main_plot <- (
        cowplot::plot_grid(
          dmso_gg,
          bar_gg,
          labels = c("", ""),
          ncol = 2,
          nrow = 1,
          rel_widths = c(0.4, 1)
        )
      )
      main_plot
    }
    )
  
  output$scatter_plot <- renderPlot(
    {
      compound_select <- compound() 
      cell_health_model_select <- cell_health_model()
      target <- rank_df %>%
        dplyr::filter(target == !!cell_health_model_select) %>%
        dplyr::pull(original_name)
      
      moa_compound_df <- moa_df %>%
          dplyr::filter(pert_iname == !!compound_select)
      
      scatter_gg <- ggplot(moa_df,
                           aes_string(x = "vb_num_live_cells",
                                      y = cell_health_model_select)) +
        xlab("# Live Cells") +
        ylab(target) +
        geom_point(aes(color = Metadata_dose_recode),
                   size = 1.25,
                   pch = 16,
                   alpha = 0.6) +
        geom_point(data = moa_compound_df,
                   aes(fill = Metadata_dose_recode),
                   color = "black",
                   size = 5,
                   pch = 23,
                   alpha = 0.7) +
        geom_point(data = dmso_df,
                   fill = "grey",
                   size = 3,
                   pch = 21,
                   alpha = 0.7) +
        scale_color_viridis_c(name = "Dose\nLevel") +
        scale_fill_viridis_c(name = "Dose\nLevel", guide = "none") +
        theme_bw()
      scatter_gg
    }
  )
  
  output$brush_info <- renderPrint({
    cell_health_model_select <- cell_health_model()
    
    output_df <- moa_df %>%
      dplyr::select(pert_iname,
                    Metadata_broad_sample,
                    Metadata_dose_recode,
                    !!cell_health_model_select,
                    vb_num_live_cells) %>%
        dplyr::arrange(desc(!!sym(cell_health_model_select)))

    brushedPoints(output_df, input$plot_brush)
    })
  
  output$rank_plot <- renderPlot({
    cell_health_model_select <- cell_health_model()
    
    rank_df$to_highlight <- ifelse(rank_df$target == cell_health_model_select, "black", "white")

    rank_gg <- ggplot(rank_df ,
           aes(x = original_name,
               y = shuffle_false)) +
      geom_bar(aes(fill = assay, color = to_highlight), lwd = 1.25, stat="identity") +
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
    
    rank_gg
  })
  })
