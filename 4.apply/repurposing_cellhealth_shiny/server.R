suppressMessages(library(shiny))
suppressMessages(library(DT))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(cowplot))
suppressMessages(source("util.R"))

set.seed(1234)

# Load profiles
moa_file <- file.path("data", "moa_cell_health_modz.tsv.gz")
rank_file <- file.path("data", "A549_ranked_models_regression_modz.tsv")

# Set column dtypes for loading with readr
moa_cols <- readr::cols(
  .default = readr::col_double(),
  Metadata_broad_sample = readr::col_character(),
  Image_Metadata_Well = readr::col_character(),
  pert_iname = readr::col_character(),
  Metadata_broad_core_id = readr::col_character(),
  pert_id = readr::col_character(),
  pert_type = readr::col_character(),
  moa = readr::col_character()
)

moa_df <- readr::read_tsv(moa_file, col_types = moa_cols)
colnames(moa_df) <- gsub("cell_health_modz_target_", "", colnames(moa_df))

# Load performance results - this will help us to gauge trust in output predictions
rank_df <- readr::read_tsv(rank_file, col_types = readr::cols()) %>%
  dplyr::filter(shuffle_false > 0)
rank_df$target <- factor(rank_df$target, levels = rev(unique(rank_df$target)))
rank_df$original_name <- factor(rank_df$original_name,
                                levels = rev(unique(rank_df$original_name)))
model_dict_df <- rank_df %>% dplyr::select(target, original_name)
model_dict_df$target <- paste(model_dict_df$target)
model_dict_df$original_name <- paste(model_dict_df$original_name)

# Subset DMSO from the MOA data - we always want to highlight where DMSO samples fall
dmso_df <- moa_df %>%
    dplyr::filter(Metadata_broad_sample == "DMSO")

# Reshape the moa dataframe for different variable plotting
melt_id_vars <- c(
  "Image_Metadata_Well",
  "Metadata_broad_core_id",
  "Metadata_broad_sample",
  "Metadata_dose_recode",
  "umap_x",
  "umap_y",
  "pert_id",
  "pert_iname",
  "pert_type",
  "moa"
)
moa_long_df <- reshape2::melt(
  moa_df,
  id.vars = melt_id_vars,
  value.name = "model_score",
  variable.name = "model"
)
moa_long_df$model <- paste(moa_long_df$model)
moa_long_df <- moa_long_df %>%
  dplyr::left_join(model_dict_df, by = c("model" = "target"))

shinyServer(function(input, output) {
  
  # Set reactive elements
  # Compound to focus on
  compound <- reactive({
    paste(input$compound)
  })
  
  # Cell Health Model to focus on in y axis
  cell_health_model_y <- reactive({
    target_name <- input$cell_health_model_yaxis
    target <- rank_df %>%
        dplyr::filter(original_name == !!target_name) %>%
        dplyr::pull(target)
    paste(target)
  })
  
  # Cell Health Model to focus on in x axis
  cell_health_model_x <- reactive({
    target_name <- input$cell_health_model_xaxis
    target <- rank_df %>%
      dplyr::filter(original_name == !!target_name) %>%
      dplyr::pull(target)
    paste(target)
  })

  # Which scatter plot to visualize
  scatter_type <- reactive({
    paste(input$scatter_type)
  })

  # Also set reactive elements for secondary tab
  compound_explorer <- reactive({
    paste(input$compound_explorer)
  })
  
  model_select_explorer <- reactive({
    models <- input$model_select_explorer
    target <- rank_df %>%
      dplyr::filter(original_name %in% !!models) %>%
      dplyr::pull(target)
    
    paste(target)
  })
  
  # Generate a series of figures
  output$bar_chart <- renderPlot(
    {
      # Load reactive elements
      compound_select <- compound() 
      cell_health_model_select_y <- cell_health_model_y()
      cell_health_model_select_x <- cell_health_model_x()
      
      # Isolate the original name (more interpretable) of the target selected
      target_y <- get_target(rank_df, cell_health_model_select_y)
      target_x <- get_target(rank_df, cell_health_model_select_x)
      
      # Isolate results for the specific compound selected
      compound_df <- moa_df %>%
          dplyr::filter(pert_iname == !!compound_select)
    
      # Determine range of y axis
      ymax_compound_y <- max(compound_df[, cell_health_model_select_y])
      ymax_dmso_y <- max(dmso_df[, cell_health_model_select_y])
      ymax_y <- max(c(ymax_compound_y, ymax_dmso_y))
      
      ymin_compound_y <- min(compound_df[, cell_health_model_select_y])
      ymin_dmso_y <- min(dmso_df[, cell_health_model_select_y])
      ymin_y <- min(c(ymin_compound_y, ymin_dmso_y))

      ymax_compound_x <- max(compound_df[, cell_health_model_select_x])
      ymax_dmso_x <- max(dmso_df[, cell_health_model_select_x])
      ymax_x <- max(c(ymax_compound_x, ymax_dmso_x))
      
      ymin_compound_x <- min(compound_df[, cell_health_model_select_x])
      ymin_dmso_x <- min(dmso_df[, cell_health_model_select_x])
      ymin_x <- min(c(ymin_compound_x, ymin_dmso_x))
      
      # Plot! 1st - Generate the dose barplot
      bar_y_gg <- ggplot(compound_df,
                       aes_string(x = "as.factor(Metadata_dose_recode)",
                                  y = cell_health_model_select_y)) +
        geom_bar(stat = "identity", aes(fill = vb_num_live_cells)) +
        theme_bw() +
        theme(plot.title = element_text(size = 14)) +
        scale_fill_continuous(name = "# Live Cells") +
        ggtitle(compound_select) +
        xlab("Dose Recoded (low to high)") +
        ylab(target_y) +
        ylim(ymin_y, ymax_y)
      
      bar_x_gg <- ggplot(compound_df,
                         aes_string(x = "as.factor(Metadata_dose_recode)",
                                    y = cell_health_model_select_x)) +
        geom_bar(stat = "identity", aes(fill = vb_num_live_cells)) +
        theme_bw() +
        theme(plot.title = element_text(size = 14)) +
        scale_fill_continuous(name = "# Live Cells") +
        ggtitle(compound_select) +
        xlab("Dose Recoded (low to high)") +
        ylab(target_x) +
        ylim(ymin_x, ymax_x)
      
      bar_gg <- cowplot::plot_grid(
        bar_y_gg,
        bar_x_gg,
        labels = c("", ""),
        ncol = 1,
        nrow = 2
      )

      # 2nd, plot the distribution of DMSO samples along the same model
      dmso_y_gg <- ggplot(dmso_df,
                        aes_string(y = cell_health_model_select_y,
                                   x = "Metadata_broad_sample")) +
        geom_jitter(stat = "identity",
                    width = 0.2,
                    size = 3,
                    pch = 21,
                    alpha = 0.7,
                    color = "black",
                    fill = "grey") +
        theme_bw() +
        ggtitle("") +
        xlab("") +
        ylab(target_y) +
        ylim(ymin_y, ymax_y)
      
      dmso_x_gg <- ggplot(dmso_df,
                          aes_string(y = cell_health_model_select_x,
                                     x = "Metadata_broad_sample")) +
        geom_jitter(stat = "identity",
                    width = 0.2,
                    size = 3,
                    pch = 21,
                    alpha = 0.7,
                    color = "black",
                    fill = "grey") +
        theme_bw() +
        ggtitle("") +
        xlab("") +
        ylab(target_x) +
        ylim(ymin_x, ymax_x)
      
      dmso_gg <- cowplot::plot_grid(
        dmso_y_gg,
        dmso_x_gg,
        labels = c("", ""),
        ncol = 1,
        nrow = 2
      )
      
      # Combine these two plots into a single figure called "barchart"
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
  
  # The second group of figures is the scatter plot
  output$scatter_plot <- renderPlot(
    {
      # Load reactive variables
      compound_select <- compound() 
      cell_health_model_select_y <- cell_health_model_y()
      cell_health_model_select_x <- cell_health_model_x()
      scatter_plot_type <- scatter_type()

      # Pull cell health model of interest
      target_y <- get_target(rank_df, cell_health_model_select_y)
      target_x <- get_target(rank_df, cell_health_model_select_x)
      
      # Isolate specific compound
      moa_compound_df <- moa_df %>%
          dplyr::filter(pert_iname == !!compound_select)
      
      # Now, visualize cell health models
      if (scatter_plot_type == "Cell Health") {
        scatter_gg <- build_cell_health_scatter(
          moa_df, moa_compound_df, dmso_df,
          cell_health_model_select_y, cell_health_model_select_x, target_y, target_x
        )
      } else {
        scatter_gg <- build_umap_scatter(
          moa_df, moa_compound_df, dmso_df, cell_health_model_select_y, target_y
        )
      }
      scatter_gg
    }
  )
  
  # Build the model ranking figure
  output$rank_plot <- renderPlot({
    # Load reactive elements
    cell_health_model_select_y <- cell_health_model_y()
    cell_health_model_select_x <- cell_health_model_x()
    
    # Highlight the selected model
    rank_df$to_highlight <- ifelse(
      rank_df$target %in% c(cell_health_model_select_y, cell_health_model_select_x),
      "black", "white"
      )
    
    # Create and output the ranking figure
    build_rank_plot(rank_df)
  })
  
  output$rank_plot_explorer <- renderPlot({
    # Highlight the selected model
    selected_models <- model_select_explorer()
    rank_df$to_highlight <- ifelse(
      rank_df$target %in% selected_models, "black", "white"
    )
    
    # Create and output the ranking figure
    build_rank_plot(rank_df)
  })
  
  output$compound_plot_explorer <- renderPlot({
    # Load reactive elements
    compound_select_explore <- compound_explorer() 
    selected_models <- model_select_explorer()
    
    build_compound_explorer_plot(
      moa_long_df,
      rank_df,
      compound_select_explore,
      selected_models)
  })
  
  # Build output text for print rendering
  output$brush_info <- renderTable({
    
    # Load reactive elements
    cell_health_model_select_y <- cell_health_model_y()
    cell_health_model_select_x <- cell_health_model_x()
    
    scatter_plot_type <- scatter_type()
  
    # Select which columns to display
    if (scatter_plot_type == "Cell Health") {
      output_df <- moa_df %>%
        dplyr::select(
          pert_iname,
          moa,
          Metadata_broad_sample,
          Metadata_dose_recode,
          !!cell_health_model_select_y,
          !!cell_health_model_select_x) %>%
          dplyr::arrange(desc(!!sym(cell_health_model_select_y)))
    } else {
      output_df <- moa_df %>%
        dplyr::select(pert_iname,
                      moa,
                      Metadata_broad_sample,
                      umap_x,
                      umap_y,
                      Metadata_dose_recode,
                      !!cell_health_model_select_y,
                      !!cell_health_model_select_x) %>%
        dplyr::arrange(desc(!!sym(cell_health_model_select_y)))
    }
    
    output_df <- output_df %>%
        dplyr::rename(perturbation = pert_iname,
                      sample = Metadata_broad_sample,
                      dose = Metadata_dose_recode) %>%
      mutate_if(is.numeric, round, 3)
    
    brushedPoints(output_df, input$plot_brush)
    })
  })
