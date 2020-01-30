library(dplyr)
library(shiny)
library(dqshiny)

# Load profiles
moa_file <- file.path("data", "moa_cell_health_modz.tsv.gz")
rank_file <- file.path("data", "A549_ranked_models_regression_modz.tsv")

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

rank_df <- readr::read_tsv(rank_file, col_types = readr::cols()) %>%
  dplyr::filter(shuffle_false > 0)
rank_df$target <- factor(rank_df$target, levels = rev(unique(rank_df$target)))
rank_df$original_name <- factor(rank_df$original_name,
                                levels = rev(unique(rank_df$original_name)))

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    # Application title
    titlePanel("Drug Repurposing Hub Cell Health Predictions"),
    
    # Setup multiple tabs
    tabsetPanel(
      tabPanel(
        "Model Explorer",
        # Sidebar with interactive layout
        sidebarLayout(
          
          sidebarPanel(
            helpText("Select compounds and cell health models"),
            selectInput("scatter_type",
                        label = "Select Scatter Plot Type",
                        choices = c("Cell Health", "UMAP"),
                        selected = "Cell Health"),
            selectInput("cell_health_model",
                        label = "Select Cell Health Variable",
                        choices = rank_df$original_name,
                        selected = "G1 - Number of Objects"),
            autocomplete_input("compound",
                               label = "Select a Compound",
                               options = sort(unique(moa_df$pert_iname)),
                               max_options = 10,
                               value = "bortezomib"),
            fluidRow(
              column(width = 12,
                     h4("Click and Drag Points"),
                     verbatimTextOutput("brush_info"))
            ),
            plotOutput("rank_plot",
                       height = 400,
                       width = 550)
          ),
          
          # Show a plot of the generated distribution
          mainPanel(
            plotOutput("scatter_plot",
                       height = 400,
                       width = 500,
                       brush = brushOpts(id = "plot_brush")),
            plotOutput("bar_chart",
                       height = 400,
                       width = 550)
          )
        )
      ),
      tabPanel(
        "Compound Explorer",
        # Sidebar with interactive layout
        sidebarLayout(
          
          sidebarPanel(
            helpText("Select compounds to explore"),
            autocomplete_input("compound_explorer",
                               label = "Select a Compound",
                               options = sort(unique(moa_df$pert_iname)),
                               max_options = 10,
                               value = "bortezomib"),
            checkboxGroupInput("model_select_explorer",
                               label = "Check Models to Visualize",
                               choices = rank_df$original_name,
                               selected = c("Live Cell Area",
                                            "G1 - Number of Objects",
                                            "edu positive - Number of Objects",
                                            "ROS Mean"))
          ),
          
          # Show a plot of the generated distribution
          mainPanel(
            plotOutput("compound_plot_explorer",
                       height = 400,
                       width = 500),
            plotOutput("rank_plot_explorer",
                       height = 400,
                       width = 550)
            )
          )
        )
      )
    )
)