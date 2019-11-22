library(dplyr)
library(shiny)
library(dqshiny)

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

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    # Application title
    titlePanel("Drug Repurposing Hub Cell Health Predictions"),
    
    # Sidebar with interactive layout
    sidebarLayout(
      
      sidebarPanel(
        helpText("Select compounds and cell health models"),
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
    )
  )