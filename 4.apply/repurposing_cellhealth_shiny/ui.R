library(dplyr)
library(shiny)
library(dqshiny)
suppressMessages(source("util.R"))

# Load profiles
data <- load_data()
moa_df <- data[["moa"]]
rank_df <- data[["rank"]]

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(theme = "bootstrap.min.css",
    # Application title
    titlePanel("Drug Repurposing Hub Cell Health Predictions"),
    
    # Setup multiple tabs
    tabsetPanel(
      tabPanel("Getting Started"),
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
            selectInput("cell_health_model_yaxis",
                        label = "Select Cell Health Variable for Y axis (UMAP Toggle)",
                        choices = rank_df$original_name,
                        selected = "G1 - Number of Objects"),
            selectInput("cell_health_model_xaxis",
                        label = "Select Cell Health Variable for X axis",
                        choices = rank_df$original_name,
                        selected = "# Live Cells"),
            autocomplete_input("compound",
                               label = "Select a Compound",
                               options = sort(unique(moa_df$pert_iname)),
                               max_options = 10,
                               value = "YM-155"),
            checkboxInput("remove_controls",
                          label = "Remove Controls from Click and Drag",
                          value = FALSE)
            ),

          plotOutput("scatter_plot",
                     height = 400,
                     width = 500,
                     brush = brushOpts("plot_brush"))
          ),
          
          # Show a plot of the generated distribution
          mainPanel(
            fluidRow(
              plotOutput("rank_plot",
                         height = 400,
                         width = 550),
              plotOutput("bar_chart",
                         height = 400,
                         width = 550)
            ),
            fluidRow(
              h3("Click and Drag Points"),
              tableOutput("brush_info")
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
                               value = "YM-155"),
            checkboxGroupInput("model_select_explorer",
                               label = "Check Models to Visualize",
                               choices = rank_df$original_name,
                               selected = c("Live Cell Area [um2]",
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