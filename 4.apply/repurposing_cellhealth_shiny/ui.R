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
      tabPanel(
        "Model Explorer",
        # Sidebar with interactive layout
        sidebarLayout(
          
          sidebarPanel(
            p(strong("Getting Started:", style="color:red"),
              a("Documentation" ,href="https://github.com/broadinstitute/cell-health/tree/master/4.apply/repurposing_cellhealth_shiny/")),
            helpText("Select compounds and cell health models. The points represent drug perturbation Consensus Profiles."),
            selectInput("scatter_type",
                        label = "Select Scatter Plot Type",
                        choices = c("Cell Health", "UMAP"),
                        selected = "Cell Health"),
            selectInput("cell_health_model_yaxis",
                        label = "Select Cell Health Model for Y axis (UMAP Toggle)",
                        choices = sort(paste(rank_df$readable_name)),
                        selected = "G1 - # cells"),
            selectInput("cell_health_model_xaxis",
                        label = "Select Cell Health Model for X axis",
                        choices = sort(paste(rank_df$readable_name)),
                        selected = "ROS"),
            autocomplete_input("compound",
                               label = "Select a Compound to Highlight",
                               options = sort(unique(moa_df$pert_iname)),
                               max_options = 10,
                               value = "HMN-214"),
            checkboxInput("remove_controls",
                          label = "Remove Controls from Click and Drag",
                          value = FALSE),
            downloadButton("downloadData", "Download Point Selection", class = "button"),
            tags$head(tags$style(".button{background-color:orange;} .button{color: black;}")) 
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
            p(strong("Getting Started:", style = "color:red"),
              a("Documentation", href="https://github.com/broadinstitute/cell-health/tree/master/4.apply/repurposing_cellhealth_shiny/")),
            helpText("Select compounds to explore"),
            autocomplete_input("compound_explorer",
                               label = "Select a Compound",
                               options = sort(unique(moa_df$pert_iname)),
                               max_options = 10,
                               value = "HMN-214"),
            checkboxGroupInput("model_select_explorer",
                               label = "Check Models to Visualize",
                               choices = rank_df$readable_name,
                               selected = c("Live Cell Area",
                                            "G1 - # cells",
                                            "S - Intensity Nucleus EdU Mean",
                                            "ROS"))
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
        ),
      tabPanel(
        "Dose Fit Explorer",
        sidebarLayout(
          sidebarPanel(
            p(strong("Getting Started:", style = "color:red"),
              a("Documentation", href="https://github.com/broadinstitute/cell-health/tree/master/4.apply/repurposing_cellhealth_shiny/")),
            helpText("Select compounds to explore"),
            autocomplete_input("dose_explorer",
                               label = "Select a Compound",
                               options = sort(unique(moa_df$pert_iname)),
                               max_options = 10,
                               value = "HMN-214"),
            selectInput("dose_model_select",
                        label = "Select Cell Health Model to Fit Dose Response Curve",
                        choices = sort(paste(rank_df$readable_name)),
                        selected = "G1 - # cells")
            ),
          mainPanel(
            plotOutput("dose_fit_curve",
                       height = 400,
                       width = 550)
            )
          ),
        fluidRow(
          h3("Dose Fitting Information"),
          tableOutput("dose_info")
          )
        )
      )
  )
)