#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(DT) 

fluidPage(
  tags$head(
    tags$title("AlzPPMap")
  ),
  titlePanel(
    div(
      "AlzPPMap: Alzheimer's disease Proteomics-PPI integrated network Map",
      style = "font-size: 36px; font-weight: bold;"
    )
  ),
  
  div(
    p("This application allows you to analyze selected gene's subnetwork for the integrated network. Select a number and genes to visualize the integrated network.",
      style = "font-size: 18px;"),
    p("The background methodology is described in...",
      style = "font-size: 18px;")
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      uiOutput("gene_selector"),
      uiOutput("slider_inputs"),
      actionButton("run_analysis", "RUN")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Network", plotlyOutput("network_plot", height = "600px")),  
        tabPanel("Table", DTOutput("community_table", width="70%"),
                 downloadButton("download_table", "Download CSV")),
        tabPanel("Graph Download", downloadButton("graph_download", "Download Graph RDS"))
      )
    )
  )
)
