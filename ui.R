
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinythemes)
library(Biobase)
library(threejs)

eset <- readRDS("./data/eset_scde.rds")
tsne <- readRDS("data/tsne_pre.rds")

textareaInput <- function(id, label, value, rows=20, cols=35, class="form-control"){
  tags$div(
    class="form-group shiny-input-container",
    tags$label('for'=id,label),
    tags$textarea(id=id,class=class,rows=rows,cols=cols,value))
}


shinyUI(
  fluidPage(
    theme = shinytheme("slate"),
    tags$head(
      tags$style(
        ".selectize-dropdown, .selectize-input, .selectize-input { 
           line-height: 23px; 
        }"
      )
    ),
    titlePanel("Heart Tube e10.5"),
  
    sidebarLayout(
      sidebarPanel(width=3,
        textareaInput(id="genes", label="Enter some genes", value="NKX2-5", rows=3),
        selectInput("color_by", "Color by:", c("None", "Genes", colnames(pData(eset))), selected = "Zone", multiple = FALSE),
        selectInput("shape_by", "Shape by:", c("None", "Genes", colnames(pData(eset))), selected = "None", multiple = FALSE),
        selectInput("size_by", "Size by:", c("None", "Genes", colnames(pData(eset))), selected = "None", multiple = FALSE),
  
        div(
          div(style="display: inline-block;vertical-align:top;", textInput("filter1", label="Filter", width=100)),
          div(style="display: inline-block;vertical-align:top;", selectInput("op1", label="Op", choices=c("<", ">"), width=60)),
          div(style="display: inline-block;vertical-align:top;", textInput("val1", label="Threshold", width=100))
        ),
        div(
          div(style="display: inline-block;vertical-align:top;", textInput("filter2", label="Filter", width=100)),
          div(style="display: inline-block;vertical-align:top;", selectInput("op2", label="Op", choices=c("<", ">"), width=60)),
          div(style="display: inline-block;vertical-align:top;", textInput("val2", label="Threshold", width=100))
        )
      ),
  
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
          tabPanel("Plot",     
             fluidRow(
               column(12, plotOutput("densityPlot"))
             ),
             fluidRow(
               column(6, plotOutput("distPlot")),
               column(6, plotOutput("distPlot2"))
             )
          ),
          tabPanel("X vs. Y",
              plotOutput("xy")
          ),
          tabPanel("3D",
               scatterplotThreeOutput("three", height="800px")
          )
          
        )
      )
    )
#      h3(textOutput("caption", container = span)),
#      plotOutput("densityPlot"),
#      plotOutput("distPlot"),
#      plotOutput("distPlot2")
#    )
#  )
))

# verbatimTextOutput("summary")

