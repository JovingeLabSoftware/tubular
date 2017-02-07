# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(scde)
library(Biobase)
library(ggthemes)
library(GGally)
library(shiny)
library(threejs)

eset <- readRDS("./data/eset_scde.rds")
tsne <- readRDS("data/tsne_pre.rds")
fd <- fData(eset)
ui_theme <- function (base_size = 16, base_family = "") 
{
  rebase <- ggthemes_data$solarized$base[c(paste0("base0", 3:0),
                                           paste0("base", 0:3))]
  ret <- (theme_foundation(base_size = base_size, base_family = base_family) + 
            theme(text = element_text(color = "#c8c8c8"), 
                  title = element_text(color = rebase["rebase0"]), 
                  line = element_line(color = rebase["rebase01"]), 
                  rect = element_rect(fill = "#272b30", color = NA), 
                  axis.ticks = element_line(color = rebase["rebase01"]), 
                  axis.line = element_line(color = "#c8c8c8", 
                                           linetype = 1), axis.title.y = element_text(angle = 90), 
                  legend.background = element_rect(fill = NULL, color = "NA"), 
                  legend.key = element_rect(fill = NULL, colour = NULL, 
                                            linetype = 0), panel.background = element_rect(fill = rebase["rebase02"], 
                                                                                           colour = NA), panel.border = element_blank(), 
                  panel.grid = element_line(color = rebase["rebase03"]), 
                  panel.grid.major = element_line(color = rebase["rebase03"]), 
                  panel.grid.minor = element_line(color = rebase["rebase03"], 
                                                  size = 0.25), plot.background = element_rect(fill = NULL, 
                                                                                               colour = NULL, linetype = 0)))
  ret
}


shinyServer(function(input, output) {
  getScheme <- reactive({
    if(colorBy() %in% colnames(pData(eset))) {
      cb <- colorBy()
    } else if(colorBy() == "Genes") {
      cb <- "gene"
    } else {
      cb <- NULL
    } 
    
    if(shapeBy() %in% colnames(pData(eset))) {
      shb <- shapeBy()
    } else if(shapeBy() == "Genes") {
      shb <- "gene"
    } else {
      shb <- NULL
    } 
    
    if(sizeBy() %in% colnames(pData(eset))) {
      sib <- sizeBy()
    } else if(sizeBy() == "Genes") {
      sib <- "gene"
    } else {
      sib <- NULL
    } 
    list(cb=cb, shb=shb, sib=sib)
  })
  
  filterBy <- reactive({
    filters <- list();
    c <- 1
    if(nchar(input$filter1) > 1 && tolower(input$filter1) %in% tolower(fd$hsapiens_homolog_associated_gene_name)) {
      ix <- which(tolower(fd$hsapiens_homolog_associated_gene_name) == tolower(input$filter1))
      filters[[c]] <- list(id=rownames(exprs(eset))[ix], op=input$op1, value=input$val1)
      c <- c+1
    }
    if(nchar(input$filter2) > 1 && tolower(input$filter2) %in% tolower(fd$hsapiens_homolog_associated_gene_name)) {
      ix <- which(tolower(fd$hsapiens_homolog_associated_gene_name) == tolower(input$filter2))
      filters[[c]] <- list(id=rownames(exprs(eset))[ix], op=input$op2, value=input$val2)
      c <- c+1
    }
    filters
  })

  colorBy <- reactive({
    if(input$color_by == "None") {
      return("")
    } else {
      return(input$color_by)
    }
  })
  
  shapeBy <- reactive({
    if(input$shape_by == "None") {
      return("")
    } else {
      return(input$shape_by)
    }
  })

  sizeBy <- reactive({
    if(input$size_by == "None") {
      return("")
    } else {
      return(input$size_by)
    }
  })

  geneList <- reactive({
    genes <- strsplit(input$genes, "\n")
    ix <- which(tolower(genes[[1]]) %in% tolower(fData(eset)$hsapiens_homolog_associated_gene_name))
    if(length(ix)) {
      return(genes[[1]][ix])
      ix2 <- which(tolower(fData(eset)$hsapiens_homolog_associated_gene_name) %in% tolower(genes[[1]]))
    } else {
      return(integer(0))
    }
  })

  getEset <- reactive({
    temp <- eset
    if(length(filterBy())) {
      for(f in filterBy()) {
        if(f$op == "<") {
          temp <- temp[, which(exprs(temp)[f$id,] < f$value)]
        } else {
          temp <- temp[, which(exprs(temp)[f$id,] > f$value)]
        }
      }
    }
    ix <- which(tolower(fData(temp)$external_gene_name) %in% tolower(geneList()))
    if(length(ix) > 1) {
      vals <- apply(exprs(temp)[ix,], 2, sum)
    } else if (length(ix) == 1) {
      vals <- exprs(temp)[ix,]  
    } else {
      vals <- rep(NA, ncol(temp))
    }
    temp$gene = vals
    temp
  })
  
  output$distPlot2 <- renderPlot({
    temp <- getEset()
    if(length(filterBy())) {
      for(f in filterBy()) {
        if(f$op == "<") {
          temp <- temp[, which(exprs(temp)[f$id,] < f$value)]
        } else {
          temp <- temp[, which(exprs(temp)[f$id,] > f$value)]
        }
      }
    }
    sch <- getScheme()
    plotReducedDim(temp, colour_by=sch$cb, shape_by=sch$shb, size_by=sch$sib) + ui_theme()
  })
  
  output$distPlot <- renderPlot({
    temp <- getEset()
    plotReducedDim(temp, colour_by="gene") + ui_theme()
    #plotReducedDim(eset, colour_by="gene") + ui_theme()   
  })

  output$densityPlot <- renderPlot({
    temp <- getEset()
    df = data.frame(values=temp$gene, groups=pData(temp)$Zone)
    ggplot(data = df, aes(x = groups, y = values)) +
      geom_violin(fill=NA, col="gray") + 
      geom_jitter(width = 0.05, size=0.1, color="orange") +
      ggtitle("Gene distribution by zone") +
      ui_theme()
  })
  
  output$xy <- renderPlot({
    genes <- geneList()
    validate(
      need((length(genes) > 1 & length(genes) < 11), paste("\n\n", length(genes), "valid gene(s) identified.  Please enter between 2 and 10 genes."))
    )
    ix <- which(tolower(fData(eset)$hsapiens_homolog_associated_gene_name) %in% tolower(genes))
    df <- as.data.frame(exprs(eset)[ix,])
    df <- as.data.frame(t(df))
    df$Zone_simple <- eset$Zone_simple

    dens <- function(data, mapping, ...) {
      ggplot(data = data, mapping=mapping) +
        geom_density(..., alpha = 0.9) 
    }
    ggpairs(as.data.frame(df), 
            ggplot2::aes(colour=Zone_simple), 
            columnLabels=fData(eset)$hsapiens_homolog_associated_gene_name[ix],
            columns=1:length(ix), 
            diag = list(continuous = dens),
            upper = list(continuous = wrap("cor", size = 5)),
            legend=c(length(ix), 1)) +
    ui_theme()
  }, height=800)

  output$three <- renderScatterplotThree({
    sch <- getScheme()
    temp <- getEset()
    p <- plotReducedDim(temp, colour_by=sch$cb, shape_by=sch$shb, size_by=sch$sib) + ui_theme()
    p <- ggplot_build(p)$data[[2]]
    cols <- p$fill
    size <- p$size
    
    scatterplot3js(tsne$Y, bg = "#272b30", 
                   color=cols, size=size/4, 
                   height = 800, 
                   axes=FALSE, 
                   flip.y = FALSE,
                   pch="@",   
                   grid = FALSE)
  })
})


## this is how tsne matrix can pre-generated
## taken directly from scater package
#  
#
#  eset <- readRDS("data/eset_scde.rds")
#  rv <- matrixStats::rowVars(exprs(eset))
#  feature_set <- order(rv, decreasing = TRUE)[seq_len(500)]
#  exprs_to_plot <- exprs(eset)[feature_set,]
#  keep_feature <- (matrixStats::rowVars(exprs_to_plot) > 0.001)
#  keep_feature[is.na(keep_feature)] <- FALSE
#  exprs_to_plot <- exprs_to_plot[keep_feature, ]
  
## Standardise expression if stand_exprs(object) is null
#  exprs_to_plot <- t(scale(t(exprs_to_plot), scale = TRUE))
  
#  set.seed(100)
#  tsne_pre <- Rtsne(t(exprs_to_plot), dims=3,
#                    initial_dims = 50, 
#                    perplexity = 8)
# saveRDS(tsne_pre, file="data/tsne_pre.rds")  
  