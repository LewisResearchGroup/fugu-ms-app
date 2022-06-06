# Fugo App

library(shiny)
library(shinydashboard)
source("fugoMS_app.R")


# Helper Functions for data analysis ----
stats_choices <- c('N/A' = 1, 'p-value' = 2)
# stats_choices <- c('N/A' = 1, 'p-value' = 2, 'Threshold & Fold Change' = 3, 'p-value, Threshold & Fold Change' = 4)

PValueCalculation <- function(dat, pThresh = 1){
  
  alpha_val <- 0.05
  pDat <- fugoStats(dat, pCalc = T, pThresh = pThresh)
  num_samples = nrow(pDat)
  idx_p <- pDat$pVal < alpha_val/(num_samples)
  pDat <- pDat[idx_p,]
  
  # Remove NA data from frame
  pDat <- na.omit(pDat)
  
  return(pDat)
}

ThresholdFoldCalculation <- function(dat, threshold = 20000, fold = 4){
  # Calculate the average threshold per species 
  mDat <- fugoStats(dat, avgRep = T)      
  
  # Remove markers if the max average across all groups < given threshold
  thresh_num <- threshold
  idx_num <- apply(abs(metaSep(mDat)$data), 1, max) > thresh_num
  
  # filter markers from data set if the threshold is less than 
  # 4-fold change in comparison to MHB
  thresh_fold <- fold
  # select MHB or 1 as reference when calling this function 
  fDat <- fugoStats(mDat, scale = 'fold')        
  idx_fold <- (apply(abs(metaSep(fDat)$data), 1, max) > thresh_fold)
  
  ## Find data with correct intensity and fold change
  idx <- which((idx_num + idx_fold) == 2)
  out <- pDat[idx,]
  
  return(out)
}

getDownloadButton <- function(inputId){
  downloadButton(outputId = inputId, label = 'Download', icon = icon("download"), style = 'float: right;')
}

# Define UI for app that draws a histogram ----
ui <- dashboardPage(
  skin = 'yellow',
  
  # App title
  dashboardHeader(
    title = "FUGO-MS",
    titleWidth = 300
  ),
  
  dashboardSidebar(
    sidebarMenu(
      div(
        style = 'padding: 10px 0px -5px;',
        fileInput(inputId = "file1", label = "Upload Data File",
                  multiple = FALSE,
                  accept = ".csv")
      ),
      
      hr(),
      fluidRow(
        h4('Statistical Analysis (Optional)'),
        style = 'padding: 10px; text-align: center'
      ),
      
      fluidRow(
        column(
          width = 6,
          numericInput('pvalueNum', 'p-value', 1, min = 0, max = 1,)
        ),
        column(
          width = 6,
          numericInput('thresholdNum', 'Threshold', 20000)
        )
      ),
      fluidRow(
        column(
          width = 6,
          numericInput('foldNum', 'Fold Change', 4, min = 1, max = 100)
        ),
        # column(
        #   width = 6,
        #   numericInput('clusteringNum', 'Clustering', 10, min = 1, max = 100)
        # )
      ),
      
      # checkboxInput("sandbox", "p-value", value = FALSE),
      
      selectInput(inputId = "statAnalysisSelect",
                  label = "Choose Statistical Analysis:",
                  choices = stats_choices,
                  # selected = "p-value",
                  selectize = TRUE,
                  multiple = FALSE),
      
      br(),
      
      # actionButton(
      #   inputId = "clearButton", 
      #   label = "Reset", 
      #   class = "btn-danger",
      #   style = 'width: 40%; margin-left: 30%; margin-right: 30%'
      # ),
      
      actionButton(
        inputId = "applyButton", 
        label = "Apply", 
        class = "btn-success",
        style = 'width: 40%; margin-left: 30%; margin-right: 30%'
      )
      
      # style='margin:20px auto; border-radius: 25px; box-shadow: rgba(0, 0, 0, 0.35) 0px 5px 15px;'
    ),
    width = 300
  ),
  
  dashboardBody(
    tags$script(HTML("$('body').addClass('fixed');")),
    
    fluidRow(
      box(
        title = 'Data Table',
        solidHeader = TRUE,
        status = 'warning',
        width = 12,
        # height = 500,
        div(
          DT::dataTableOutput('contents', width = '100%'),
          style = 'overflow-y: scroll;overflow-x: scroll;'
        ),
        br(),
        getDownloadButton('downloadDataTable')
      )
    ),
    
    # PCA, Bar Plot
    fluidRow(
      # div(
      #   class = "col-md-6 col-lg-4",
      # ),
      box(
        title = 'Heat Map',
        solidHeader = TRUE,
        status = 'warning',
        width = 6,
        # height = 500,
        textOutput('textHeatMap'),
        plotOutput("heatmap", width = "100%"),
        getDownloadButton('downloadHeatMap')
      ),
      box(
        title = 'Violin Plot',
        solidHeader = TRUE,
        status = 'warning',
        width = 6,
        # height = 500,
        textOutput('textViolin'),
        plotOutput("violin", width = "100%"),
        getDownloadButton('downloadViolin')
      )
    ),
    fluidRow(
      box(
        title = 'Bar Plot',
        solidHeader = TRUE,
        status = 'warning',
        width = 6,
        # height = 500,
        textOutput('textBar'),
        plotOutput("barplot", width = "100%"),
        getDownloadButton('downloadBar')
      ),
      box(
        title = 'PCA Plot',
        solidHeader = TRUE,
        status = 'warning',
        width = 6,
        # height = 500,
        textOutput('textPca'),
        plotOutput("pca", width = "100%"),
        getDownloadButton('downloadPCA')
      )
    )
  )
)

# Define Server for app ----
server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2)
  
  observeEvent(input$applyButton,{
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    data <- read.csv(file$datapath, header = TRUE)
    plotData <- data
    threshold <- input$thresholdNum
    pvalue <- input$pvalueNum
    fold <- input$foldNum
    
    statSelection <- input$statAnalysisSelect
    if(statSelection == 1) {
      plotData <- data
    } else if(statSelection == 2) {
      plotData <- PValueCalculation(data)
    } else if(statSelection == 3) {
      plotData <- ThresholdFoldCalculation(data, pvalue)
    } else if(statSelection == 4) {
      pdata <- PValueCalculation(data)
      plotData <- ThresholdFoldCalculation(pdata, threshold, fold)
    }
    
    output$contents <- DT::renderDataTable({
      plotData
    })
    
    output$heatmap <- renderPlot({
      fugoPlot(plotData, heatMap = T, scale = 'row', rCst = T, grid = F, cCst = F)
    },
    res = 96
    )
    
    output$violin <- renderPlot({
      # fugoPlot(plotData,vioPlot = T)
    },
    res = 96
    )
    
    underConstructionMessage <- 'Not currently available. Refer to FugoCMDScript for command line functions to visualize data.'
    
    output$textPca <- renderText({
      paste(underConstructionMessage)
    })
    
    output$textBar <- renderText({
      paste(underConstructionMessage)
    })
    
    output$textViolin <- renderText({
      paste(underConstructionMessage)
    })
    
    output$downloadHeatMap <- downloadHandler(
      
      filename = "fugo_HeatMap.png",
      content = function(file) {
        png(file)
        print(fugoPlot(plotData, heatMap = T, scale = 'row', rCst = T, grid = F, cCst = F))
        dev.off()
      })
    
    output$downloadDataTable <- downloadHandler(
      filename = "fugo_HeatMap.csv",
      content = function(file) {
        write.csv(plotData, file)
      }) 
    
  })
}

# Deploy app ----
shinyApp(ui = ui, server = server)
