# Fugo App

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinybusy)

source("fugoMS_app.R")


# Helper Functions for data analysis ----
# stats_choices <- c('P-Value' = 1, 'Threshold' = 2, 'Fold Change' = 3)
stats_choices <- c('P-Value' = 1, 'Threshold' = 2)

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

ThresholdCalculation <- function(dat, threshold = 20000){
  # Calculate the average threshold per species 
  mDat <- fugoStats(dat, avgRep = T)      
  
  # Remove markers if the max average across all groups < given threshold
  thresh_num <- threshold
  idx_num <- apply(abs(metaSep(mDat)$data), 1, max) > thresh_num
  
  ## Find data with correct intensity and fold change
  idx <- which((idx_num) == 1)
  out <- mDat[idx,]
  
  return(out)
}

FoldCalculation <- function(dat, fold = 4, reference = 1){
  # Calculate the average threshold per species 
  mDat <- fugoStats(dat, avgRep = T)      
  
  # filter markers from data set if the threshold is less than 
  # 4-fold change in comparison to MHB
  thresh_fold <- fold

  # select MHB or 1 as reference when calling this function
  fDat <- fugoStats(mDat, scale = 'fold')
  idx_fold <- (apply(abs(metaSep(fDat)$data), 1, max) > thresh_fold)

  ## Find data with correct intensity and fold change
  idx <- which((idx_fold) == 1)
  out <- mDat[idx,]
  
  return(out)
}

getDownloadButton <- function(inputId){
  downloadButton(outputId = inputId, label = 'Download', icon = icon("download"), style = 'float: right;')
}

sendWarningAlert <- function(session, headerText, messageText, isSuccess = FALSE){
  alertType <- 'warning'
  if(isSuccess){
    alertType <- 'success'
  }
  
  sendSweetAlert(
    session = session,
    title = headerText,
    text = messageText,
    type = alertType
  )
}

# Define UI for app that draws a histogram ----
ui <- dashboardPage(
  skin = 'yellow',

  # App title
  dashboardHeader(
    title = "FUGU-MS",
    titleWidth = 350
  ),
  
  dashboardSidebar(
    sidebarMenu(
      div(
        style = 'padding-top: 10px;',
        div(
          style = 'padding: 20px',
          imageOutput(
            outputId = "fugoLogo",
            height = '200px',
          )
        ),
        fileInput(inputId = "file1", label = "Upload Data File",
                  width = '100%',
                  multiple = FALSE,
                  accept = ".csv")
      ),
      
      hr(),
      fluidRow(
        h3('Statistical Analysis'),
        style = 'text-align: center; width: auto;'
      ),
      
      div(
        style = 'padding: 0 15px',
        fluidRow(
          column(
            width = 6,
            style = 'padding: 0;',
            numericInput('pvalueNum', 'P-Value', 1, min = 0, max = 1, step = 0.1, width = '100%'),
          ),
          column(
            width = 6,
            style = 'padding: 0;',
            numericInput('thresholdNum', 'Threshold', 20000, step = 10000, width = '100%'),
          )
        ),
        fluidRow(
          column(
            width = 6,
            style = 'padding: 0;',
            numericInput('foldNum', 'Fold Change', 4, min = 1, max = 100, width = '100%'),
          ),
          # column(
          #   width = 6,
          #   style = 'padding: 0;',
          #   selectInput(
          #     inputId = 'referencePicker',
          #     label = 'Reference',
          #     choices = c("MHB" = 1, "Transmission" = 2, "Gears" = 3),
          #     width = '100%',
          #     # options = list(style = "btn-primary")
          #   ),
          # )
        ),
      ),

      checkboxGroupButtons(
        inputId = "statAnalysisSelect",
        label = "Choose analysis to apply:",
        choices = stats_choices,
        justified = TRUE,
        checkIcon = list(
          yes = tags$i(class = "fa fa-check-square", 
                       style = "color: orange"),
          no = tags$i(class = "fa fa-square-o", 
                      style = "color: steelblue"))
      ),
      
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
        class = "btn-warning",
        style = 'width: 40%; margin-left: 30%; margin-right: 30%'
      )
    ),
    width = 350
  ),
  
  dashboardBody(
    useSweetAlert(),
    # add_busy_spinner(
    #   spin = "fading-circle",
    #   position = 'top-left',
    #   margins =  c('30%', '50%'),
    #   height = '70px',
    #   width = '70px'
    # ),
    add_busy_bar(
      timeout = 1200,
      color = '#FF0000',
      height = '5px'
    ),
    
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
      box(
        title = 'Heat Map',
        solidHeader = TRUE,
        status = 'warning',
        width = 12,
        # height = 500,
        textOutput('textHeatMap'),
        plotOutput("heatmap", width = "100%"),
        getDownloadButton('downloadHeatMap')
      ),
      # box(
      #   title = 'Violin Plot',
      #   solidHeader = TRUE,
      #   status = 'warning',
      #   width = 6,
      #   # height = 500,
      #   textOutput('textViolin'),
      #   plotOutput("violin", width = "100%"),
      #   getDownloadButton('downloadViolin')
      # )
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
server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2)
  
  output$fugoLogo <- renderImage({
    
    list(
      src = "assets/Logo.png",
      contentType = "image/png",
      width = '300px',
      # height = '100%',
      alt = "Fugo"
    )
    
  }, deleteFile = FALSE)
  
  observeEvent(input$applyButton,{

    file <- input$file1
    if (is.null(file)) {
      sendWarningAlert(session, 'No File Provided!', 'Please upload a data file to continue.')
    }

    ext <- tools::file_ext(file$datapath)
    req(file)
    if (ext != "csv") {
      sendWarningAlert(session, 'Invalid File Type!', 'Please upload a csv file.')
      return(NULL)
    }
    validate(need(ext == "csv", "Please upload a csv file"))
    
    data <- read.csv(file$datapath, header = TRUE)
    plotData <- data
    pThresh <- input$pvalueNum
    threshold <- input$thresholdNum
    foldValue <- input$foldNum
    referenceMarker <- input$referencePicker
    
    statSelection <- input$statAnalysisSelect

    if('1' %in% statSelection)
      plotData <- PValueCalculation(data, pThresh)
    if('2' %in% statSelection)
      plotData <- ThresholdCalculation(data, threshold)
    if('3' %in% statSelection)
      plotData <- FoldCalculation(data, foldValue, referenceMarker)
    
    output$contents <- DT::renderDataTable({
      plotData
      # fugoStats(plotData, avgRep = T)
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
    
    underConstructionMessage <- 'Not currently available. Refer to FugoCMDScript for information on using command line functions to visualize data.'
    
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
