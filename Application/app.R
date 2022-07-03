# Fugu-MS App

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinybusy)

fuguMSScript <- paste(dirname(getwd()), 'R/fuguMS.R', sep = '/')
source(fuguMSScript)

# Global Variables ----
stats_choices <-
  c('P-Value' = 1,
    'Threshold' = 2,
    'Fold Change' = 3)

sidePanelWidth <- 360
logoFileLocation <- 'www/assets/Logo.png'

# Helper Functions ----
PValueCalculation <- function(dat, pThresh) {
  alpha_val <- 0.05
  pDat <- fuguStats(dat, pCalc = T, pThresh = pThresh)
  num_samples = nrow(pDat)
  idx_p <- pDat$pVal < alpha_val / (num_samples)
  pDat <- pDat[idx_p,]
  
  # Remove NA data from frame
  pDat <- na.omit(pDat)
  
  return(pDat)
}

ThresholdCalculation <- function(dat, threshold) {
  # Calculate the average threshold per species
  mDat <- fuguStats(dat, avgRep = T)
  
  # Remove markers if the max average across all groups < given threshold
  idx_num <- apply(abs(metaSep(mDat)$data), 1, max) > threshold
  
  ## Find data with correct intensity and fold change
  idx <- which((idx_num) == 1)
  out <- mDat[idx,]
  
  return(out)
}

FoldCalculation <- function(dat, fold, reference) {
  # Calculate the average threshold per species
  mDat <- fuguStats(dat, avgRep = T)
  
  # filter markers from data set if the threshold is less than
  # 4-fold change in comparison to MHB
  # Select MHB or 1 as reference when calling this function
  fDat <-
    fuguStats(
      mDat,
      scale = 'fold',
      isFuguApp = TRUE,
      selectedReferenceList = reference
    )
  idx_fold <- (apply(abs(metaSep(fDat)$data), 1, max) > fold)
  
  ## Find data with correct intensity and fold change
  idx <- which((idx_fold) == 1)
  out <- mDat[idx,]
  
  return(out)
}

getDownloadButton <- function(inputId) {
  downloadButton(
    outputId = inputId,
    label = 'Download',
    icon = icon("download"),
    style = 'float: right;'
  )
}

sendWarningAlert <-
  function(session,
           headerText,
           messageText,
           isSuccess = FALSE) {
    alertType <- 'error'
    if (isSuccess) {
      alertType <- 'success'
    }
    
    sendSweetAlert(
      session = session,
      title = headerText,
      text = messageText,
      type = alertType,
      btn_labels = "Okay"
    )
  }

# Define UI for app ----
ui <- dashboardPage(
  skin = 'yellow',
  
  # App title
  dashboardHeader(title = "FUGU-MS",
                  titleWidth = sidePanelWidth),
  
  dashboardSidebar(
    sidebarMenu(
      div(
        style = 'padding-top: 10px',
        div(style = 'padding: 20px',
            imageOutput(outputId = "fuguLogo",
                        height = '200px' )),
        fileInput(
          inputId = "inputFile",
          label = "Upload Data File",
          width = '100%',
          multiple = FALSE,
          accept = c(".csv",
                     "text/csv",
                     "text/comma-separated-values,text/plain")
        )
      ),
      
      # hr(),
      fluidRow(h3('Statistical Analysis'),
               style = 'text-align: center; width: auto;'),
      
      div(
        style = 'padding: 0 15px',
        fluidRow(
          column(
            width = 6,
            style = 'padding: 0;',
            numericInput(
              'pvalueNum',
              'P-Value',
              1,
              min = 0,
              max = 1,
              step = 0.1,
              width = '100%'
            ),
          ),
          column(
            width = 6,
            style = 'padding: 0;',
            numericInput(
              'thresholdNum',
              'Threshold',
              20000,
              step = 10000,
              width = '100%'
            ),
          )
        ),
        fluidRow(
          column(
            width = 6,
            style = 'padding: 0;',
            numericInput(
              'foldNum',
              'Fold Change',
              4,
              min = 1,
              max = 100,
              width = '100%'
            ),
          ),
          column(
            width = 6,
            style = 'padding: 0;',
            # selectInput(
            #   inputId = 'referencePicker',
            #   label = 'Reference',
            #   choices = NULL,
            #   width = '100%',
            #   multiple = T,
            #   selected = 1,
            #   # options = list(style = "btn-primary")
            # ),
            pickerInput(
              inputId = "referencePicker",
              label = "Reference",
              choices = NULL,
              multiple = TRUE,
              width = '100%',
              options = list(
                `none-selected-text` = "Upload data first..."
              )
            ),
          )
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
                      style = "color: steelblue")
        )
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
    width = sidePanelWidth
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "apple-touch-icon", sizes = "180x180", href = "/apple-touch-icon.png"),
      tags$link(rel = "icon", type = "image/png", sizes = "32x32", href = "/favicon-32x32.png"),
      tags$link(rel = "icon", type = "image/png", sizes = "16x16", href = "/favicon-16x16.png"),
      tags$link(rel = "manifest", href = "/site.webmanifest"),
      tags$link(rel = "mask-icon", href = "/safari-pinned-tab.svg", color="#5bbad5"),
      tags$meta(name = "msapplication-TileColor", content = "#da532c"),
      tags$meta(name = "theme-color", content = "#ffffff")
    ),
    # tags$script(HTML("$('body').addClass('fixed');")),

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
      color = '#337AB7',
      height = '5px'
    ),
    
    fluidRow(
      box(
        title = 'Data Table',
        solidHeader = TRUE,
        status = 'warning',
        width = 12,
        # height = 500,
        div(DT::dataTableOutput('contents', width = '100%'),
            style = 'overflow-y: scroll;overflow-x: scroll;'),
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
        width = 6,
        # height = 500,
        # textOutput('textHeatMap'),
        plotOutput("heatmap", width = "100%"),
        getDownloadButton('downloadHeatMap')
      ),
      box(
        title = 'Violin Plot',
        solidHeader = TRUE,
        status = 'warning',
        width = 6,
        # height = 500,
        # textOutput('textViolin'),
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
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  
  output$fuguLogo <- renderImage({
    list(
      src = logoFileLocation,
      contentType = "image/png",
      width = '300px',
      # height = '100%',
      alt = "Fugu"
    )
    
  }, deleteFile = FALSE)
  
  values <- reactiveValues(data = NULL)
  
  observeEvent(input$inputFile, {
    file <- input$inputFile
    req(file)
    
    ext <- tools::file_ext(file$datapath)
    if (ext != "csv") {
      sendWarningAlert(session, 'Invalid File Type!', 'Please upload a csv file.')
      return(NULL)
    }
    
    values$data <- fuguRead(file$datapath)
  })
  
  choicesReferenceList <- reactive({
    req(values$data)
    by = "\\."
    sam <- names(values$data)
    sLab <- sapply(strsplit(sam, split = by), function (x)
      x[1])
    sLab <- sLab[-c(1, 2, 3)]
    unique(sLab)
  })
  
  observeEvent(choicesReferenceList(), {
    updateSelectInput(
      session = session,
      inputId = 'referencePicker',
      choices = choicesReferenceList(),
      selected = choicesReferenceList()[1]
    )
    updatePickerInput(
      session = session,
      inputId = 'referencePicker',
      selected = choicesReferenceList()[1],
      choices = choicesReferenceList(),
      options = list(
        `actions-box` = TRUE,
        # `deselect-all-text` = "Deselect All",
        # `select-all-text` = "Select All",
        `none-selected-text` = "Select marker...",
        `selected-text-format` = "count > 3",
        `size` = 5
      ),
      clearOptions = TRUE
    )
  })
  
  observeEvent(input$applyButton, {
    if (is.null(values$data)) {
      sendWarningAlert(session,
                       'No File Provided',
                       'Please upload a data file to continue.')
      return(NULL)
    }
    
    plotData <- values$data
    pThresh <- input$pvalueNum
    threshold <- input$thresholdNum
    foldValue <- input$foldNum
    referenceMarker <- input$referencePicker
    
    statSelection <- input$statAnalysisSelect
    
    if ('1' %in% statSelection)
      plotData <- PValueCalculation(plotData, pThresh)
    if ('2' %in% statSelection)
      plotData <- ThresholdCalculation(plotData, threshold)
    if ('3' %in% statSelection) {
      if (is.null(referenceMarker)) {
        sendWarningAlert(session,
                         'No Marker Provided',
                         'Please select a marker from the Reference drop-down to perform a Fold Change analysis.')
        return(NULL)
      }
        
      plotData <- FoldCalculation(plotData, foldValue, referenceMarker)
    }
      
    
    output$contents <- DT::renderDataTable({
      # plotData
      fuguStats(plotData, avgRep = T)
    })
    
    output$heatmap <- renderPlot({
      fuguPlot(
        plotData,
        heatMap = T,
        isFuguApp = TRUE,
        scale = 'row',
        rCst = T,
        # grid = F,
        # cCst = F
      )
    },
    res = 96)
    
    output$violin <- renderPlot({
      # fuguPlot(plotData,vioPlot = T)
    },
    res = 96)
    
    underConstructionMessage <-
      'Not currently available. Refer to FuguCMDScript for information on using command line functions to visualize data.'
    
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
      filename = "fugu_HeatMap.png",
      content = function(file) {
        png(file)
        print(fuguPlot(
          plotData,
          heatMap = T,
          scale = 'row',
          rCst = T,
          grid = F,
          cCst = F
        ))
        dev.off()
      }
    )
    
    output$downloadDataTable <- downloadHandler(
      filename = "fugu_HeatMap.csv",
      content = function(file) {
        write.csv(plotData, file)
      }
    )
  })
}

# Deploy app ----
shinyApp(ui = ui, server = server)
