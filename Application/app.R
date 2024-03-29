# Fugu-MS App

library(plotly)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinybusy)

# fuguMSScript <- paste(dirname(getwd()), "R/fuguMS.R", sep = "/")
fuguMSScript <- "FuguMSLibrary/fuguMS.R"
source(fuguMSScript)

# Global Variables ----
stats_choices <-
  c("P-Value" = 1,
    "Threshold" = 2,
    "Fold Change" = 3)

scaleChoices = c(
  "None" = "none",
  "Log" = "log",
  "Durbin" = "durbin",
  "Row" = "row",
  "Column" = "column",
  "Fold" = "fold",
  "Ratio" = "ratio",
  "ZScore" = "zScore"
)

sidePanelWidth <- 360
logoFileLocation <- "www/assets/Logo.png"

compoundErrorMessage <-
  "At least one compound must be selected to plot the"
sampleErrorMessage <-
  "At least one sample must be selected to plot the"

downloadButtonLabel <- "Download"
plotButtonLabel <- "Generate Plot"

# Helper Functions ----
PValueCalculation <- function(dat, pThresh) {
  alpha_val <- 0.05
  pDat <- fuguStats(dat, pCalc = T, pThresh = pThresh)
  num_samples = nrow(pDat)
  idx_p <- pDat$pVal < alpha_val / (num_samples)
  pDat <- pDat[idx_p, ]
  
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
  out <- mDat[idx, ]
  
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
      scale = "fold",
      isFuguApp = TRUE,
      selectedReferenceList = reference
    )
  idx_fold <- (apply(abs(metaSep(fDat)$data), 1, max) > fold)
  
  ## Find data with correct intensity and fold change
  idx <- which((idx_fold) == 1)
  out <- mDat[idx, ]
  
  return(out)
}

getDownloadButton <- function(inputId) {
  div(
    style = "padding-top: 4px",
    downloadButton(
      outputId = inputId,
      label = downloadButtonLabel,
      icon = icon("download"),
      style = "float: right;"
    )
  )
}

sendWarningAlert <-
  function(session,
           headerText,
           messageText,
           isSuccess = FALSE) {
    alertType <- "error"
    if (isSuccess) {
      alertType <- "success"
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
  skin = "yellow",
  
  # App title
  dashboardHeader(
    title = "FUGU-MS",
    titleWidth = sidePanelWidth,
    tags$li(
      div(
        tags$a(tags$i(class = "fab fa-github fa-2x"), href="https://github.com/LewisResearchGroup/fugu-ms-app", target="_blank"),
        tags$head(tags$style(HTML("a {color: black;}"))),
        style = "padding:11px; padding-right:15px;"
      ),
      class = "dropdown"
    )
  ),
  
  dashboardSidebar(
    sidebarMenu(
      div(
        style = "padding-top: 10px",
        div(style = "padding: 20px",
            tags$a(imageOutput(outputId = "fuguLogo", height = "200px"), href="https://github.com/LewisResearchGroup/fugu-ms-app", target="_blank")
            ),
        fileInput(
          inputId = "inputFile",
          label = "Upload Data File",
          width = "100%",
          multiple = FALSE,
          accept = c(".csv",
                     "text/csv",
                     "text/comma-separated-values,text/plain")
        )
      ),
      
      # Statistical Analysis
      fluidRow(
        h3("Statistical Analysis"),
        style = "text-align: center; width: auto;"
      ),
      
      div(
        style = "padding: 0 15px 20px",
        fluidRow(
          column(
            width = 6,
            style = "padding: 0;",
            numericInput(
              "pvalueNum",
              "P-Value",
              1,
              min = 0,
              max = 1,
              step = 0.1,
              width = "100%"
            )
          ),
          column(
            width = 6,
            style = "padding: 0;",
            numericInput(
              "thresholdNum",
              "Threshold",
              20000,
              step = 10000,
              width = "100%"
            )
          )
        ),
        fluidRow(
          column(
            width = 4,
            style = "padding: 0; margin: -15px 0",
            numericInput("foldNum",
                         "Fold Change",
                         4,
                         min = 1,
                         # max = 10000,
                         width = "100%")
          ),
          column(
            width = 8,
            style = "padding: 0; margin: -15px 0",
            pickerInput(
              inputId = "referencePicker",
              label = "Reference",
              choices = NULL,
              multiple = TRUE,
              width = "100%",
              options = list(`none-selected-text` = "Upload data first...")
            )
          )
        ),
      ),
      
      checkboxGroupButtons(
        inputId = "statAnalysisSelect",
        label = "Select analysis to apply:",
        choices = stats_choices,
        justified = TRUE,
        # individual = TRUE,
        checkIcon = list(
          yes = tags$i(class = "fa fa-check-square",
                       style = "color: orange"),
          no = tags$i(class = "fa fa-square-o",
                      style = "color: steelblue")
        )
      ),

      # Clustering
      br(),
      hr(),
      fluidRow(
        h3("Clustering"),
        style = "text-align: center; width: auto;"
      ),
      
      div (
        style = "padding: 0px 15px;",
        sliderInput("sliderMz", label = "Mass", min = 0, max = 1, step = 0.1, value = 1, tick = TRUE, width = "100%"),
        sliderInput("sliderRT", label = "Retention Time", min = 0, max = 1, step = 0.1, value = 1, tick = TRUE, width = "100%"),
        sliderInput("sliderCov", label = "Covariance", min = 0, max = 1, step = 0.1, value = 1, tick = TRUE, width = "100%"),
        sliderInput("sliderPthres", label = "P-Threshold", min = 0, max = 1, step = 0.1, value = 0.8, tick = TRUE, width = "100%")
      ),

      div(
        style = "padding-top: 4px; float: right;",
        materialSwitch(
          inputId = "clusteringToggle", 
          label = "Apply Clustering", 
          status = "warning",
          inline = TRUE
        ),
      ),

      br(),
      br(),
      hr(),
      # Button Group
      div (
        style = "padding: 40px 0;",
        # actionButton(
        #   inputId = "clearButton",
        #   label = "Reset",
        #   class = "btn-danger",
        #   style = "width: 40%; margin-left: 30%; margin-right: 30%"
        # ),
        
        actionButton(
          inputId = "applyButton",
          label = "Plot",
          class = "btn-warning",
          style = "width: 40%; margin-left: 30%; margin-right: 30%"
        )
      )
    ),
    width = sidePanelWidth
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "apple-touch-icon", sizes = "180x180", href = "/apple-touch-icon.png"),
      tags$link(
        rel = "icon",
        type = "image/png",
        sizes = "32x32",
        href = "/favicon-32x32.png"
      ),
      tags$link(
        rel = "icon",
        type = "image/png",
        sizes = "16x16",
        href = "/favicon-16x16.png"
      ),
      tags$link(rel = "manifest", href = "/site.webmanifest"),
      tags$link(rel = "mask-icon", href = "/safari-pinned-tab.svg", color =
                  "#5bbad5"),
      tags$meta(name = "msapplication-TileColor", content = "#da532c"),
      tags$meta(name = "theme-color", content = "#ffffff")
    ),
    # tags$script(HTML("$("body").addClass("fixed");")),
    
    useSweetAlert(),
    # add_busy_spinner(
    #   spin = "fading-circle",
    #   position = "top-left",
    #   margins =  c("30%", "50%"),
    #   height = "70px",
    #   width = "70px"
    # ),
    add_busy_bar(
      timeout = 1200,
      color = "#337AB7",
      height = "5px"
    ),
    
    fluidRow(
      box(
        title = "Data Table",
        solidHeader = TRUE,
        status = "warning",
        width = 12,
        # height = 500,
        collapsible = TRUE,
        div(DT::dataTableOutput("contents", width = "100%"),
            style = "overflow-y: scroll;overflow-x: scroll;"),
        br(),
        getDownloadButton("downloadDataTable")
      )
    ),
    
    # Dot, Bar Plot
    fluidRow(
      div(
        class = "col-sm-12 col-md-12 col-lg-6",
        box(
          title = "Heat Map",
          solidHeader = TRUE,
          status = "warning",
          width = NULL,
          collapsible = TRUE,
          # height = 500,
          fluidRow(
            column(
              width = 6,
              pickerInput(
                inputId = "heatmapScalePicker",
                label = "Scale",
                choices = scaleChoices,
                multiple = FALSE,
                width = "100%"
              )
            ),
            column(
              width = 6,
              ## If fold, ratio, or zScore, scale by:
              conditionalPanel(
                condition = 'input.heatmapScalePicker == "fold" | input.heatmapScalePicker == "ratio" | input.heatmapScalePicker == "zScore"',
                pickerInput(
                  inputId = "heatmapScaleByPicker",
                  label = "Samples to Scale By",
                  choices = NULL,
                  multiple = TRUE,
                  width = "100%",
                  options = list(`none-selected-text` = "Upload data first...")
                )
              )
            )
          ),
          actionButton(
            inputId = "heatmapApplyButton",
            label = plotButtonLabel,
            class = "btn-warning",
            # width = "100%"
            style = "float: right;"
          ),
          fluidRow(column(
            width = 12,
            plotOutput("heatmap", width = "100%")
          )),
          getDownloadButton("downloadHeatMap")
        )
      ),
      div(
        class = "col-sm-12 col-md-12 col-lg-6",
        box(
          title = "Violin Plot",
          solidHeader = TRUE,
          status = "warning",
          width = NULL,
          collapsible = TRUE,
          # height = 500,
          fluidRow(column(
            width = 6,
            pickerInput(
              inputId = "violinSamplePicker",
              label = "Sample Groups",
              choices = NULL,
              multiple = TRUE,
              width = "100%",
              options = list(`none-selected-text` = "Upload data first...")
            )
          ),
          column(
            width = 6,
            pickerInput(
              inputId = "violinCompoundPicker",
              label = "Compounds",
              choices = NULL,
              multiple = FALSE,
              width = "100%",
              options = list(`none-selected-text` = "Upload data first...")
            )
          )),
          actionButton(
            inputId = "violinApplyButton",
            label = plotButtonLabel,
            class = "btn-warning",
            # width = "100%"
            style = "float: right;"
          ),
          fluidRow(column(
            width = 12,
            plotOutput("violin", width = "100%")
          )),
          getDownloadButton("downloadViolin")
        )
      )
    ),
    
    fluidRow(
      div(
        class = "col-sm-12 col-md-12 col-lg-6",
        box(
          title = "Bar Plot",
          solidHeader = TRUE,
          status = "warning",
          width = NULL,
          collapsible = TRUE,
          # height = 500,
          fluidRow(column(
            width = 6,
            pickerInput(
              inputId = "barSamplePicker",
              label = "Sample Groups",
              choices = NULL,
              multiple = TRUE,
              width = "100%",
              options = list(`none-selected-text` = "Upload data first...")
            )
          ),
          column(
            width = 6,
            pickerInput(
              inputId = "barCompoundPicker",
              label = "Compounds",
              choices = NULL,
              multiple = FALSE,
              width = "100%",
              options = list(`none-selected-text` = "Upload data first...")
            )
          )),
          actionButton(
            inputId = "barApplyButton",
            label = plotButtonLabel,
            class = "btn-warning",
            # width = "100%"
            # style = "width: 40%; margin-left: 30%; margin-right: 30%"
            style = "float: right;"
          ),
          fluidRow(column(
            width = 12,
            plotOutput("barplot", width = "100%")
          )),
          getDownloadButton("downloadBar")
        )
      ),
      div(
        class = "col-sm-12 col-md-12 col-lg-6",
        box(
          title = "Dot Plot",
          solidHeader = TRUE,
          status = "warning",
          width = NULL,
          collapsible = TRUE,
          # height = 500,
          fluidRow(column(
            width = 6,
            pickerInput(
              inputId = "dotSamplePicker",
              label = "Sample Groups",
              choices = NULL,
              multiple = TRUE,
              width = "100%",
              options = list(`none-selected-text` = "Upload data first...")
            ),
          ),
          column(
            width = 6,
            pickerInput(
              inputId = "dotCompoundPicker",
              label = "Compounds",
              choices = NULL,
              multiple = FALSE,
              width = "100%",
              options = list(`none-selected-text` = "Upload data first...")
            )
          )),
          actionButton(
            inputId = "dotApplyButton",
            label = plotButtonLabel,
            class = "btn-warning",
            # width = "100%"
            style = "float: right;"
          ),
          fluidRow(column(width = 12,
                          plotOutput("dot", width = "100%"))),
          getDownloadButton("downloadDot")
        )
      )
    ),
    
    fluidRow(
      div(
        class = "col-sm-12 col-md-12 col-lg-6",
        box(
          title = "PCA Plot",
          solidHeader = TRUE,
          status = "warning",
          width = NULL,
          collapsible = TRUE,
          fluidRow(
            column(
              width = 12,
              plotlyOutput("pca")
              ),
            ),
        )
      )
    )
  )
)


# Server Helper Functions ----

updateUIList <- function(session, id, updatedList) {
  updatePickerInput(
    session = session,
    inputId = id,
    selected = updatedList,
    choices = updatedList,
    options = list(
      `actions-box` = TRUE,
      # `deselect-all-text` = "Deselect All",
      # `select-all-text` = "Select All",
      `none-selected-text` = "Select reference marker(s)...",
      `selected-text-format` = "count > 3",
      `size` = 5
    ),
    clearOptions = TRUE
  )
}

processData <- function(input, output, session, values) {
  verifyDataFile(session, values$data)
  
  plotData <- values$data
  pThresh <- input$pvalueNum
  threshold <- input$thresholdNum
  foldValue <- input$foldNum
  referenceMarker <- input$referencePicker

  mz <- input$sliderMz
  rt <- input$sliderRT
  cov <- input$sliderCov
  pThreshCluster <- input$sliderPthres
  
  statSelection <- input$statAnalysisSelect
  if ("1" %in% statSelection)
    plotData <- PValueCalculation(plotData, pThresh)
  if ("2" %in% statSelection)
    plotData <- ThresholdCalculation(plotData, threshold)
  if ("3" %in% statSelection) {
    if (is.null(referenceMarker)) {
      sendWarningAlert(
        session,
        "No Marker Provided",
        "Please select a marker from the Reference drop-down to perform a Fold Change analysis."
      )
      return(NULL)
    }
    plotData <-
      FoldCalculation(plotData, foldValue, referenceMarker)
  }

  if (input$clusteringToggle) {
    mz_ad <- c(0,
           0.0005,        #(e-)
           1.0078,        #(H+)         https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           1.0062,        #(H neutron), https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           1.0034,        #(C neutron), https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           0.9694,        #(N neutron), https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           18.0106,       #(H2O)        https://hmdb.ca/metabolites/HMDB0002111
           34.9689,       #(Cl)         https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           38.9637,       # K           https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           22.9898,       #(Na)         https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           44.9977,       #(Formate),   https://hmdb.ca/metabolites/HMDB0304356
           1.0042,        #(O neutrons) https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           0.9994,        #(S neutrons) https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           1.9958,        #(S neutrons) https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
           3.9950         #(S neutrons) https://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
    )

    cDat <-
      fuguCluster(
        plotData,
        pW = c(mz, rt, cov), #### Mass, Retention Time (0,1), Covariance
        pThresh = pThreshCluster,##### 
        adList = mz_ad, #####
        collapse = FALSE,
        clusterOnly = FALSE,
      )
    plotData <- cDat[!duplicated(cDat$grp), ]
  }
  
  return(plotData)
}

verifyDataFile <- function(session, data) {
  if (is.null(data)) {
    sendWarningAlert(session,
                     "No File Provided",
                     "Please upload a data file to continue.")
    return(NULL)
  }
}

verifyInput <- function(inputParameter, errorMessage, session) {
  if (is.null(inputParameter)) {
    sendWarningAlert(session,
                     "Missing field(s)",
                     errorMessage)
    return(FALSE)
  } else {
    return (TRUE)
  }
}

plotHeatMap <- function(input, output, session, plotData) {
  rowCluster <- TRUE
  
  heatmapScale <- input$heatmapScalePicker
  sampleSelected <- input$heatmapScaleByPicker
  heatmapErrorMessage <-
    paste(sampleErrorMessage, "heat map.", sep = " ")
  if (!verifyInput(sampleSelected, heatmapErrorMessage, session))
    return(FALSE)
  
  output$heatmap <- renderPlot({
    fuguPlot(
      plotData,
      heatMap = T,
      isFuguApp = TRUE,
      scale = heatmapScale,
      shiny_grp_list = sampleSelected,
      rCst = rowCluster
      # cCst = F
    )
  },
  res = 96)
  
  output$downloadHeatMap <- downloadHandler(
    filename = "fugu_HeatMap.png",
    content = function(file) {
      png(file)
      print(
        fuguPlot(
          plotData,
          heatMap = T,
          isFuguApp = TRUE,
          scale = heatmapScale,
          shiny_grp_list = sampleSelected,
          rCst = rowCluster
          # cCst = F
        )
      )
      dev.off()
    }
  )
  return(TRUE)
}

plotBar <- function(input, output, session, plotData) {
  sampleSelected <- input$barSamplePicker
  compoundSelected <- input$barCompoundPicker
  
  barFlag <- "bar graph."
  barSampleErrorMessage <-
    paste(sampleErrorMessage, barFlag, sep = " ")
  barCompoundErrorMessage <-
    paste(compoundErrorMessage, barFlag, sep = " ")
  
  if (!verifyInput(sampleSelected, barSampleErrorMessage, session) ||
      !verifyInput(compoundSelected, barCompoundErrorMessage, session))
    return(FALSE)
  
  output$barplot <- renderPlot({
    fuguPlot(
      plotData,
      BoxWhisker = T,
      isFuguApp = T,
      shiny_grp_list = sampleSelected,
      shiny_cmp_list = compoundSelected
    )
  },
  res = 96)
  
  output$downloadBar <- downloadHandler(
    filename = "fugu_BarPlot.png",
    content = function(file) {
      png(file)
      print(
        fuguPlot(
          plotData,
          BoxWhisker = T,
          isFuguApp = T,
          shiny_grp_list = sampleSelected,
          shiny_cmp_list = compoundSelected
        )
      )
      dev.off()
    }
  )
  return(TRUE)
}

plotViolin <- function(input, output, session, plotData) {
  sampleSelected <- input$violinSamplePicker
  compoundSelected <- input$violinCompoundPicker
  
  violinFlag <- "violin graph."
  violinSampleErrorMessage <-
    paste(sampleErrorMessage, violinFlag, sep = " ")
  violinCompoundErrorMessage <-
    paste(compoundErrorMessage, violinFlag, sep = " ")
  
  if (!verifyInput(sampleSelected, violinSampleErrorMessage, session) ||
      !verifyInput(compoundSelected, violinCompoundErrorMessage, session))
    return(FALSE)
  
  output$violin <- renderPlot({
    fuguPlot(
      plotData,
      vioPlot = T,
      isFuguApp = T,
      shiny_grp_list = sampleSelected,
      shiny_cmp_list = compoundSelected
    )
  },
  res = 96)
  
  output$downloadViolin <- downloadHandler(
    filename = "fugu_ViolinPlot.png",
    content = function(file) {
      png(file)
      print(
        fuguPlot(
          plotData,
          vioPlot = T,
          isFuguApp = T,
          shiny_grp_list = sampleSelected,
          shiny_cmp_list = compoundSelected
        )
      )
      dev.off()
    }
  )
  return(TRUE)
}

plotDot <- function(input, output, session, plotData) {
  sampleSelected <- input$dotSamplePicker
  compoundSelected <- input$dotCompoundPicker
  
  dotFlag <- "dot plot."
  dotSampleErrorMessage <-
    paste(sampleErrorMessage, dotFlag, sep = " ")
  dotCompoundErrorMessage <-
    paste(compoundErrorMessage, dotFlag, sep = " ")
  
  if (!verifyInput(sampleSelected, dotSampleErrorMessage, session) ||
      !verifyInput(compoundSelected, dotCompoundErrorMessage, session))
    return(FALSE)
  
  output$dot <- renderPlot({
    fuguPlot(
      plotData,
      dotPlot = T,
      isFuguApp = T,
      shiny_grp_list = sampleSelected,
      shiny_cmp_list = compoundSelected
    )
  },
  res = 96)
  
  output$downloadDot <- downloadHandler(
    filename = "fugu_DotPlot.png",
    content = function(file) {
      png(file)
      print(
        fuguPlot(
          plotData,
          dotPlot = T,
          isFuguApp = T,
          shiny_grp_list = sampleSelected,
          shiny_cmp_list = compoundSelected
        )
      )
      dev.off()
    }
  )
  return(TRUE)
}

plotPca <- function(input, output, session, plotData) {
  # p <- plot_ly(iris, x = ~get(input$choice), y = ~Sepal.Length, type = 'scatter', mode = 'markers')
  p <- fuguPlot(plotData, pcaPlot = T,pcaLabels = T,pcaVectors = T,pcaNumVectors= 100,pcaVectorNames = F)
  output$pca <- renderPlotly({
    p
  })

  return(TRUE)
}

# Define Server for app ----
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)

  output$fuguLogo <- renderImage({
    list(
      src = logoFileLocation,
      contentType = "image/png",
      width = "300px",
      # height = "100%",
      alt = "Fugu"
    )
    
  }, deleteFile = FALSE)
  
  values <- reactiveValues(data = NULL, plotData = NULL)
  
  observeEvent(input$inputFile, {
    file <- input$inputFile
    req(file)
    
    ext <- tools::file_ext(file$datapath)
    if (ext != "csv") {
      sendWarningAlert(session, "Invalid File Type!", "Please upload a csv file.")
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
    updateUIList(session, "referencePicker", choicesReferenceList())
    updateUIList(session, "barSamplePicker", choicesReferenceList())
    updateUIList(session, "barCompoundPicker", values$data$compound)
    updateUIList(session, "dotSamplePicker", choicesReferenceList())
    updateUIList(session, "dotCompoundPicker", values$data$compound)
    updateUIList(session, "violinSamplePicker", choicesReferenceList())
    updateUIList(session, "violinCompoundPicker", values$data$compound)
    updateUIList(session, "heatmapScaleByPicker", choicesReferenceList())
  })
  
  observeEvent(input$applyButton, {
    values$plotData <- processData(input, output, session, values)
    
    output$contents <- DT::renderDataTable({
      fuguStats(values$plotData, avgRep = T)
    })
    
    output$downloadDataTable <- downloadHandler(
      filename = "fugu_DataTable.csv",
      content = function(file) {
        write.csv(values$plotData, file)
      }
    )
    
    # Heat Map
    isSuccess <- plotHeatMap(input, output, session, values$plotData)
    if (!isSuccess)
      return (NULL)
    
    # Violin Plot
    isSuccess <- plotViolin(input, output, session, values$plotData)
    if (!isSuccess)
      return (NULL)
    
    # Bar Plot
    isSuccess <- plotBar(input, output, session, values$plotData)
    if (!isSuccess)
      return (NULL)
    
    # Dot Plot
    isSuccess <- plotDot(input, output, session, values$plotData)
    if (!isSuccess)
      return (NULL)
    
    # Pca Plot
    isSuccess <- plotPca(input, output, session, values$plotData)
    if (!isSuccess)
      return (NULL)
    
  })
  
  observeEvent(input$heatmapApplyButton, {
    verifyDataFile (session, values$plotData)
    
    isSuccess <- plotHeatMap(input, output, session, values$plotData)
    if (!isSuccess)
      return (NULL)
  })
  
  observeEvent(input$barApplyButton, {
    verifyDataFile (session, values$plotData)
    
    isSuccess <- plotBar(input, output, session, values$plotData)
    if (!isSuccess)
      return (NULL)
  })
  
  observeEvent(input$dotApplyButton, {
    verifyDataFile (session, values$plotData)
    
    isSuccess <- plotDot(input, output, session, values$plotData)
    if (!isSuccess)
      return (NULL)
  })
  
  observeEvent(input$violinApplyButton, {
    verifyDataFile (session, values$plotData)
    
    isSuccess <- plotViolin(input, output, session, values$plotData)
    if (!isSuccess)
      return (NULL)
  })
}

# Deploy app ----
shinyApp(ui = ui, server = server)
