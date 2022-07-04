################################################################################
##  FUGU-MS Tools for processing and visualization of mass spectroscpy files  ##
################################################################################
##    Original Author: Ian Lewis											                        ##
##    Last Edit: Raied Aburashed, Saad Luqman                                 ##
## 	  Version: 1.1.0														                              ##
##   																		                                      ##
##    This program is free software: you can redistribute it and/or modify    ##
##    it under the terms of the GNU General Public License as published by    ##
##    the Free Software Foundation, either version 3 of the License, or       ##
##    any later version.                                                      ##
##                                                                            ##
##    This program is distributed in the hope that it will be useful,         ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of          ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           `##
##    GNU General Public License for more details.                            ##
##                                                                            ##
##    A copy of the GNU General Public License can be found at:               ##
##    www.r-project.org/Licenses/GPL-3                                        ##
##                                                                            ##
################################################################################

################################################################################
##                                                                            ##
##                  Simplified Reading and writing files                      ##
##                                                                            ##
################################################################################
# Read Files
fuguRead <- function(inPath = file.choose()) {
  return(read.csv(inPath, head = TRUE, stringsAsFactors = FALSE))
}

# Write Files
fuguWrite <- function(table) {
  write.csv(table,
            file.choose(new = TRUE),
            quote = FALSE,
            row.names = FALSE)
}

################################################################################
isInputListValid <- function(inDat) {
  ## check format of peak list
  if (!is.data.frame(inDat) || nrow(inDat) < 1 || ncol(inDat) < 2 || !any(names(inDat) == 'compound'))
    stop("A dataframe with a 'compound' and at least one data column is required")

  ## Check if in compound list there is an NA
  com <- inDat$compound
  if (any(is.na(com)) &&  any(names(inDat) == 'medMz') && any(names(inDat) == 'medRt')) {
    com[which(is.na(com))] <- paste(inDat$medMz[which(is.na(com))],
                                    inDat$medRt[which(is.na(com))], sep = "@")
  }
  
  for (i in unique(com)) {
    tCpd <- com[com == i]
    if (length(tCpd) == 1)
      next
    com[com == i] <- paste(tCpd, 1:length(tCpd), sep = '_')
  }
  inDat$compound <- com
  
  invisible(inDat)
}

## Internal function for separating metadata from data sets
## excludeName - character - names of columns to be excluded from analysis
## meta - character string - default names of Maven format metadata
metaSep <- function(inDat,
                    meta = c(
                      'label',
                      'metaGroupId',
                      'groupId',
                      'goodPeakCount',
                      'medMz',
                      'medRt',
                      'maxQuality',
                      'note',
                      'compound',
                      'compoundId',
                      'expectedRtDiff',
                      'ppmDiff',
                      'parent',
                      'adductName',
                      'isotopeLabel',
                      'formula',
                      'ctName',
                      'pVal'
                    ),
                    excludeName) {
  ## Check input data set
  if (missing(inDat) || !is.data.frame(inDat) || (nrow(inDat) < 1))
    stop('A METEN format dataset is required')
  ## Combine default and user-provided names
  if (!missing(excludeName))
    meta <- c(meta, excludeName)
  ## Remove all standard METen meta data
  mIdx <- na.omit(match(meta, names(inDat)))
  ## Find metadaFta
  if (length(mIdx) == 1) {
    metaDat <- data.frame(compound = inDat[, mIdx], stringsAsFactors = F)
  } else
    metaDat <- inDat[, mIdx]
  ## remove the metaData
  outDat <- inDat[, -na.omit(match(meta, names(inDat)))]
  return(list(metaData = metaDat, data = outDat))
}

print_app <- function(widget) {
  # Generate random file name
  temp <- paste(tempfile('plotly'), 'html', sep = '.')
  
  # Save. Note, leaving selfcontained=TRUE created files that froze my browser
  htmlwidgets::saveWidget(widget, temp, selfcontained = FALSE)
  
  # Launch with desired application
  if (.Platform$OS.type == "windows") {
    system(sprintf("chromium-browser -app=file://%s", temp))
  }
  else{
    system(sprintf("open -a 'safari'  /%s", temp))
  }
  
  # Return file name if it's needed for any other purpose
  temp
}

################################################################################
########################## Math Operations #####################################
################################################################################
## Log transform function
## Bioinformatics: Rocke and Durbin 2003
durbin = function(y, c = 0) {
  if (!is.null(c))
    y = log((y + sqrt(y ^ 2 + c ^ 2)) / 2)
  tmp = !is.finite(as.matrix(y))
  if (sum(tmp))
    y[tmp] = NA
  y
}

foldChange <-
  function(inTable,
           normVec,
           zeroOffset = F,
           normMin = .01,
           scale = 'fold') {
    ## Calculate standard deviation for normalization vector
    if (scale == "zScore") {
      if (length(normVec) < 1)
        stop('Z scales require a standard deviation from multiple replicates')
      sdTab <- apply(inTable[, normVec], 1, mySd)
      normTab <- apply(inTable[, normVec], 1, myMean)
    }
    else{
      if (length(normVec) > 1)
        normTab <- apply(inTable[, normVec], 1, myMean)
      else
        normTab <- inTable[, normVec]
    }
    
    ## Make norm matrix
    normTab <- as.vector(unlist(normTab))
    normTab <- data.frame(matrix(rep(normTab, ncol(inTable)),
                                 ncol = ncol(inTable)))
    
    ## Make sd matrix
    if (scale == "zScore") {
      sdTab <- as.vector(unlist(sdTab))
      sdTab <-
        data.frame(matrix(rep(sdTab, ncol(inTable)), ncol = ncol(inTable)))
    }
    
    ## Set all 0 values to 1% of min observed value
    if (zeroOffset) {
      zSet <- inTable == 0
      if (mySum(zSet)) {
        tTab <- inTable
        tTab[zSet] <- NA
        mSet <- as.vector(unlist(apply(tTab, 1, myMin)))
        for (i in 1:length(mSet)) {
          if (any(inTable[i, ] == 0, na.rm = T)) {
            idx <- which(inTable[i, ] == 0)
            inTable[i, idx] <- mSet[i]
          }
        }
      }
    }
    
    ## Calculate z score
    if (scale == 'zScore') {
      out <- (inTable - normTab) / sdTab
      #out <- apply( out, 2, myNaN )
      return(out)
    }
    
    
    ## Normalize
    out <- myDiv(inTable,  normTab)
    
    if (scale == 'ratio') {
      out <- apply(out, 2, myNaN)
      return(out)
    }
    
    if (any(out < 1)) {
      foldFun <- function(x) {
        tIdx <- which(x < 1)
        x[tIdx] <- abs(1 / x[tIdx]) * -1
        return(x)
      }
      out <- apply(out, 2, foldFun)
      out <- apply(out, 2, myNaN)
    }
    
    return(out)
    
  }

## Function for dividing by zero
## x - dataframe with numerator data
## y - dataframe with denominator data
## num - numeric, the value to be assigned for non-zero numerators dived by 0
## den - numeric, the value to be assigned for 0/0
## returns dataframe of dim(x)
myDiv <- function(x, y, num = 0, den = 1) {
  dX <- dim(x)
  dY <- dim(y)
  if (all(dX != dY))
    stop("x and y must have the same number of rows and columns")
  vX <- as.numeric(as.vector(unlist(x)))
  vY <- as.numeric(as.vector(unlist(y)))
  if (any(vY == 0))
    vY[vY == 0] <- den
  out <- data.frame(matrix(vX / vY, ncol = ncol(x)))
  names(out) <- names(x)
  return(out)
}
## Internal function for computing MAX
myMax <- function(x) {
  max(x, na.rm = T)
}
## Internal function for computing MIN
myMin <- function(x) {
  min(x, na.rm = TRUE) * normMin
}
## Internal function for computing SUM of a vector
mySum <- function(x) {
  sum(x, na.rm = TRUE)
}
## Internal function for computing Standard Deviation
mySd <- function(x) {
  sd(x, na.rm = TRUE)
}
## Internal function for computing MEAN
myMean <- function(x) {
  mean(x, na.rm = TRUE)
}
## Make NaN values zeros
myNaN <- function(x) {
  if (any(is.nan(x)))
    x[is.nan(x)] <- 0
  return(x)
}
## Internal function to find where the max iis indexed
myWhich.max <-
  function(x) {
    out <- myMax(x)
    return(which(x == out)[1])
  }
# return a selection from user
mySelect <- function(inVector, index = T) {
  uSel <- select.list(inVector, multiple = T)
  if (!index)
    return(uSel)
  
  return(match(uSel, inVector))
}
## Internal function for finding a mode of a dataset; for multimode data
## this function only returns the first index
myMode <- function(x) {
  ux <- sort(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}
## Utility function for matching multiples with multiples
## x input list || y index list
## returns match index in units of Y
myMatch <- function(x, y) {
  outI <- which(!is.na(match(y, x)))
  out <- rep(NA, length(y))
  uX <- unique(x)
  for (i in 1:length(uX)) {
    if (is.na(uX[i]))
      next
    if (any(y == uX[i]))
      out[y == uX[i]] <- i
  }
  
  if (any(is.na(out)))
    out <- out[-which(is.na(out))]
  return(outI[order(out)])
}
## Internal function for splitting text vectors
mySplit <- function(inDat, by = '\\.', index = 1) {
  tDat <- strsplit(inDat, split = by)
  
  for (i in 1:length(tDat)) {
    if (length(tDat[[i]]) < 1) {
      tDat[i] <- NA
    }
    if (length(tDat[[i]]) > 1) {
      tDat[i] <- tDat[[i]][index]
    }
    if (length(tDat[[i]]) == 1) {
      tDat[i] <- tDat[[i]]
    }
    
  }
  return(unlist(tDat))
}
## Standard error function that removes NAs
stdErr <- function(x) {
  std.err <- function(y) {
    if (any(is.na(y)))
      y <- y[-which(is.na(y))]
    return(sqrt(var(y) / length(y)))
  }
  if (is.null(nrow(x)) || nrow(x) < 2)
    return(std.err(x))
  return(apply(x, 2, std.err))
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

## This function takes a METen-format list and makes several types of plots
## inDat        - el-maven format table
## excludeName  - character - names of columns to be excluded from analysis
## na.rm  - Logical - TRUE converts NA to 0, INF to positive max, -INF to negative max
## rmZero - logical - removes rows with all zero intensities
##
## rCst   - Logical - TRUE clusters by row
## cCst   - Logical - TRUE clusters by column
## scale  - Controls how data is normalized; must be one of the following: 
##          "log", "durbin", "row", "column", or "none"
##					  'log' applies base 2 log transform
##            'durbin' applies durbin transform
##            'row' normalizes by row
##            'column' normalizes by column
## 					  'zScore' normalizes relative to a control, sets equal variance
##            'none' does not apply any normalization
##
## pca        - Logical - TRUE creates PCA scores and returns plot
## returnPCA  - logical - TRUE returns PCA scores
## returnMeta - Logical - TRUE returns all input meta data
##
## heatMap    - Logical - TRUE makes a classic heat map
## vioPlot    - Logical - TRUE makes violin plot
## dotPlot    - Logical - TRUE makes a dot plot
## BoxWhisker - Logical - TRUE makes a box and whisker plot
##
## n            - Number of levels used in heat map
## jitter       - Logical - TRUE will offset the x scale
## durbinScale  - integer - Constant for durbin trasform (0 is ln)
##
## zlim     - Numeric vector of length 2 for the c(min, max) plotted value.
##            If zlim is missing or NULL, then all data is plotted
## logScale - Logical - TRUE uses log2 intervals and disable n
##            **Note: This cannot be used with column/p clustering or thresholding
##
## isFuguApp - Logical - True when function call is from the FUGU-app
## shiny_cmp_list - character - List of compound names. Used when isFuguApp = True
## shiny_grp_list - character - List of reference names. Used when isFuguApp = True.
##
## returns the data frame used to plot the data
fuguPlot <-
  function(inDat,
           by = "\\.",
           excludeName = NULL,
           na.rm  = T,
           rmZero = F,
           rCst = F,
           cCst = F,
           scale = c("log",
                     "durbin",
                     "row",
                     "column",
                     "fold",
                     "ratio",
                     "zScore",
                     "none"),
           pca = F,
           returnPCA = F,
           returnMeta = T,
           heatMap = F,
           sLab,
           vioPlot = F,
           dotPlot = F,
           BoxWhisker = F,
           plot3D = F,
           n = 10,
           saveFig = F,
           jitter = 0.2,
           durbinScale = 1,
           zlim,
           logScale = F,
           isFuguApp = F,
           shiny_cmp_list = NULL,
           shiny_grp_list = NULL,
           ...) {
    ## check format of peak list
    inDat <- isInputListValid(inDat)
    
    scale <- if (missing(scale))
      "none"
    else
      match.arg(scale)
    
    ## Take compound list
    com <- inDat$compound
    
    ## Combine any metadata fields to be excluded
    tDat <- metaSep(inDat, excludeName = excludeName)
    metaDat <- tDat[[1]]
    inDat <- tDat[[2]]
    
    ## Store compound list in dataframe label
    sam <- names(inDat)
    rownames(inDat) <- metaDat$compound
    
    ## Get sample labels
    if (missing (sLab)) {
      sLab <- sapply(strsplit(sam, split = by), function (x)
        x[1])
    }
    
    ## Scale data
    if (scale == 'log') {
      inDat <- log(inDat)
    } else if (scale == 'durbin') {
      inDat <- durbin(inDat, c = durbinScale)
    } else if (scale == "row") {
      inDat <-
        sweep(inDat, 1L, rowMeans(inDat, na.rm = TRUE), check.margin = FALSE)
      sx <- apply(inDat, 1L, sd, na.rm = TRUE)
      inDat <- sweep(inDat, 1L, sx, "/", check.margin = FALSE)
    } else if (scale == "column") {
      inDat <-
        sweep(inDat, 2L, colMeans(inDat, na.rm = TRUE), check.margin = FALSE)
      sx <- apply(inDat, 2L, sd, na.rm = TRUE)
      inDat <- sweep(inDat, 2L, sx, "/", check.margin = FALSE)
    } else if (scale == "fold" ||
               scale == "ratio" || scale == 'zScore') {
      if (!isFuguApp) {
        userSel <- select.list(names(inDat), multiple = T,
                               title = "Select reference:")
      }
      else{
        userSel <- shiny_grp_list
      }
      if (!length(userSel) || !nzchar(userSel)) {
        stop('No file selected')
      }
      normCol <- match(userSel, names(inDat))
      inDat <-
        foldChange(inTable = inDat,
                   normVec = normCol,
                   scale = scale)
    }
    
    ## Remove NA
    if (na.rm) {
      printList <- newList <- NULL
      
      if (any(is.na(inDat))) {
        inDat[is.na(inDat)] <- 0
        printList <- c(printList, "NA")
        newList <- c(newList, 0)
      }
      
      if (any(inDat == Inf)) {
        inDat[inDat == Inf] <- range(unlist(inDat), finite = T)[2]
        newList <- c(printList, "Inf")
        newList <- c(newList, max(inDat, na.rm = T))
      }
      
      if (any(inDat == -Inf)) {
        inDat[inDat == -Inf] <- range(unlist(inDat), finite = T)[1]
        printList <- c(printList, "-Inf")
        newList <- c(newList, min(inDat, na.rm = T))
      }
      if (!is.null(printList))
        print(paste(c(
          paste(printList, collapse = ", "),
          "data values were changed to:",
          paste(newList, collapse = ", ")
        ), collapse = " "))
    }
    
    ## Remove all zero data
    if (any(!apply(inDat, 1, function(x)
      all(x == 0))) && rmZero) {
      inDat <- inDat[-which(apply(inDat, 1, function(x)
        all(x == 0))), ]
      
    }
    
    ##################################################
    #################### PLOTTING ####################
    ##################################################
    
    ## Get user selections
    if (BoxWhisker || vioPlot || dotPlot) {
      if (!isFuguApp) {
        userSel <-
          select.list(com, multiple = T, title = "Select compounds:")
      }
      else{
        userSel <- shiny_cmp_list
      }
      if (!length(userSel) || !nzchar(userSel)) {
        print('No compounds selected')
        break
      }
      usrRow <- match(userSel, com)
    }
    
    require(ggplot2)
    require(tidyr)
    require(dplyr)
    
    ## Compare operating systems for visualization
    if (.Platform$OS.type == "windows") {
      quartz <- function(width = 2)
        windows()
    }
    
    if (vioPlot) {
      if (!isFuguApp) {
        userSel <- select.list(unique(sLab), multiple = T)
      }
      else{
        userSel <- shiny_grp_list
      }
      # make sure list isnt empty
      if (!length(userSel) || !nzchar(userSel)) {
        print('No files selected')
        break
      }
      
      usrCol <- myMatch(userSel, sLab)
      temp_dat <- inDat[usrRow, usrCol]
      
      for (numCmp in 1:nrow(temp_dat)) {
        vio_dat <- temp_dat[numCmp, ]
        
        # Need to convert data from wide format to long format for ggplot
        vio_dat <- vio_dat %>%
          gather(key = "Group", value = "Intensity")
        
        # Remove dot index from label
        vio_dat$Group <-
          sapply(strsplit(vio_dat$Group, split = by), function (x)
            x[1])
        
        # Setup window for plotting
        if (!isFuguApp) quartz()
        plt <- ggplot(data = vio_dat,
                      aes(x = Group,
                          y = Intensity,
                          fill = Group)) +
          geom_violin(trim = FALSE,
                      alpha = 0.5, ) +
          geom_jitter(shape = 16,
                      position = position_jitter(jitter),
                      alpha = .3) +
          geom_boxplot(width = .05) +
          labs(title = paste('Compound:', numCmp),
               x = 'Group',
               y = '') +
          theme_classic() +
          theme(title =  element_text(size = 18, face = 'bold'),
                text = element_text(size = 16))
        
        print(plt)
        
        if (saveFig) {
          ggsave(
            paste('Violin_Compound', numCmp, '.pdf', sep = "_"),
            plt,
            height = 8.5,
            width = 11,
            units = "in"
          )
        }
      }
      # Return selected data
      invisible(temp_dat)
    }
    
    if (BoxWhisker) {
      if (!isFuguApp) {
        userSel <- select.list(unique(sLab), multiple = T)
      }
      else{
        userSel <- shiny_grp_list
      }
      # make sure list isnt empty
      if (!length(userSel) || !nzchar(userSel)) {
        print('No files selected')
        break
      }
      
      usrCol <- myMatch(userSel, sLab)
      temp_dat <- inDat[usrRow, usrCol]
      
      for (numCmp in 1:nrow(temp_dat)) {
        bar_dat <- temp_dat[numCmp, ]
        
        # Need to convert data from wide format to long format for ggplot
        bar_dat <- bar_dat %>%
          gather(key = "Group", value = "Intensity")
        
        # Remove dot index from label
        bar_dat$Group <-
          sapply(strsplit(bar_dat$Group, split = by), function (x)
            x[1])
        
        # Setup window for plotting
        if (!isFuguApp) quartz()
        plt <- ggplot(data = bar_dat,
                      aes(x = Group,
                          y = Intensity,
                          fill = Group)) +
          geom_boxplot(width = .4) +
          geom_jitter(shape = 16,
                      position = position_jitter(jitter),
                      alpha = .3) +
          labs(title = paste('Compound:', numCmp),
               x = 'Group',
               y = '') +
          theme_classic() +
          theme(title =  element_text(size = 18, face = 'bold'),
                text = element_text(size = 16))
        
        print(plt)
        
        if (saveFig) {
          ggsave(
            paste('BoxWhisker_Compound', numCmp, '.pdf', sep = "_"),
            plt,
            height = 8.5,
            width = 11,
            units = "in"
          )
        }
      }
      # Return selected data
      invisible(temp_dat)
    }
    
    if (dotPlot) {
      if (!isFuguApp) {
        userSel <- select.list(unique(sLab), multiple = T)
      }
      else{
        userSel <- shiny_grp_list
      }
      # make sure list isnt empty
      if (!length(userSel) || !nzchar(userSel)) {
        print('No files selected')
        break
      }
      
      usrCol <- myMatch(userSel, sLab)
      temp_dat <- inDat[usrRow, usrCol]
      
      for (numCmp in 1:nrow(temp_dat)) {
        dotPlot_dat <- temp_dat[numCmp, ]
        
        # Need to convert data from wide format to long format for ggplot
        dotPlot_dat <- dotPlot_dat %>%
          gather(key = "Group", value = "Intensity")
        
        # Remove dot index from label
        dotPlot_dat$Group <-
          sapply(strsplit(dotPlot_dat$Group, split = by), function (x)
            x[1])
        
        # Setup window for plotting
        if (!isFuguApp) quartz()
        plt <- ggplot(data = dotPlot_dat,
                      aes(x = Group,
                          y = Intensity,
                          fill = Group)) +
          geom_jitter(
            shape = 20,
            position = position_jitter(jitter),
            alpha = 1.0,
            aes(colour = Group)
          ) +
          labs(title = paste('Compound:', numCmp),
               x = 'Group',
               y = '') +
          theme_classic() +
          theme(title =  element_text(size = 18, face = 'bold'),
                text = element_text(size = 16))
        
        print(plt)
        
        if (saveFig) {
          ggsave(
            paste('Dotplot_Compound', numCmp, '.pdf', sep = "_"),
            plt,
            height = 8.5,
            width = 11,
            units = "in"
          )
        }
      }
      # Return selected data
      invisible(temp_dat)
    }
    
    if (heatMap) {
      lDat <- levelFun(inDat,
                       n = n,
                       zlim = zlim,
                       logScale = logScale)
      metImage(
        as.matrix(inDat[nrow(inDat):1, ]),
        col = lDat$col,
        breaks = lDat$breaks,
        rCst = rCst,
        cCst = cCst,
        grid = grid,
        isFuguApp = isFuguApp
      )
      
      invisible(inDat)
    }
    
    if (plot3D) {
      if (!isFuguApp) {
        userSel <- select.list(unique(sLab), multiple = T)
        # make sure list isnt empty
        if (!length(userSel) || !nzchar(userSel)) {
          print('No files selected')
          break
        }
        usrCol <- myMatch(userSel, sLab)
      }
      else{
        usrCol <- shiny_grp_list
      }
      
      temp_dat <- inDat[, usrCol]
      
      if (!("medMz" %in% colnames(metaDat)) ||
          !("medRt" %in% colnames(metaDat))) {
        stop("3D plotting requires medMz and medRt\n")
      }
      
      rt <- metaDat$medRt
      mz <- metaDat$medMz
      
      compound <- row.names(temp_dat)
      
      temp_dat <- cbind(compound, temp_dat)
      
      dat3D <- cbind(rt, mz, temp_dat)
      
      grpsList <- colnames(dat3D)
      grpsList <-  grpsList[-c(1, 2, 3)]
      
      df_temp = data.frame()
      
      for (i in grpsList) {
        df = data.frame(rt = dat3D$rt,
                        mz = dat3D$mz,
                        intensity = dat3D[, i])
        
        df$grpID <- gsub("\\.[0-9]*$", "", i)
        df_temp = rbind(df_temp, df)
      }
      
      require(plotly)
      fig <-
        plot_ly(
          x = df_temp$mz,
          y = df_temp$rt,
          z = df_temp$intensity,
          marker = list(size = 2,  mode = 'scatter3d'),
          color = df_temp$grpID,
          colors = "Dark2"
        )
      print_app(fig)
      return(dat3D)
    }
  }

## This function takes a Maven format list and makes several types of plots
## inDat  - Input data
## pCst   - Logical - TRUE ranks rows by P value
## nCst   - Logical - TRUE sorts samples by name
## scale  - Controls how data is normalized, must be one of the following: 
##          "log", "durbin", "row", "column", or "none".
##					  'log'     applies base 2 log transform
##            'durbin'  applies durbin transform
##            'row'     normalizes by row
##            'column'  normalizes by column
## 					  'zScore'  normalizes relative to a control, sets equal variance
##            'none'    does not apply any normalization
## durbScale    - integer - constant for durbin transform (0 is ln)
## excludeName  - character, names of columns to be excluded from analysis
## na.rm  - Logical - TRUE converts NA to 0, INF to positive max, -INF to negative max
## maxDup - Logical, finds the maximum record in cases of duplicated compounds
## avgRep - Logical, returns the mean intensity for each group
## sdRep  - Logical, returns the standard deviation for each group
## 				  **Note: avgRep and sdRep cannot both be set to T
## inverse - logical, inverts the incoming data (1/data)
## pThresh - logical - excludes data above the pValue threshold
## returnMeta - Logical - TRUE will return all input meta data
##             note: this cannot be used with column/p clustering or thresholding
## rmZero     - logical - TRUE removes rows with all zero intensities
##
## isFuguApp - Logical - True when function call is from the FUGU-app
## selectedReferenceList - character - List of reference names. Used when isFuguApp = True.
##
## returns the data frame used to plot the data
fuguStats <- function(inDat,
                      pCst = F,
                      nCst = F,
                      usrSelect = F,
                      sdRep = F,
                      by = "\\.",
                      excludeName = NULL,
                      scale = c("log",
                                "durbin",
                                "row",
                                "column",
                                "fold",
                                "ratio",
                                "zScore",
                                "none"),
                      inverse = F,
                      durbScale = 1,
                      n = 10,
                      maxDup = F,
                      avgRep = F,
                      na.rm  = T,
                      sLab,
                      dThresh,
                      pCalc = F,
                      pThresh = 1,
                      returnMeta = T,
                      rmZero = F,
                      isFuguApp = F,
                      selectedReferenceList = NULL,
                      ...) {
  ## check format of peaklist
  inDat <- isInputListValid(inDat)
  
  ## Set Scale argument
  scale <- if (missing(scale))
    "none"
  else
    match.arg(scale)
  
  ## Check pThresh
  if (avgRep && pThresh < 1)
    stop("avgRep and pThresh can not be used in combination")
  if (avgRep && sdRep)
    stop("avgRep and sdRep can not be used in combination")
  
  
  ## Take compound list
  com <- inDat$compound
  
  ## Combine any metadata fields to be excluded
  tDat <- metaSep(inDat, excludeName = excludeName)
  metaDat <- tDat[[1]]
  inDat <- tDat[[2]]
  
  ## Store compound list in dataframe label
  sam <- names(inDat)
  rownames(inDat) <- metaDat$compound
  
  ## Get sample labels
  if (missing (sLab)) {
    sLab <- sapply(strsplit(sam, split = by), function (x)
      x[1])
  }
  
  ## Remove duplicate metabolite records and find max for each
  if (maxDup)
    inDat <- findMaxRecord(inDat, excludeName = excludeName)
  
  
  ## Order samples by name
  if (nCst) {
    inDat <- inDat[, order(names(inDat))]
  }
  
  ## Zero any data below the data threshold
  if (!missing(dThresh))
    inDat[inDat < dThresh] <- 0
  
  ## Remove all zero data
  if (any(!apply(inDat, 1, function(x)
    all(x == 0))) && rmZero) {
    inDat <- inDat[-which(apply(inDat, 1, function(x)
      all(x == 0))), ]
  }
  
  ## Get sample labels
  if ((scale == 'zScore' ||
       pCst || pCalc || pThresh < 1) && missing (sLab))
    sLab <- sapply(strsplit(sam, split = by), function (x)
      x[1])
  
  
  ## Average sample replicates
  if (avgRep || sdRep) {
    tLab <- list(factor(sapply(strsplit(names(inDat), split = by), function (x)
      x[1])))
    
    if (avgRep)
      inDat <-
        aggregate(
          t(inDat),
          FUN = function(x)
            mean(x, na.rm = T),
          by = tLab
        )
    else
      inDat <-
        aggregate(
          t(inDat),
          FUN = function(x)
            sd(x, na.rm = T),
          by = tLab
        )
    
    
    uName <- levels(inDat[, 1])
    inDat <- data.frame(t(inDat[, -1]))
    names(inDat) <- uName
    
    ## Reorder colums to match input list
    sam <-
      unique(sapply(strsplit(sam, split = by), function (x)
        x[1]))
    inDat <- inDat[, c(which(names(inDat) == 'compound'),
                       match(sam, names(inDat)))]
  }
  
  
  ## Scale data
  if (scale == 'log') {
    inDat <- log(inDat)
  } else if (scale == 'durbin') {
    inDat <- durbin(inDat, c = durbScale)
  } else if (scale == "row") {
    inDat <-
      sweep(inDat, 1L, rowMeans(inDat, na.rm = TRUE), check.margin = FALSE)
    sx <- apply(inDat, 1L, sd, na.rm = TRUE)
    inDat <- sweep(inDat, 1L, sx, "/", check.margin = FALSE)
  } else if (scale == "column") {
    inDat <-
      sweep(inDat, 2L, colMeans(inDat, na.rm = TRUE), check.margin = FALSE)
    sx <- apply(inDat, 2L, sd, na.rm = TRUE)
    inDat <- sweep(inDat, 2L, sx, "/", check.margin = FALSE)
  } else if (scale == "fold" ||
             scale == "ratio" || scale == 'zScore') {

    userSel <- ifelse(isFuguApp,
                      selectedReferenceList,
                      select.list(names(inDat), multiple = T, title = "Select reference:"))

    if (!length(userSel) || !nzchar(userSel))
      stop('No file selected')
    normCol <- match(userSel, names(inDat))
    inDat <-
      foldChange(inTable = inDat,
                 normVec = normCol,
                 scale = scale)
  }
  
  ## Remove NA and Inf data
  if (na.rm) {
    printList <- newList <- NULL
    
    if (any(is.na(inDat))) {
      inDat[is.na(inDat)] <- 0
      printList <- c(printList, "NA")
      newList <- c(newList, 0)
    }
    
    if (any(inDat == Inf)) {
      inDat[inDat == Inf] <- range(unlist(inDat), finite = T)[2]
      newList <- c(printList, "Inf")
      newList <- c(newList, max(inDat, na.rm = T))
    }
    
    if (any(inDat == -Inf)) {
      inDat[inDat == -Inf] <- range(unlist(inDat), finite = T)[1]
      printList <- c(printList, "-Inf")
      newList <- c(newList, min(inDat, na.rm = T))
    }
    if (!is.null(printList))
      print(paste(c(
        paste(printList, collapse = ", "),
        "data values were changed to:",
        paste(newList, collapse = ", ")
      ), collapse = " "))
  }
  if (inverse)
    inDat <- 1 / inDat
  
  ## Rank values according to ANOVA p value
  if (pCst || pCalc || pThresh < 1) {
    metaDat$pVal <- findPval(sLab = sLab, inDat, type = 'p')
    if (pCst) {
      sidx <- order(metaDat$pVal)
      inDat <- inDat[sidx, ]
      if (ncol(metaDat) > 1)
        metaDat <- metaDat[sidx, ]
      else
        metaDat <-
        data.frame(compound = metaDat[sidx, ], stringsAsFactors = F)
    }
    
    if (pThresh < 1) {
      idxP <- which(metaDat$pVal < pThresh)
      inDat <- inDat[idxP, ]
      if (ncol(metaDat) > 1)
        metaDat <- metaDat[idxP, ]
      else
        metaDat <-
        data.frame(compound = metaDat[idxP, ], stringsAsFactors = F)
    }
  }
  
  ## Allow users to select sub groups of data
  if (usrSelect) {
    userChoice <- mySelect(names(inDat))
    if (length(userChoice) > 0) {
      inDat <- inDat[, userChoice]
      sLab <- sLab[userChoice]
    }
  }
  
  ## Add metadata
  if (returnMeta) {
    inDat <- cbind(metaDat, inDat)
  } else{
    inDat <-
      data.frame(compound = row.names(inDat),
                 inDat,
                 stringsAsFactors = F)
    names(inDat) <- c('compound', sam)
  }
  invisible(inDat)
}

## Function for clustering untargeted mz data
## mzTable  - dataframe - METen format data frame from untargeted csv METen output
## mzThresh - numeric - Mass tolerance in ppm
## rtThresh - numeric - Retention time tolerance in minutes
## pW       - numeric vector, 0 to 1 - weight given to (m/z, rt and cov)
## pThresh  - numeric, 0 to 1 - overall probability threshold
## adList   - numeric vector - list of known adducts
## collapse - logical - T collapses groups to highest intensity signal
## clusterOnly - logical - T returns collapsed groups with multiple signals
fuguCluster <-
  function(mzTable,
           mzThresh = 10,
           rtThresh = .5,
           pW = c(1, 1, 1),
           pThresh = 0.75,
           adList = c(0,
                      1.003355,
                      1.007825,
                      22.9897692809,
                      18.010565,
                      43.99038,
                      34.968853),
           collapse = FALSE,
           clusterOnly = FALSE) {
    mzList <- mzTable$medMz
    rtList <- mzTable$medRt
    pW <- pW / sum(pW)
    
    ## Remove all standard METen meta data
    METNam <-
      c(
        'label',
        'metaGroupId',
        'groupId',
        'goodPeakCount',
        'medMz',
        'medRt',
        'maxQuality',
        'note',
        'compound',
        'compoundId',
        'expectedRtDiff',
        'ppmDiff',
        'parent',
        'adductName',
        'isotopeLabel',
        'formula',
        'pVal'
      )
    mzData <-
      mzTable[, -na.omit(match(METNam, names(mzTable)))] ## remove the non numeric data
    
    ## Find mz diff matrix
    mzMat <- matrix(NA, ncol = length(mzList), nrow = length(mzList))
    for (i in 1:length(mzList)) {
      mzMat[i, ] <- abs(mzList[i] - mzList)
    }
    
    ## Calculate mass prob matrix
    myMin <- function(x) {
      return(min(x, na.rm = TRUE))
    }
    mzProb <- matrix(NA, ncol = ncol(mzMat), nrow = nrow(mzMat))
    for (i in adList) {
      mzTmp <- abs(mzMat - i)
      mzTmp <- (mzTmp / mzList) * 10 ^ 6
      
      ## Keep the best mz match
      mzIdx <- which(mzTmp < mzThresh)
      tmpDat <- cbind(mzTmp[mzIdx], mzProb[mzIdx])
      mzProb[mzIdx] <- apply(tmpDat, 1, myMin)
    }
    
    mzProb <- (abs(mzProb - mzThresh) / mzThresh)
    mzProb[which(is.na(mzProb))] <- 0
    
    ## Find rt diff matix
    rtMat <- matrix(NA, ncol = length(rtList), nrow = length(rtList))
    for (i in 1:length(rtList)) {
      rtMat[i, ] <- abs(rtList[i] - rtList)
    }
    
    ## Create linear probability matix with a threshold
    rtMat[rtMat > rtThresh] <- NA
    rtMat <- abs((rtMat - rtThresh) / rtThresh)
    rtMat[which(is.na(rtMat))] <- 0
    
    ## Find covarance
    corMat <- abs(cor(t(mzData)))
    
    
    ## Make weighted prob matrix
    probMat <- (mzProb * pW[1]) + (rtMat * pW[2]) + (corMat * pW[3])
    
    ## Find compound grouping
    group <- apply(probMat, 1, function(x)
      which(x > pThresh))
    
    oLen <- sapply(group, length)
    nLen <- 1
    while (any(oLen != nLen)) {
      oLen <- sapply(group, length)
      for (i in 1:length(group)) {
        #i == 2201
        group[[i]] <- sort(unique(unlist(group[sapply(group,
                                                      function(x)
                                                        any(match(x, group[[i]]), na.rm = T))])))
      }
      nLen <- sapply(group, length)
    }
    
    group <- unique(group)
    group <-
      unique(group[order(sapply(group, length), decreasing = T)])
    grp <- nGrp <-  rep(NA, length(unlist(group)))
    
    for (i in 1:length(group)) {
      grp[group[[i]]] <- i
      nGrp[group[[i]]]  <- length(group[[i]])
    }
    
    ## Reorder list
    mzTable <- cbind(grp, nGrp, gPar = NA, mzTable)[order(grp), ]
    grp <- grp[order(grp)]
    
    ## Find parent isotope for each group
    uGrp <- unique(grp)
    datMax <-
      as.vector(apply(mzTable[, grep('_', names(mzTable))], 1, sum))
    
    for (i in 1:length(uGrp)) {
      idx <- which(grp == uGrp[i])
      if (any(mzTable[idx, ]$nGrp == 1))
        next
      idx1 <- idx[which.max(datMax[idx])]
      idx2 <- idx[which.min(mzTable$medMz[idx])]
      
      if (idx1 == idx2)
        mzTable[idx1, ]$gPar <- 'Parent'
      else{
        mzTable[idx1, ]$gPar <- 'MaxI'
      }
    }
    
    mzTable <- unique(mzTable)
    
    ## Collapse groups to parent masses
    if (collapse) {
      sTab <- mzTable[mzTable$nGrp > 1,]
      sTab <- sTab[!is.na(sTab$gPar),]
      if (clusterOnly)
        mzTable <- sTab
      else
        mzTable <- rbind(sTab, mzTable[-(which(mzTable$nGrp > 1)),])
      mzTable$nGrp <- 1
    }
    
    return(mzTable)
  }

levelFun <-
  function(x,
           n = 10,
           r = c(0, 1),
           g = c(0, 0),
           b = c(1, 0),
           zlim,
           logScale = T,
           returnInterval = FALSE) {
    x <- suppressWarnings(as.numeric(as.vector(unlist(x))))
    ## Remove NA and Inf data
    if (any(is.na(x)))
      x[is.na(x)] <- 0
    if (any(abs(x) == Inf))
      x[abs(x) == Inf] <- 0
    
    ## Threshold data
    if (!missing(zlim) && length(zlim) == 2 && all(is.numeric(zlim))) {
      zlim <- sort(abs(zlim))
    } else
      zlim <- c(min(abs(x)),  max(abs(x)))
    
    ## Set breaks
    if (logScale) {
      scale <- 2 ^ (0:50)
      xSeq <-
        scale[which(scale >= zlim[1])[1]:which(scale >= zlim[2])[1]]
      n <- length(xSeq)
      xSeq <- unique(c(-rev(xSeq), xSeq))
      
    } else{
      xSeq <- (c(-rev(seq(
        zlim[1], zlim[2], length.out = n
      )),
      seq(zlim[1], zlim[2], length.out = n)))
    }
    
    
    ## Set colors
    col <- c(rgb(rep(r[1], (n - 1)),
                 rep(g[1], (n - 1)),
                 rep(b[1], (n - 1)), alpha = rev((1:n)[-1] / n)),
             'transparent',
             rgb(rep(r[2], (n - 1)),
                 rep(g[2], (n - 1)),
                 rep(b[2], (n - 1)),  alpha = (1:n)[-1] / n))
    
    if (!returnInterval)
      return(list(breaks = xSeq, col = col))
    return(list(
      breaks = xSeq,
      col = col,
      interval = findInterval(x, xSeq)
    ))
  }

metImage <-
  function(x,
           col,
           breaks,
           cCst = FALSE,
           rCst = FALSE,
           grid = TRUE,
           isFuguApp = FALSE) {
    ## Generate levels and colors
    if (missing(col) || missing(breaks)) {
      lDat <- levelFun(x)
      col <- lDat$col
      breaks <- lDat$breaks
    }
    
    ## Set col/row dendrogram
    Colv <- Rowv <- NA
    if (cCst)
      Colv <- NULL
    if (rCst)
      Rowv <- NULL

    if (!isFuguApp) quartz(width = 2)
    w <- (3 + par('mar')[2L]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    par(mar = c(5.1, 1, 4.1, 4.1), las = 1)
    plot.new()
    plot.window(
      xlim = c(0, 1),
      ylim = range(breaks, na.rm = T),
      xaxs = "i",
      yaxs = "i"
    )
    rect(0, breaks[-length(breaks)], 1, breaks[-1L], col = col)
    axis(4, at = breaks)

    #######################################
    if (!isFuguApp) quartz()
    heatmap(
      as.matrix(x),
      Colv = Colv,
      Rowv = Rowv,
      na.rm = T,
      col = col,
      breaks = breaks,
      grid = grid
    )
  }

## Internal function for compressing multiple comound records
findMaxRecord <- function(inDat, excludeName = NULL) {
  ## Get compound name
  com <- inDat$compound
  uCpd <- unique(com)
  
  if (length(com) != nrow(inDat))
    stop("compound column required")
  
  dat <- metaSep(inDat, excludeName = excludeName)
  
  inDat <- dat[[2]]
  inDat <- cbind(compound = com, inDat, stringsAsFactors = F)
  
  ## Remove duplicate metabolite records and find max for each
  for (i in 1:length(uCpd)) {
    idx <- which(com == uCpd[i])
    
    if (nrow(inDat[idx, ]) < 2)
      next
    
    subRow <- apply(inDat[idx, -1], 2, max)
    
    for (j in 1:nrow(inDat[idx, ]))
      inDat[idx, -1][j, ] <- subRow
  }
  
  inDat <- unique(inDat)
  return(inDat)
}

findPval <- function(sLab, mzData, type = 'p') {
  nr <- nrow(mzData)
  pVal <- list()
  sLab <- factor(sLab)
  
  for (i in 1:nr) {
    tVal <-
      tryCatch(
        pVal[[i]] <- aov(as.vector(unlist(mzData[i,])) ~ sLab),
        error = function(er)
          return(NA)
      )
    
    if (is.na(tVal[1])) {
      pVal[[i]] <- NA
      next
    }
    
    if (type == 'p')
      pVal[[i]] <- unclass(summary.aov(tVal))[[1]]$Pr[1]
    
    if (type == 'f')
      pVal[[i]] <- unclass(summary.aov(tVal))[[1]]$"F value"[1]
    
    if (type == 'df')
      pVal[[i]] <- unclass(summary.aov(tVal))[[1]]$"rank"
    
  }
  return(as.vector(unlist(pVal)))
}
