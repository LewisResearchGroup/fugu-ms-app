################################################################################
##                                                                            ##
##               Command-line script to process BSI samples                   ##
##                                                                            ##
################################################################################

# Set Working Directory
setwd("~/Desktop/Fugo")

#Load in source code library, Fugu
source(file.choose())

#Load S1 - Data
dat <- fuguRead()


#Heatmap
fuguPlot(dat, heatMap = T, rCst = T)
#3D plot
plot3D_data <- fuguPlot(dat, plot3D = T)
#ViolinPlot
vio <- fuguPlot(dat, vioPlot = T)
#BoxPlot
bar <- fuguPlot(dat, BoxWhisker = T)
#DotPlot
dot <- fuguPlot(dat, dotPlot = T)


############################
by = "\\."
sam <- names(dat)
sLab <- sapply(strsplit(sam, split = by), function (x)
  x[1])
sLab <- sLab[-c(1, 2, 3)]
grpSel <- select.list(unique(sLab), multiple = T)


com <- dat$compound
cmpSel <-
  select.list(com, multiple = T, title = "Select compounds:")

#3D plot
plot3D_data <- fuguPlot(dat,
                        plot3D = T,
                        isFuguApp = T,
                        shiny_grp_list = grpSel)
#ViolinPlot
vio <- fuguPlot(
  dat,
  vioPlot = T,
  isFuguApp = T,
  shiny_grp_list = grpSel,
  shiny_cmp_list = cmpSel
)
#BoxPlot
bar <- fuguPlot(
  dat,
  BoxWhisker = T,
  isFuguApp = T,
  shiny_grp_list = grpSel,
  shiny_cmp_list = cmpSel
)
#DotPlot
dot <- fuguPlot(
  dat,
  dotPlot = T,
  isFuguApp = T,
  shiny_grp_list = grpSel,
  shiny_cmp_list = cmpSel
)

# Run ANOVA on all markers and filter insignificant
alpha_val <- 0.05
pDat <- fuguStats(dat, pCalc = T)
num_samples = nrow(pDat)
idx_p <-
  pDat$pVal < alpha_val / (num_samples)      # bonferroni Correction and filter
pDat <- pDat[idx_p, ]

alpha_val / (num_samples)
# Remove NA data from frame
pDat <- na.omit(pDat)

#
fuguPlot(pDat, heatMap = T, rCst = T)

# Calculate the average threshold per species
mDat <- fuguStats(pDat, avgRep = T)
fuguPlot(mDat, heatMap = T, rCst = T)


# Remove markers if the max average across all groups < 20000
thresh_num <- 20000
idx_num <- apply(abs(metaSep(mDat)$data), 1, max) > thresh_num
length(idx_num[idx_num == TRUE])


# filter markers from data set if the threshold is less than
# 4-fold change in comparison to MHB
thresh_fold <- 4
# select MHB or 1 as reference when calling this function
fDat <- fuguStats(mDat, scale = "fold")
idx_fold <- (apply(abs(metaSep(fDat)$data), 1, max) > thresh_fold)

## Find data with correct intensity, pValue and fold change (533 markers)
idx <- which((idx_num + idx_fold) == 2)
out <- pDat[idx, ]

fuguPlot(out, heatMap = T, rCst = T)

fuguWrite(out)

################################################################################
##                                                                            ##
##               Command-line script for Clustering                           ##
##                                                                            ##
################################################################################
# resolution of masses to delta 5ppm between mass range of 100-500 m/z
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
  METCluster(
    out,
    pThresh = 0.8,
    collapse = FALSE,
    clusterOnly = FALSE,
    adList = mz_ad,
    pW = c(1, 1, 1)
  )
cDatCleaned <- cDat[!duplicated(cDat$grp), ]
METWrite(cDatCleaned)				#2
METPlot(
  cDatCleaned,
  heat = T,
  scale = "row",
  rCst = T,
  grid = F,
  cCst = F
)
################################################################################
##                                                                            ##
##               Command-line script to Visualize Data                        ##
##                                                                            ##
################################################################################

vio <- fuguPlot(out, vioPlot = T, saveFig = T)
bar <- fuguPlot(out, barPlot = T, saveFig = T)
dot <- fuguPlot(out, dotPlot = T, saveFig = T)

################################################################################
