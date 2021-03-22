###################################################################################
### SpectreMAP 1 - Setup spatial data object
###################################################################################

    ### Load packages

        library(Spectre)
        library(SpectreMAP)

        package.check(type = 'spatial')
        package.load(type = 'spatial')

    ### Set PrimaryDirectory

        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set InputDirectory

        setwd(PrimaryDirectory)
        setwd("data")
        InputDirectory <- getwd()

    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("Output - SpectreMAP 1 - setup")
        setwd("Output - SpectreMAP 1 - setup")
        OutputDirectory <- getwd()

###################################################################################
### Read in TIFF files and masks
###################################################################################

    ### Initialise the spatial data object with channel TIFF files

        setwd(InputDirectory)
        setwd("ROIs")

        rois <- list.dirs(full.names = FALSE, recursive = FALSE)
        as.matrix(rois)

        spatial.dat <- read.spatial.files(rois = rois, roi.loc = getwd())

        str(spatial.dat, 3)

###################################################################################
### Read in masks files
###################################################################################

    ### Define cell mask extension for different mask types

        setwd(InputDirectory)
        setwd("Masks")

        all.masks <- list.files()

    ### Import CELL masks
        
        as.matrix(all.masks)
        cell.mask.ext <- '_ilastik_s2_Object Identities.tif'
        
        cell.masks <- list.files(pattern = cell.mask.ext)
        cell.masks

        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = cell.masks,
                                    mask.label = 'cell.mask',
                                    mask.ext = cell.mask.ext)

        str(spatial.dat, 3)

    ### Import CELL TYPE masks
        
        as.matrix(all.masks)
        cell.type.ext <- '_ilastik_s2_Object Predictions.tif'
        
        cell.type.masks <- list.files(pattern = cell.type.ext)
        cell.type.masks

        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = cell.type.masks,
                                    mask.label = 'cell.type',
                                    mask.ext = cell.type.ext)

        str(spatial.dat, 3)

    ### Import REGION masks
        
        as.matrix(all.masks)
        region.ext <- '_ilastik_s2_Simple Segmentation.tif'

        region.masks <- list.files(pattern = region.ext)
        region.masks

        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = region.masks,
                                    mask.label = 'region',
                                    mask.ext = region.ext)

        str(spatial.dat, 3)

###################################################################################
### Generate polygons and outlines
###################################################################################

    ### Generate polygons and outlines

        library(dplyr) ## Need to add to the package functions, as well as the do.create.outlines function

        spatial.dat <- do.create.outlines(spatial.dat, 'cell.mask')
        spatial.dat <- do.create.outlines(spatial.dat, 'cell.type')
        spatial.dat <- do.create.outlines(spatial.dat, 'region')

    ### Checks

        str(spatial.dat, 3)

###################################################################################
### Extract cellular data
###################################################################################

    ### Extract cellular data
        
        spatial.dat <- do.extract(spatial.dat, 'cell.mask', 'CellData')
        
        str(spatial.dat, 3)
        spatial.dat[[1]]$DATA
        
###################################################################################
### Save spatial.dat object as RDS file
###################################################################################
        
    ### Save as quick serial (qs) file
        
        setwd(OutputDirectory)
        qsave(spatial.dat, "spatial.dat.qs")

###################################################################################
### Some quick QC spatial plots
###################################################################################

    ### Choose an ROI to use
        as.matrix(names(spatial.dat))

        roi.plot <- names(spatial.dat)[c(1)]
        roi.plot

    ### Choose some markers to plot
        as.matrix(names(spatial.dat[[1]]$RASTERS))

        exp.plot <- names(spatial.dat[[1]]$RASTERS)[c(13:25)]
        exp.plot

    ### Choose some cell data factors to plot
        
        as.matrix(names(spatial.dat[[1]]$DATA$CellData))
        
        factor.plot <- names(spatial.dat[[1]]$DATA$CellData)[c(30,31)]
        factor.plot
        
        as.matrix(names(spatial.dat[[1]]$RASTERS))
        
        background <- names(spatial.dat[[1]]$RASTERS)[15]
        background
        
    ### Make plots - no data points

        setwd(OutputDirectory)
        dir.create("QC spatial plots")
        setwd("QC spatial plots")
        
        for(i in exp.plot){
          make.spatial.plot(spatial.dat,
                            image.roi = roi.plot,
                            image.channel = i,
                            mask.outlines = 'cell.mask')
          rm(i)
        }
        
    ### Make plots - with data points showing expression
        
        setwd(OutputDirectory)
        dir.create("QC spatial plots - with data points")
        setwd("QC spatial plots - with data points")
        
        for(i in exp.plot){
            make.spatial.plot(spatial.dat,
                              image.roi = roi.plot,
                              image.channel = i,
                              mask.outlines = 'cell.mask', 
                              cell.dat = 'CellData', 
                              cell.col = i)
            rm(i)
        }
        
    ### Make plots - with data points showing factors
        
        setwd(OutputDirectory)
        dir.create("QC spatial plots - with factor points")
        setwd("QC spatial plots - with factor points")
        
        for(i in factor.plot){
            make.spatial.plot(spatial.dat,
                              image.roi = roi.plot,
                              image.channel = background,
                              mask.outlines = 'cell.mask', 
                              cell.dat = 'CellData', 
                              cell.col = i, 
                              cell.col.type = 'factor')
            rm(i)
        }
        

