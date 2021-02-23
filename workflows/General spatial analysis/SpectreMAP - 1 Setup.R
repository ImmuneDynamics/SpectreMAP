###################################################################################
### SpectreMAP 1 - Setup spatial data object
###################################################################################

    ### Install packages if not previously installed

        # ## Remove previous installations
        # remove.packages('Spectre')
        # remove.packages('SpectreMAP')
        #
        # ## Install (if not already installed)
        # if(!require('devtools')) {install.packages('devtools')}
        # library('devtools')
        #
        # ## Install Spectre
        # install_github("sydneycytometry/spectre")
        # install_github("tomashhurst/SpectreMAP")
        #
        # ## Install any uninstalled dependencies
        # Spectre::package.install(type = 'spatial')
        # Spectre::package.check(type = 'spatial')
        # Spectre::package.load(type = 'spatial')
        #
        # ## At this point, you should restart RStudio

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

    ### Define cell mask extension for different mask types

        setwd(InputDirectory)
        setwd("Masks")

        all.masks <- list.files()
        as.matrix(all.masks)

        cell.mask.ext <- '_ilastik_s2_Multicut Segmentation.tif'
        cell.type.ext <- '_ilastik_s2_Object Predictions.tif'
        # region.ext <- '_ilastik_s2_Simple Segmentation_regions.tif'

        cell.masks <- list.files(pattern = cell.mask.ext)
        cell.masks

        cell.type.masks <- list.files(pattern = cell.type.ext)
        cell.type.masks

        # region.masks <- list.files(pattern = region.ext)
        # region.masks

    ### Import CELL masks

        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = cell.masks,
                                    mask.label = 'cell.mask',
                                    mask.ext = cell.mask.ext)

        str(spatial.dat, 3)

    ### Import CELL TYPE masks

        spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
                                    mask.loc = getwd(),
                                    masks = cell.type.masks,
                                    mask.label = 'cell.type',
                                    mask.ext = cell.type.ext)

        str(spatial.dat, 3)

    ### Import REGION masks

        # spatial.dat <- do.add.masks(spatial.dat = spatial.dat,
        #                             mask.loc = getwd(),
        #                             masks = region.masks,
        #                             mask.label = 'region',
        #                             mask.ext = region.ext)
        #
        # str(spatial.dat, 3)

###################################################################################
### Generate polygons and outlines
###################################################################################

    ### Generate polygons and outlines

        library(dplyr) ## Need to add to the package functions, as well as the do.create.outlines function

        spatial.dat <- do.create.outlines(spatial.dat, 'cell.mask')
        spatial.dat <- do.create.outlines(spatial.dat, 'cell.type')
        # spatial.dat <- do.create.outlines(spatial.dat, 'region')

    ### Checks

        str(spatial.dat, 3)

###################################################################################
### Some quick QC spatial plots
###################################################################################

    ### Setup to do some quick spatial plots

        setwd(OutputDirectory)
        dir.create("QC spatial plots")
        setwd("QC spatial plots")

        ## Choose an ROI to use
        as.matrix(names(spatial.dat))

        plot.roi <- names(spatial.dat)[c(1)]
        plot.roi

        ## Choose some markers to plot
        as.matrix(names(spatial.dat[[1]]$RASTERS))

        to.plot <- names(spatial.dat[[1]]$RASTERS)[c(13:25)]
        to.plot

    ### Make plots

        for(i in to.plot){
          make.spatial.plot(spatial.dat,
                            image.roi = plot.roi,
                            image.channel = i,
                            mask.outlines = 'cell.mask')
          rm(i)
        }

###################################################################################
### Extract cellulra data from images using masks
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
