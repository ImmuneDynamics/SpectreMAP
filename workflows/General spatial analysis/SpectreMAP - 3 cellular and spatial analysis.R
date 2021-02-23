###################################################################################
### SpectreMAP 3 - cellular and spatial analysis
###################################################################################

    ### Load libraries

        library('Spectre')
        library('SpectreMAP')

        Spectre::package.check(type = 'spatial')
        Spectre::package.load(type = 'spatial')

    ### Set PrimaryDirectory

        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Set InputDirectory

        setwd(PrimaryDirectory)
        setwd("Output - SpectreMAP 1 - setup/")
        InputDirectory <- getwd()

    ### Set metadata directory

        setwd(PrimaryDirectory)
        setwd("metadata")
        MetaDirectory <- getwd()

    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("Output - SpectreMAP 3 - cellular and spatial analysis")
        setwd("Output - SpectreMAP 3 - cellular and spatial analysis")
        OutputDirectory <- getwd()

###################################################################################
### Read in spatial data object
###################################################################################

    ### Read in spatial data object

        setwd(InputDirectory)

        spatial.dat <- qread('spatial.dat.qs')

        str(spatial.dat, 3)

###################################################################################
### Calculate area of each 'region'
###################################################################################

    ### Calculate the area of each region in each ROI

        for(i in names(spatial.dat)){
          spatial.dat[[i]]$DATA$CellData$Region <- "Total"
        }

        area.table <- do.calculate.area(spatial.dat, region = 'Region')
        area.table

        area.table <- data.frame('ROI' = names(spatial.dat), 'Total' = rep(500, length(names(spatial.dat))))
        area.table <- as.data.table(area.table)
        area.table

###################################################################################
### Pull cellular data out of spatial data object and add annotations
###################################################################################

    ### Pull

        cell.dat <- do.pull.data(spatial.dat, 'CellData')
        cell.dat

        as.matrix(names(cell.dat))

        cellular.cols <- names(cell.dat)[c(18:30)]
        as.matrix(cellular.cols)

        cell.dat <- do.asinh(cell.dat, cellular.cols, cofactor = 1)
        cell.dat

        cellular.cols <- paste0(cellular.cols, "_asinh")
        as.matrix(cellular.cols)

    ### Read in metadata files

        # setwd(MetaDirectory)
        #
        # sample.meta <- fread("ROIs and samples.csv")
        # sample.meta
        #
        # cell.dat <- do.add.cols(cell.dat, base.col = 'ROI', add.dat = sample.meta, add.by = 'ROI')
        # cell.dat

    ### Add annotations

        # cell.type.meta <- fread("cell types.csv")
        # cell.type.meta

        # cell.dat <- do.add.cols(cell.dat, base.col = 'cell.type', add.dat = cell.type.meta, add.by = 'CellTypeNum')
        # cell.dat

    ### Add annotations

        # region.meta <- fread("regions.csv")
        # region.meta
        #
        # cell.dat <- do.add.cols(cell.dat, base.col = 'region', add.dat = region.meta, add.by = 'RegionNum')
        # cell.dat

    ## Also adjust area table

        # region.meta
        # area.table
        #
        # names(area.table) <- c("ROI", "White pulp", "Red pulp", "Other")
        # area.table

###################################################################################
### x
###################################################################################

    setwd(OutputDirectory)

    ### Expression heatmap

        exp <- do.aggregate(cell.dat, cellular.cols, 'cell.type', 'mean')
        exp

        make.pheatmap(exp, 'cell.type', cellular.cols)

    ### Run area analysis on the DATA.TABLE

        ### test
        # cell.dat$Region <- "Total"
        # area.table <- data.frame(ROI = '20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac', 'Total' = 250500)
        # area.table <- as.data.table(area.table)

        reg.dat <- SpectreMAP::do.region.analysis(dat = cell.dat, 'ROI', 'cell.type', 'Region', area.table = area.table)
        reg.dat

        as.matrix(names(reg.dat))
        to.plot <- names(reg.dat)[c(8:13)]
        to.plot

    ### Pheatmap

        reg.dat.z <- do.zscore(reg.dat, use.cols = to.plot)
        reg.dat.z

        # names(reg.dat.z) <- gsub("Cells per region", "", names(reg.dat.z))
        # names(reg.dat.z) <- gsub("Cells per 100 um^2 of region", "", names(reg.dat.z))
        # names(reg.dat.z) <- gsub("Percent of cell type in sample", "", names(reg.dat.z))
        # names(reg.dat.z) <- gsub("Percent of cells in region", "", names(reg.dat.z))

        as.matrix(names(reg.dat.z))

        make.pheatmap(reg.dat.z, sample.col = 'ROI', plot.cols = to.plot, is.fold = TRUE, cutree_rows = 4, cutree_cols = 3)

    ### AutoGraphs










