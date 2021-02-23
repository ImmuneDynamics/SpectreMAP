###################################################################################
### Load packages and set directories
###################################################################################

    ### Load packages
        library(Spectre)
        library(SpectreMAP)

        package.check('spatial')
        package.load('spatial')

    ### Set PrimaryDirectory
        dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
        setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
        getwd()
        PrimaryDirectory <- getwd()
        PrimaryDirectory

    ### Create output directory
        setwd(PrimaryDirectory)
        dir.create("Output - SpectreMAP demo")
        setwd("Output - SpectreMAP demo")
        OutputDirectory <- getwd()

###################################################################################
### Extract 'cellular' data using masks
###################################################################################

    ### Load demo data
        
        spatial.dat <- SpectreMAP::demo.spatial
        # spatial.dat <- demo.spatial
        
        str(spatial.dat, 3)
        
    ### Remove calculated components for the purpose of the demo
        
        for(i in names(spatial.dat)){
            spatial.dat[[i]]$MASKS$cell.mask$polygons <- NULL
            spatial.dat[[i]]$MASKS$cell.mask$outlines <- NULL
            spatial.dat[[i]]$MASKS$cell.mask$centroids <- NULL
            
            spatial.dat[[i]]$MASKS$cell.type$polygons <- NULL
            spatial.dat[[i]]$MASKS$cell.type$outlines <- NULL
            spatial.dat[[i]]$MASKS$cell.type$centroids <- NULL
            
            spatial.dat[[i]]$DATA <- NULL
        }
        
    ### Review demo datq
        
        str(spatial.dat, 3)
        str(spatial.dat[[1]]$MASKS, 3)
        
###################################################################################
### Calculate polygons and outlines from mask rasters
###################################################################################        
        
    ### Create polygons, outlines, and centroids
        
        str(spatial.dat[[1]]$MASKS, 2)
        
        spatial.dat <- do.create.outlines(spatial.dat = spatial.dat, mask.name = "cell.mask")
        
        str(spatial.dat[[1]]$MASKS, 2)
        
        spatial.dat <- do.create.outlines(spatial.dat = spatial.dat, mask.name = "cell.type")
        
        str(spatial.dat[[1]]$MASKS, 2)
        
    ### Make some plots
        
        
        
        
        
###################################################################################
### Calculate 'cellular' data from masks
###################################################################################   
        
    ### Calculate cellular data

        spatial.dat <- do.extract(dat = spatial.dat, mask = "cell.mask", name = "CellData", fun = "mean")

        str(spatial.dat, 3)
        spatial.dat[[1]]$DATA

###################################################################################
### Make some spatial plots
###################################################################################

    ### Review ROI names

        setwd(OutputDirectory)
        dir.create("Spatial plots")
        setwd("Spatial plots")

        as.matrix(names(spatial.dat))

    ### Review channel names

        as.matrix(names(spatial.dat[[1]]$RASTERS))

    ### Plots

        make.spatial.plot(spatial.dat = spatial.dat,
                          image.roi = '20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac',
                          image.channel = "CD20_Dy161")


        make.spatial.plot(spatial.dat = spatial.dat,
                          image.roi = '20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac',
                          image.channel = "CD20_Dy161",
                          mask.outlines = "cell_mask")

        make.spatial.plot(spatial.dat = spatial.dat,
                          image.roi = '20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac',
                          image.channel = "CD20_Dy161",
                          mask.outlines = "cell_mask",
                          cell.dat = "CellData",
                          cell.col = "CD20_Dy161")

###################################################################################
### Pull and merge cell data from each ROI
###################################################################################

    ### Pull cell data from each ROI and merge

        setwd(OutputDirectory)
        dir.create("Cellular plots")
        setwd("Cellular plots")

        cell.dat <- do.pull.data(spatial.dat = spatial.dat,
                                 target.dat = "CellData")

        cell.dat

    ### Define cellular and clustering columns

        as.matrix(names(cell.dat))
        cellular.cols <- names(cell.dat)[c(17:29)]

        as.matrix(names(cell.dat))
        clustering.cols <- names(cell.dat)[c(18:24,29)]

    ### Arcsinh

        cell.dat <- do.asinh(cell.dat, cellular.cols, cofactor = 1)

        cellular.cols <- paste0(cellular.cols, "_asinh")
        clustering.cols <- paste0(clustering.cols, "_asinh")

    ### Clustering and DR

        cell.dat <- run.flowsom(cell.dat, clustering.cols, meta.k = 10)
        cell.dat <- run.umap(cell.dat, clustering.cols)

    ### Make some data plots

        make.colour.plot(cell.dat, "UMAP_X", "UMAP_Y", col.axis = "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
        make.multi.plot(cell.dat, "UMAP_X", "UMAP_Y", cellular.cols)

    ### Make some spatial plots

        make.spatial.plot(spatial.dat = spatial.dat,
                          image.roi = '20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac',
                          image.channel = "CD20_Dy161",
                          mask.outlines = "cell_mask",

                          cell.dat = cell.dat[cell.dat[["ROI"]] == '20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac',],
                          cell.col = "CD20_Dy161")

        make.spatial.plot(spatial.dat = spatial.dat,
                          image.roi = '20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac',
                          image.channel = "CD20_Dy161",
                          mask.outlines = "cell_mask",

                          cell.dat = cell.dat[cell.dat[["ROI"]] == '20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac',],
                          cell.col = "FlowSOM_metacluster",
                          cell.col.type = 'factor')

###################################################################################
### Expression heatmap
###################################################################################

    exp <- do.aggregate(dat = cell.dat, use.cols = cellular.cols, by = 'FlowSOM_metacluster', func = 'mean')

    make.pheatmap(dat = exp, sample.col = 'FlowSOM_metacluster', plot.cols = cellular.cols)


###################################################################################
### Save CSV and FCS files
###################################################################################

    ### Set directory

        setwd(OutputDirectory)
        dir.create("Data")
        setwd("Data")

    ### Create a flipped y-axis

        all.neg <- function(test) -1*abs(test)

        y_invert <- cell.dat[['y']]
        y_invert <- all.neg(y_invert)
        cell.dat[['y_invert']] <- y_invert

        cell.dat

    ### Save CSV and FCS files

        write.files(cell.dat, "Sp_AllData", write.csv = TRUE, write.fcs = TRUE)
        write.files(cell.dat, "Sp", divide.by = "ROI", write.csv = TRUE, write.fcs = TRUE)

    ### Save RDS file

        saveRDS(spatial.dat, 'spatial.data.rds')
        save(spatial.dat, file = 'demo.spatial.RData')




