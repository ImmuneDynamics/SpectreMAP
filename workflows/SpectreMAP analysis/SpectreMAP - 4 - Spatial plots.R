###################################################################################
### SpectreMAP 4 - spatial plots
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
        setwd("Output - SpectreMAP 2 - clustering/")
        InputDirectory <- getwd()
        
    ### Set metadata directory
        
        setwd(PrimaryDirectory)
        setwd("metadata")
        MetaDirectory <- getwd()
        
    ### Create output directory
        
        setwd(PrimaryDirectory)
        dir.create("Output - SpectreMAP 4 - spatial plots")
        setwd("Output - SpectreMAP 4 - spatial plots")
        OutputDirectory <- getwd()

###################################################################################
### Read in data
###################################################################################
        
    ### Read in spatial data object
        
        setwd(InputDirectory)
        setwd("../Output - SpectreMAP 1 - setup/")
        
        list.files(getwd(), '.qs')
        
        spatial.dat <- qread('spatial.dat.qs')
        
        str(spatial.dat, 3)
        
    ### Read in data.table
        
        setwd(InputDirectory)
        list.files(getwd(), '.csv')
        
        cell.dat <- fread("cell.dat.csv")
        cell.dat

###################################################################################
### Setup
###################################################################################

    ### Specify masks and data
        
        as.matrix(names(spatial.dat[[1]]$MASKS))
        mask <- "cell.mask"
        mask
        
    ### Specify rasters for plotting
        
        as.matrix(names(spatial.dat[[1]]$RASTERS))
        to.plot <- names(spatial.dat[[1]]$RASTERS)[c(13:25)]
        to.plot

    ### Cellular data
        
        dat <- cell.dat
        
        ext <- "_asinh_rescaled"
        
        # as.matrix(names(spatial.dat[[1]]$DATA))
        # dat <- "CellData"
        # dat
        
    ### Specify plotting for factors
        
        as.matrix(names(spatial.dat[[1]]$RASTERS))
        main.plot <- names(spatial.dat[[1]]$RASTERS)[c(15)]
        main.plot
        
        as.matrix(names(spatial.dat[[1]]$DATA[[dat]]))
        
        as.matrix(names(dat))
        
        factor.plots <- names(dat)[c(63,34,35)]
        factor.plots

###################################################################################
### Plotting loop
###################################################################################
        
    ### Spatial plots loop (FACTOR data points)
        
        for(i in names(spatial.dat)){
            # i <- names(spatial.dat)[1]
            
            setwd(OutputDirectory)
            dir.create("Spatial plots with data points - factors")
            setwd("Spatial plots with data points - factors")
            dir.create(i)
            setwd(i)
            
            for(a in factor.plots){
                
                make.spatial.plot(spatial.dat,
                                  image.roi = i,
                                  image.channel = main.plot,
                                  mask.outlines = mask,
                                  cell.dat = dat[dat[['ROI']] == i,],
                                  cell.col = a,
                                  cell.col.type = 'factor')
            }
        }
        
    ### Spatial plots loop (WITH data points)
        
        for(i in names(spatial.dat)){
            # i <- names(spatial.dat)[1]
            
            setwd(OutputDirectory)
            dir.create("Spatial plots with data points")
            setwd("Spatial plots with data points")
            dir.create(i)
            setwd(i)
            
            for(a in to.plot){
                
                make.spatial.plot(spatial.dat,
                                  image.roi = i,
                                  image.channel = a,
                                  mask.outlines = mask,
                                  cell.dat = dat[dat[['ROI']] == i,],
                                  cell.col = paste0(a, ext))
            }
        }
        
    ### Spatial plots loop (no data points)
        
        for(i in names(spatial.dat)){
            # i <- names(spatial.dat)[1]
            
            setwd(OutputDirectory)
            dir.create("Spatial plots")
            setwd("Spatial plots")
            dir.create(i)
            setwd(i)
            
            for(a in to.plot){
                
                make.spatial.plot(spatial.dat,
                                  image.roi = i,
                                  image.channel = a,
                                  mask.outlines = mask)
            }
        }
        
