###################################################################################
### Convert TIFF files to FCS and CSV files
###################################################################################

    ### Thomas Ashhurst
    ### https://github.com/ImmuneDynamics/SpectreMAP

###################################################################################
### USER INPUT
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
        setwd("ROIs")
        InputDirectory <- getwd()
        InputDirectory

    ### Create a list of ROIs
        
        rois <- list.dirs(full.names = FALSE, recursive = FALSE)
        as.matrix(rois)
        
    ### Define mask extension
        
        setwd(InputDirectory)
        setwd(rois[[1]])
        as.matrix(list.files())
        
        mask.pattern <- '_MASK'
        mask.ext <- '.tiff'
        
        as.matrix(list.files()[grepl(paste0(mask.pattern, mask.ext), list.files())])
        
    ### Create output directory

        setwd(PrimaryDirectory)
        dir.create("Output - TIFF to FCS")
        setwd("Output - TIFF to FCS")
        OutputDirectory <- getwd()
        
    ### Value modifier
        
        value.modifier <- 65535
        correct.extent <- TRUE
        flip.y <- TRUE
        cofactor <- 1
        
###################################################################################
### END USER INPUT -- RUN ALL OF THE CODE BELOW
###################################################################################


    for(i in rois){
        # i <- rois[[1]]
        setwd(InputDirectory)
        
        ### Initialise the spatial data object with channel TIFF files
        
            spatial.dat <- read.spatial.files(rois = i, roi.loc = getwd(), value.modifier = value.modifier)
            str(spatial.dat, 3)
            as.matrix(names(spatial.dat[[i]]$RASTERS))
            
        ### Add mask file
            
            mask.nme <- list.files()[grepl(paste0(mask.pattern, mask.ext), list.files())]
            length(mask.nme)
            
            if(length(mask.nme) == 1){
                mask.img <- readTIFF(mask.nme)
                mask.img <- raster(mask.img)
                
                if(correct.extent == TRUE){
                    extent(mask.img) <- c(0, dim(mask.img)[2], 0,dim(mask.img)[1]) # Y axis - X axis
                }
                
                if(flip.y == TRUE){
                    mask.img <- flip(mask.img, 'y')
                }
                
                values(mask.img) <- values(mask.img)*value.modifier
                names(mask.img) <- 'cell.mask'
                
                spatial.dat[[i]]$MASKS[['cell.mask']]$maskraster <- mask.img
                
                str(spatial.dat, 3)
                
            } else {
                
                stop('Error - more than 1 mask found')
                
            }
         
        ### Generate polygons and outlines

            spatial.dat <- do.create.outlines(spatial.dat, 'cell.mask')
        
        ### Extract cellular data
        
            spatial.dat <- do.extract(spatial.dat, 'cell.mask', 'CellData')
            
            str(spatial.dat, 3)
            spatial.dat[[1]]$DATA
        
        ### Prep cellular data
        
            cell.dat <- do.pull.data(spatial.dat = spatial.dat,
                                     target.dat = "CellData")
            
            cell.dat
            
            temp <- names(cell.dat)[grepl(mask.pattern, names(cell.dat))]
            temp
            
            cell.dat[[temp]] <- NULL
            cell.dat
            
            as.matrix(names(cell.dat))
            cellular.cols <- names(cell.dat)[c(6:ncol(cell.dat))]
            as.matrix(cellular.cols)
        
        ### Arcsinh transformation
        
            cell.dat <- do.asinh(cell.dat, cellular.cols, cofactor = cofactor)
            
            all.neg <- function(test) -1*abs(test)
            
            y_invert <- cell.dat[['y']]
            y_invert <- all.neg(y_invert)
            cell.dat[['y_invert']] <- y_invert
            
            cell.dat
        
        ### Save CSV and FCS files
        
            setwd(OutputDirectory)
            write.files(cell.dat, i, write.csv = TRUE, write.fcs = TRUE)
            
        ### Cleanup
            
            rm(spatial.dat)
            rm(cell.dat)
            rm(mask.nme)
            rm(mask.img)
            rm(temp)
    }
    