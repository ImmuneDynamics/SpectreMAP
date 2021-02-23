###################################################################################
### SpectreMAP 2 - Spatial plots
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
        dir.create("Output - SpectreMAP 2 - spatial plots")
        setwd("Output - SpectreMAP 2 - spatial plots")
        OutputDirectory <- getwd()

###################################################################################
### Read in spatial data object
###################################################################################

    ### Read in spatial data object

        setwd(InputDirectory)

        spatial.dat <- qread('spatial.dat.qs')

        str(spatial.dat, 3)

###################################################################################
### Big loop to generate spatial plots for all channels and all ROIs
###################################################################################

    ### Specify the name of the mask to use for outlines

        as.matrix(names(spatial.dat[[1]]$RASTERS))
        to.plot <- names(spatial.dat[[1]]$RASTERS)[c(13:25)]
        to.plot

        as.matrix(names(spatial.dat[[1]]$RASTERS))
        main.plot <- names(spatial.dat[[1]]$RASTERS)[c(22)]
        main.plot

        as.matrix(names(spatial.dat[[1]]$MASKS))
        mask <- "cell.mask"
        mask

        as.matrix(names(spatial.dat[[1]]$DATA))
        dat <- "CellData"
        dat

        as.matrix(names(spatial.dat[[1]]$DATA[[dat]]))
        factor.plots <- names(spatial.dat[[1]]$DATA[[dat]])[c(30:31)]
        factor.plots

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
                              cell.dat = dat,
                              cell.col = a)

          }
        }

    ### Spatial plots loop (WITH data points)

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
                              cell.dat = dat,
                              cell.col = a,
                              cell.col.type = 'factor')

          }
        }

