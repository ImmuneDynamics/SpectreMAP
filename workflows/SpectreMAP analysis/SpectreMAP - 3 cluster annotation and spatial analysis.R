###################################################################################
### SpectreMAP 3 - cluster annotation and spatial analysis
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
        dir.create("Output - SpectreMAP 3 - cluster annotation and spatial analysis")
        setwd("Output - SpectreMAP 3 - cluster annotation and spatial analysis")
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
### Annotate clusters
###################################################################################

    setwd(OutputDirectory)
        
    ### Create list of annotations
        
        
    ### Add annotations to cell.dat


###################################################################################
### Calculate area of each 'region'
###################################################################################

    ### Calculate the area of each region in each ROI
        
        area.table <- do.calculate.area(spatial.dat)
        area.table
        
    ### Add region metadata
        
        # cell.dat <- do.add.cols(cell.dat, base.col = 'region', add.dat = region.meta, add.by = 'RegionNum')
        # cell.dat
        
        ## If needed, add 'total' region column
        
        cell.dat$Region <- 'Total'
        cell.dat
        
###################################################################################
### Spatial analysis
###################################################################################

    setwd(OutputDirectory)

    ### Define columns
        
        as.matrix(names(cell.dat))
        
        roi.col <- 'ROI'
        # sample.col <- 'Sample'
        group.col <- 'Group'
        
        pop.col <- 'CellType'
        region.col <- 'Region'
        
    ### Some checks
        
        as.matrix(unique(cell.dat[[roi.col]]))
        as.matrix(unique(cell.dat[[group.col]]))
        as.matrix(unique(cell.dat[[region.col]]))
        area.table
        
        as.matrix(unique(cell.dat[[pop.col]]))
    
    ### Run area analysis on the DATA.TABLE

        reg.dat <- run.spatial.analysis(dat = cell.dat, 
                                        sample.col = roi.col, 
                                        pop.col = pop.col, 
                                        annot.cols = group.col, 
                                        region.col = region.col, 
                                        area.table = area.table) ## Also calculate on total by default
        reg.dat

        as.matrix(names(reg.dat))
        to.plot <- names(reg.dat)[c(3:ncol(reg.dat))]
        to.plot

    ### Pheatmap

        reg.dat.z <- do.zscore(reg.dat, use.cols = to.plot)
        reg.dat.z
        
        reg.dat.z <- reg.dat.z[,colSums(is.na(reg.dat.z))<nrow(reg.dat.z), with = FALSE]
        reg.dat.z
        
        as.matrix(names(reg.dat.z))
        to.plot <- names(reg.dat.z)[c(99:188)]
        
        # names(reg.dat.z) <- gsub("Cells per region", "", names(reg.dat.z))
        # names(reg.dat.z) <- gsub("Cells per 100 um^2 of region", "", names(reg.dat.z))
        # names(reg.dat.z) <- gsub("Percent of cell type in sample", "", names(reg.dat.z))
        # names(reg.dat.z) <- gsub("Percent of cells in region", "", names(reg.dat.z))

        as.matrix(names(reg.dat.z))

        make.pheatmap(reg.dat.z, 
                      sample.col = 'ROI', 
                      plot.cols = to.plot, 
                      is.fold = TRUE, 
                      dendrograms = 'column',
                      row.sep = 2,
                      #cutree_rows = 4, 
                      cutree_cols = 3)

    ### AutoGraphs

        setwd(OutputDirectory)
        dir.create("Autographs")
        setwd("Autographs")
        
        # meas.type <- unique(sub(" -- .*", "", names(sum.dat[,..plot.cols])))
        # meas.type
        
        unique(cell.dat$Group)
        as.matrix(names(reg.dat))
        
        for(i in names(reg.dat)[c(3:ncol(reg.dat))]){
            
            # pop <- sub(".* -- ", "", i) # population
            # meas <- sub(" -- .*", "", i) # measurement
            
            make.autograph(reg.dat,
                           x.axis = 'Group',
                           y.axis = i,
                           y.axis.label = i,
                           
                           grp.order = c('A', 'B'),
                           my_comparisons = list(c('A', 'B')),
                           
                           Variance_test = 'kruskal.test',
                           Pairwise_test = 'wilcox.test',
                           
                           title = i#,
                           #subtitle = meas
            )
        }
        

