#' do.calculate.area
#'
#' @import data.table
#'
#' @export

do.calculate.area <- function(dat,
                              region = NULL){
  
  ### Setup
  
  # dat <- spatial.dat
  # region <- 'regions'
  
  ### Preparation
  
  roi.names <- names(dat)
  
  ### Processing
  
  area.list <- list()
  
  for(i in roi.names){
    # i <- roi.names[[1]]

    ## If region masks are present
        if(!is.null(region)){
          poly.names <- dat[[i]]$MASKS[[region]]$polygons@data
          areas <- area(dat[[i]]$MASKS[[region]]$polygons)
          
          areas <- as.data.table(t(areas))
          names(areas) <- as.character(poly.names[[1]])
          
          total.area <- sqrt(dim(dat[[i]]$RASTERS[[1]])[1] * dim(dat[[i]]$RASTERS[[1]])[2])
          total.area <- as.data.table(total.area)
          names(total.area) <- 'Total'
          areas <- rbind(areas, total.area)
          area.list[[i]] <- areas
        }
    
    ## If region masks are NOT present
        if(is.null(region)){
          total.area <- sqrt(dim(dat[[i]]$RASTERS[[1]])[1] * dim(dat[[i]]$RASTERS[[1]])[2])
          total.area <- as.data.table(total.area)
          
          names(total.area) <- 'Total'
          area.list[[i]] <- total.area
        }
  }
  
  area.list
  
  area.res <- rbindlist(area.list, fill = TRUE)
  area.res <- cbind(roi.names, area.res)
  names(area.res)[1] <- "ROI"
  
  ### Return
  
  return(area.res)
  
}