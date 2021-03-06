#' make.spatial.plot
#'
#' @import data.table
#'
#' @export

make.spatial.plot <- function(spatial.dat, # spatial data object
                              image.roi, # name of ROI
                              image.channel, # name of channel

                              ## Options for adding cell outlines
                              mask.outlines = NULL, # character -- the outlines in spatial.dat object

                              ## Options for adding cellular data
                              cell.dat = NULL, # can be character (if it's data within spatial.dat) or a data.table
                              cell.col = NULL, # column for colouration

                              ## Other settings (with defaults)
                              image.y.flip = TRUE,
                              image.mask.size = 0.1,
                              image.mask.colour = "gold",
                              image.min.threshold = 0.00,
                              image.max.threshold = 0.99,
                              image.blank = FALSE,

                              cell.x = "x",
                              cell.y = "y",
                              cell.col.type = "numeric",
                              cell.colours = "spectral",
                              cell.col.min.threshold = 0.01,
                              cell.col.max.threshold = 0.995,

                              title = paste0(image.roi),
                              dot.size = 1,
                              dot.alpha = 1,
                              align.xy.by = cell.dat, # choose a data frame to set absolute limits for X/Y/colour
                              align.col.by = cell.dat,
                              save.to.disk = TRUE,
                              path = getwd(),
                              plot.width = 9,
                              plot.height = 7,
                              blank.axis = FALSE)
{

  ### TESTING
  # library(raster)
  # library(data.table)
  # library(tiff)
  # library(ggplot2)
  #
  # spatial.dat = spatial.dat
  #
  # spatial.dat$meta.data
  #
  # roi = "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac"
  # roi.marker = "CD20_Dy161"

  # cell.dat <- spatial.dat$cell.dat.means.filtered
  # cell.dat <- cell.dat[cell.dat[["ImageName"]] == "20171228_spleen315_500x500_editedforFAS_s1_p9_r2_a2_ac_ilastik_s2_Probabilities_mask.tiff",]
  # cell.dat = cell.dat
  # cell.x = "X"
  # cell.y = "Y"
  # cell.colour = 'CD20'
  #
  # add.outlines = TRUE
  # flip.y.axis = TRUE

  message("This is a developmental Spectre-spatial function that is still in testing phase with limited documentation. We recommend only using this function if you know what you are doing.")

  ### Check that necessary packages are installed
  if(!is.element('Spectre', installed.packages()[,1])) stop('Spectre is required but not installed')
  if(!is.element('ggplot2', installed.packages()[,1])) stop('ggplot2 is required but not installed')
  if(!is.element('scales', installed.packages()[,1])) stop('scales is required but not installed')
  if(!is.element('colorRamps', installed.packages()[,1])) stop('colorRamps is required but not installed')
  if(!is.element('ggthemes', installed.packages()[,1])) stop('ggthemes is required but not installed')
  if(!is.element('RColorBrewer', installed.packages()[,1])) stop('RColorBrewer is required but not installed')
  if(!is.element('raster', installed.packages()[,1])) stop('raster is required but not installed')
  if(!is.element('rgeos', installed.packages()[,1])) stop('rgeos is required but not installed')

  ### Require packages
  require(Spectre)
  require(ggplot2)
  require(scales)
  require(colorRamps)
  require(ggthemes)
  require(RColorBrewer)
  require(raster)
  require(rgeos)

  ### Compatability conversions

  roi <- image.roi
  roi.marker <- image.channel

  #cell.dat
  cell.colour <- cell.col

  add.outlines <- image.outlines <- mask.outlines
  flip.y.axis <- image.y.flip

  cell.colour.type <- cell.col.type

  raster.mask.size <- image.mask.size
  raster.mask.colour <- image.mask.colour
  raster.min.threshold <- image.min.threshold
  raster.max.threshold <- image.max.threshold

  col.min.threshold <- cell.col.min.threshold
  col.max.threshold <- cell.col.max.threshold

  colours <- cell.colours

  ### Colour setup

  # Jet
  if(colours == "jet"){
    colour.scheme <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  }

  # Spectral
  if(colours == "spectral"){
    spectral.list <- colorRampPalette(brewer.pal(11,"Spectral"))(50)
    spectral.list <- rev(spectral.list)
    colour.scheme <- colorRampPalette(c(spectral.list))
  }

  # Viridis
  if(colours == "viridis"){
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "viridis")(50)))
  }

  # Inferno
  if(colours == "inferno"){
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "inferno")(50)))
  }

  #Magma
  if(colours == "magma"){
    colour.scheme <- colorRampPalette(c(viridis_pal(option = "magma")(50)))
  }

  ### cell.dat setup

  if(!is.null(cell.dat)){

    if(is.character(cell.dat) == TRUE){
      temp <- spatial.dat[[roi]]$DATA[[cell.dat]]
      cell.dat <- temp
    }


    if(cell.colour.type == "numeric"){

      # Dot point colouration
      if(is.null(align.col.by) == TRUE){
        ColrMin <- quantile(cell.dat[[cell.colour]], probs = c(col.min.threshold))
        ColrMax <- quantile(cell.dat[[cell.colour]], probs = c(col.max.threshold))
      }

      if(is.null(align.col.by) == FALSE){
        ColrMin <- quantile(align.col.by[[cell.colour]], probs = c(col.min.threshold))
        ColrMax <- quantile(align.col.by[[cell.colour]], probs = c(col.max.threshold))
      }

    }

  }


  ### Preparat the raster data

  ## Image prep

  raster.image <- spatial.dat[[roi]]$RASTERS[[roi.marker]]

  tiff.p <- rasterToPoints(raster.image)
  tiff.df <- data.frame(tiff.p)
  raster.label <- names(tiff.df)[3]
  colnames(tiff.df) <- c("x_axis", "y_axis", raster.label)

  ## Create cell outlines
  if(!is.null(mask.outlines)){
    outline <- spatial.dat[[roi]]$MASKS[[mask.outlines]]$outlines
    centroids <- spatial.dat[[roi]]$MASKS[[mask.outlines]]$centroids

    centroid.xmin <- centroids@bbox[1]
    centroid.xmax <- centroids@bbox[3]

    centroid.ymin <- centroids@bbox[2]
    centroid.ymax <- centroids@bbox[4]
  }

  ## Flip y-axis values

  # if(flip.y.axis == TRUE){
  #   dat <- invert.y.axis(dat, y.axis)
  # }

  ## Normalise XY for cell centroids
  plot.normalize <- function(dat, min, max){
    return(((dat- min(dat)) / (max(dat)-min(dat))) * (max - min) + min)
  }

  # if(!is.null(cell.dat)){
  #   # X AXIS
  #   cell.dat[[cell.x]] <- plot.normalize(cell.dat[[cell.x]], min = centroid.xmin, max = centroid.xmax)
  #
  #   # Y AXIS
  #   cell.dat[[cell.y]] <- plot.normalize(cell.dat[[cell.y]], min = centroid.ymin, max = centroid.ymax)
  #
  # }

  ## Raster colour limits

  RastMin <- quantile(tiff.df[[3]], probs = c(raster.min.threshold))
  RastMax <- quantile(tiff.df[[3]], probs = c(raster.max.threshold))

  ###############################################
  ### Add a check to see if centroids line up ###
  ###############################################

  ### Generate and show coloured plot

  if(image.blank == FALSE){
    p <- ggplot(data=tiff.df, aes(x=tiff.df[[1]], y=tiff.df[[2]])) +

      ## Plot the raster (IMC image)
      geom_raster(aes(fill=tiff.df[[3]])) +
      scale_fill_gradient(raster.label,
                          low = "black",
                          high = "white",
                          limits=c(RastMin,RastMax),
                          oob=squish)
  }

  if(image.blank == TRUE){
    p <- ggplot(data=tiff.df, aes(x=tiff.df[[1]], y=tiff.df[[2]])) +

      ## Plot the raster (IMC image)
      geom_raster(aes(fill=tiff.df[[3]])) +
      scale_fill_gradient(raster.label,
                          low = "black",
                          high = "black",
                          limits=c(RastMin,RastMax),
                          oob=squish)
  }



  ### Plot the cell mask boundaries

  if(!is.null(image.outlines)){
    p <- p + geom_path(aes(x = long, y = lat, group = group),
                       data = outline,
                       size = raster.mask.size,
                       col = raster.mask.colour)
  }

  ## Plot the cellular data

  if(!is.null(cell.dat)){
    if(cell.colour.type == "numeric"){
      p <- p + geom_point(data=cell.dat,
                          aes(x=cell.dat[[cell.x]], y=cell.dat[[cell.y]], color = cell.dat[[cell.colour]]),  #as.numeric(as.character(col))
                          size = dot.size, #dot.size
                          alpha = dot.alpha # shape = 1
      ) +

        scale_color_gradientn(colours = colour.scheme(50),
                              limits = (c(ColrMin,ColrMax)),
                              oob=squish,
                              name = cell.colour)
    }

    if(cell.colour.type != "numeric"){
      p <- p + geom_point(data=cell.dat,
                          aes(x=cell.dat[[cell.x]], y=cell.dat[[cell.y]], color = as.factor(cell.dat[[cell.colour]])),  #as.numeric(as.character(col))
                          size = dot.size, #dot.size
                          alpha = dot.alpha # shape = 1
      ) +

        scale_colour_discrete(name = cell.colour)
    }

  }

  ## Setup some themes
  p <- p + theme_bw() +
    coord_equal() +
    xlab(cell.x)+
    ylab(cell.y)+
    ggtitle(title)

  ## More themes
  p <- p + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5), # change 'colour' to black for informative axis
                 axis.title.x=element_text(color="Black", face="bold", size=18),
                 axis.title.y=element_text(color="Black", face="bold", size=18),
                 legend.text=element_text(size=12), # large = 30 # small = 8
                 legend.key.height=unit(1,"cm"), # large = 3 # small = 1.2
                 legend.key.width=unit(0.4,"cm"), # large = 1 # small = 0.4
                 #legend.title=element_blank(),
                 plot.title = element_text(color="Black", face="bold", size=16, hjust=0) # size 70 for large, # 18 for small
  )

  if(flip.y.axis == TRUE){

    p <- p + scale_y_reverse()

  }

  ### Save ggplot to disk if desired
  if(save.to.disk == TRUE){
    ggsave(filename = paste0(title, "_ROI_", roi.marker, "_marker_", cell.colour,".png"),
           plot = p,
           path = path,
           width = plot.width,
           height = plot.height,
           limitsize = FALSE)
  }

  print(p)

}





