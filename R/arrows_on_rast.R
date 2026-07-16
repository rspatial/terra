#' Visualize Directional Arrows on a SpatRaster
#'
#' This function overlays directional arrows on a plotted \code{\link[terra]{SpatRaster}} object, 
#' typically used to represent vector fields such as wind directions, drainage flow directions, or other angular data.
#'
#' @param x A \code{\link[terra]{SpatRaster}} object representing angles or directions (in degrees or radians or flowdirs). 
#' @param x_base Optional. A \code{\link[terra]{SpatRaster}} object to use as a basemap if no plot is already present.
#' @param unit Character. Unit of the angle values in \code{x}. Either \code{"degrees"} or \code{"radians"} or \code{"flowdir"}. See \code{\link[terra]{terrain}} for \code{flowdir} option.
#' @param arrow_length Numeric. Length of the arrows to be plotted.
#' @param clockwise Logical. Default is \code{FALSE}. If \code{TRUE}, angles measured in radians or degrees are considered clockwise. 
#' @param angle_x_axis,angle_h_axis Anti-clockwise angle from x (or horizontal) axis to the reference direction. Default is 0, e.g., in the case of northing wind direction it is often used \code{angle_x_axis=90, unit="deg", clockwise=TRUE}.
#' @param length,... Additional graphical parameters passed to \code{\link[graphics]{arrows}}.
#'
#' @return Invisibly returns the coordinates and vector components used for plotting
#' 
#' @export
#' @importFrom graphics arrows
#' @author Emanuele Cordano
#' @examples
#'
#'
#' f <- system.file("ex/elev_vinschgau.tif", package="terra")
#' r <- rast(f)  #|> aggregate(fact=5,fun=min)
#' d <- terrain(r, "flowdir")
#' dl0 <- flowdirD8lad(r,lambda=0)
#' dl1 <- flowdirD8lad(r,lambda=1)
#' 
#' plot(r,col=terrain.colors(10))
#' 
#' arrows_on_rast(dl0, unit="flowdir",col="blue")
#' arrows_on_rast(dl1, unit="flowdir",col="red")
#' arrows_on_rast(d, unit="flowdir",col="black")
#' ###


#'
#'
#'

arrows_on_rast <- function(x="flowdir",x_base=NULL,unit=c("deg","rad","flowdir"),arrow_length=NA,clockwise=FALSE,
                           angle_from_h_axis=0,angle_from_x_axis=angle_from_h_axis,length=0,...) {
  
  unit <- unit[1]
  if (is.character(x)) if (x[1]=="flowdir") {
    x <- terrain(x_base,"flowdir")
    unit <- "flowdir"
  }
  if (unit=="flowdir") {
    x[x==0] <- NA
    clockwise <- FALSE
    x <- -log(x)/log(2)*(pi/4)
    angle_from_x_axis <- 0 
    
    
  } else if (unit %in% c("deg","degree","degrees")) {
    
    x <- x/180*pi
    angle_from_x_axis <- angle_from_x_axis/180*pi
    
  }
  names(x) <- "angle"
  
  
  if (clockwise) x <- (-x)   
  x <- x+angle_from_x_axis
  xdf <- as.data.frame(x,xy=TRUE)
  names(xdf)[1] <- "x0"
  names(xdf)[2] <- "y0"
  dr <- (xres(x)+yres(x))/4
  
  if (is.na(arrow_length)) {
   xdf$length <-  2*dr/apply(cbind(abs(cos(xdf$angle)),abs(sin(xdf$angle))),FUN=max,MARGIN=1)

  } else {
    
    xdf$length <- arrow_length
  }  
  
  
  xdf$x1 <- xdf$x0+xdf$length*cos(xdf$angle)
  xdf$y1 <- xdf$y0+xdf$length*sin(xdf$angle)
  

  
  out <- arrows(x0=xdf$x0,y0=xdf$y0,x1=xdf$x1,y1=xdf$y1,length=length,...) #### y_centers, x_centers + dx, y_centers + dy, length = arrow_length, col = "black")
  return(out)  
}  
    
 