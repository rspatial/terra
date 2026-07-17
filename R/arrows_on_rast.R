# Visualize Directional Arrows on a SpatRaster

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
    
 