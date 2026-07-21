# Visualize Directional Arrows on a SpatRaster

as.arrows <- function(x, unit=c("degrees", "radians", "flowdir"), 
			arrow_length=NA, clockwise=FALSE, angle_from_h_axis=0,
			angle_from_x_axis=angle_from_h_axis, length=0, ...) {
  
	dev <- grDevices::dev.cur()
	if (names(dev) == "null device") {
		error("as.arrows", "there is no plot to add arrows to")
	}
  
	unit <- match.arg(tolower(unit[1]), c("degrees","radians", "flowdir"))
	if (unit == "flowdir") {
		x <- classify(x, cbind(0, NA))
		clockwise <- FALSE
		x <- -log(x)/log(2)*(pi/4)
		angle_from_x_axis <- 0		 
	} else if (unit %in% c("deg","degree","degrees")) {
		x <- x/180*pi
		angle_from_x_axis <- angle_from_x_axis/180*pi
	}
	names(x) <- "angle"
	 
	if (clockwise) x <- (-x)	 
	x <- x + angle_from_x_axis
	xdf <- as.data.frame(x, xy=TRUE)
	names(xdf)[1] <- "x0"
	names(xdf)[2] <- "y0"
	dr <- (xres(x)+yres(x))/4
	
	if (is.na(arrow_length)) {
	  max_val <- pmax(abs(cos(xdf$angle)), abs(sin(xdf$angle)))
	  xdf$length <-	2*dr/max_val
	} else {
		xdf$length <- arrow_length
	}	
	xdf$x1 <- xdf$x0+xdf$length*cos(xdf$angle)
	xdf$y1 <- xdf$y0+xdf$length*sin(xdf$angle)
	
 #### y_centers, x_centers + dx, y_centers + dy, length = arrow_length, col = "black")	
	graphics::arrows(x0=xdf$x0, y0=xdf$y0, x1=xdf$x1, y1=xdf$y1, length=length,...)
}	
		
 
