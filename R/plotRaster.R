

setMethod("plot", signature(x="SpatRaster", y="numeric"), 
	function (x, y, cols, maxpixels = 10000, zlim, z.mar=3, useRaster = TRUE, ...) {
		if (missing(cols)) cols = rev(terrain.colors(25))
		object <- sampleRegular(x[[y]], maxpixels)
		Y <- yFromRow(object, nrow(object):1)
		Z <- t(as.matrix(object, TRUE)[nrow(object):1, , drop = FALSE])
		X <- xFromCol(object, 1:ncol(object))

		if (missing(zlim)) {
			zlim <- range(Z, na.rm=TRUE)
		} else {
			Z[Z < zlim[1]] <- zlim[1]
			Z[Z > zlim[2]] <- zlim[2]
		}
		plot.new()
		mars <- par("mar")
		par(mar = mars + c(0, 0, 0, z.mar))
		image(X, Y, Z, col = cols, useRaster = useRaster, ...)
		usr <- par()$usr
		dx <- par()$cxy[1] * par("cex")
		p <- c(usr[2] + 1 * dx, usr[2] + 2 * dx, usr[3], usr[4])
		
		rect(p[1], p[3], p[2], p[4], xpd=TRUE, col="red")
		
		setHook("before.plot.new", 
			function(...) {
				m <- graphics::par()$mar
				m[4] <- mars[4]
				graphics::par(mar=m)
				setHook("before.plot.new", NULL, action="replace")
			}, action="replace")		
	}
)

