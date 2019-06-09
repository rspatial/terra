

.contLegend <- function(p, zlim, cols=terrain.colors(25), 
    digits=0, leg.levels=5, leg.shrink=c(0.1, 0.1), 
	leg.main = NULL, leg.main.cex = 1, ...) {

    nc <- length(cols)
	xmin = p[1]
	xmax = p[2]
	shrk <- (p[4]-p[3])
	leg.shrink <- rep_len(leg.shrink,2)
	if (!is.null(leg.main)) {
		n <- length(leg.main)		
		leg.shrink[2] <- max(leg.shrink[2], (.05*n)) 
	}
	ymin = p[3] + shrk * leg.shrink[1]
	ymax = p[4] - shrk * leg.shrink[2]

    Y <- seq(ymin, ymax, length.out = nc + 1)
    rect(xmin, Y[-(nc + 1)], xmax, Y[-1], col = rev(cols), border = NA, xpd=TRUE)
    rect(xmin, ymin, xmax, ymax, border ="black", xpd=TRUE)
	
    dx <- xmax - xmin
	dy <- ymax - ymin
	zz <- pretty(zlim, n = (leg.levels + 1))	
	zz <- zz[zz >= zlim[1] & zz <= zlim[2]]
    ypos <- ymin + (zz - zlim[1])/(zlim[2] - zlim[1]) * dy
    segments(xmin, ypos, xmax + dx * 0.25, ypos, xpd=TRUE)
    text(xmax, ypos, formatC(zz, digits = digits, format = "f"), pos=4, xpd=TRUE, ...)
    if (!is.null(leg.main)) {
		n <- length(leg.main)
		ymax = ymax + 0.05 * dy
        for (i in 1:n) {
			text(x = xmax,
				y = ymax + (n-i) * 0.05 * dy,
				labels = leg.main[i], 
				#adj = c(0.5, 0.5), 
				cex = leg.main.cex,
				xpd=TRUE)
		}
    } 	 	
}



setMethod("plot", signature(x="SpatRaster", y="numeric"), 
	function (x, y, cols, maxpixels = 10000, leg.mar=3, 
		leg.levels=5, leg.shrink=c(0,0), leg.main=NULL, leg.main.cex=1,
		useRaster = TRUE, zlim, asp, ...) {
		
		if (missing(asp)) {
			if (couldBeLonLat(x)) {
				e <- as.vector(ext(x))
				ym <- mean(e[3:4])
				asp <- 1/cos((ym * pi)/180)
			} else {
				asp <- 1
			}		
		}
		
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
		par(mar = mars + c(0, 0, 0, leg.mar))
		image(X, Y, Z, col = cols, useRaster = useRaster, asp=asp, ...)
		usr <- par()$usr
		dx <- par()$cxy[1] * par("cex")
		p <- c(usr[2] + 1 * dx, usr[2] + 2 * dx, usr[3], usr[4])
		
		.contLegend(p, zlim, leg.levels=leg.levels, leg.shrink=leg.shrink,
		leg.main=leg.main, leg.main.cex=leg.main.cex)
		
		setHook("before.plot.new", 
			function(...) {
				m <- graphics::par()$mar
				m[4] <- mars[4]
				graphics::par(mar=m)
				setHook("before.plot.new", NULL, action="replace")
			}, action="replace")		
	}
)

plot(x, leg.main=c("hell", "fellow", "citizens"))