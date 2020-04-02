

setMethod("plot", signature(x="SpatRaster", y="numeric"), 
	function (x, y, col, maxcell = 100000, leg.mar=NULL, leg.levels=5, leg.shrink=c(0,0), leg.main=NULL, leg.main.cex=1, leg.ext=NULL, digits, useRaster = TRUE, zlim, xlab="", ylab="", axes=TRUE, add=FALSE, ...) {

		x <- x[[y[1]]]
		if (!hasValues(x)) {
			warning("SpatRaster has no cell values")
			return()
		}
		if (missing(col)) {
			col <- rev(grDevices::terrain.colors(100))
		}
		if (add) {
			image(x, maxcell=maxcell, col=col, add=TRUE, ...) 
			return(invisible(NULL))
		}

		if (couldBeLonLat(x, warn=FALSE)) {
			asp <- 1/cos((mean(as.vector(ext(x))[3:4]) * pi)/180)
		} else {
			asp <- 1
		}
		
		object <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		
		Y <- yFromRow(object, nrow(object):1)
		Z <- t(as.matrix(object, TRUE)[nrow(object):1, , drop = FALSE])
		X <- xFromCol(object, 1:ncol(object))

		fact = FALSE
		if (missing(zlim)) {
			uvals <- unique(stats::na.omit(as.vector(Z)))
			if (length(uvals) == 0) { return(invisible(NULL)) }
			if (length(uvals) < 10 & (!is.factor(x))) {
				fact = TRUE
				Z[is.nan(Z)] = NA
				fz <- as.factor(Z)
				Z[] = as.numeric(fz)
				lvs <- list(levels=sort(unique(as.vector(Z))), labels=levels(fz)) 
				uvals <- lvs$levels
				zlim <- range(lvs$levels)
			} else {
				zlim <- range(uvals, na.rm=TRUE)
			}
		} else {
			zlim <- sort(zlim)
			Z[Z < zlim[1]] <- zlim[1]
			Z[Z > zlim[2]] <- zlim[2]
			uvals <- sort(unique(stats::na.omit(as.vector(Z))))
			if (length(uvals) == 0) { return(invisible(NULL)) }			
		}
		
#		dv <- dev.list()
#		if (!is.null(dv)) {
#			# excluding pdf, png, ....
#			if (names(dv[length(dv)]) %in% 
#			c("windows", "cairox11", "RStudioGD")){
#				graphics::plot.new()
#			}
#		}

		if (is.null(leg.mar)) {
			if (is.null(leg.ext)) {
				leg.mar=3
			} else {
				leg.mar=0
			}
		}
		old.mar <- graphics::par()$mar 
		on.exit(graphics::par(mar=old.mar))
		leg.hor <- FALSE
		if (leg.hor) {
			graphics::par(mar=.getMar(c(leg.mar, 0, 0, 0)))
		} else {
			graphics::par(mar=.getMar(c(0, 0, 0, leg.mar)))
		}		
		graphics::image(X, Y, Z, col=col, useRaster=useRaster, asp=asp, xlab=xlab, ylab=ylab, axes=axes, ...)

		if (missing(digits)) {
			dif <- diff(zlim)
			if (dif == 0) {
				digits = 0;
			} else {
				digits <- max(0, -floor(log10(dif/10)))
			}
		}
		
		usr <- graphics::par()$usr
		dx <- graphics::par()$cxy[1] * graphics::par("cex")	

		if (!is.null(leg.ext)) {
			xex <- as.vector(leg.ext)
			leg.ext <- .getLegCoords(NULL, xex, leg.shrink, leg.main)
			leg.ext.set <- TRUE
		} else {
			p <- c(usr[2]+dx, usr[2]+2*dx, usr[3], usr[4])
			xex <- as.vector(ext(object))
			leg.ext <- .getLegCoords(p, xex, leg.shrink, leg.main)
			leg.ext.set <- FALSE
		} 


		if (is.factor(x)) {
			lvs <- levels(x)[[1]]
			if (length(lvs$labels) > 0) {
				fact = TRUE
			}
		}
		if (fact) {
			levs <- lvs$levels
			labs <- lvs$labels
			i <- levs %in% uvals
			levs <- levs[i]
			labs <- labs[i]
			n <- ifelse(leg.ext.set, length(labs), 20)
			col <- .sampleColors(col, length(labs))
			.factorLegend(leg.ext, levs, col, labs, n)
		} else {			
			.contLegend(leg.ext, col, zlim, digits, leg.levels)
		}	
		.legMain(leg.main, leg.ext$xmax, leg.ext$ymax, leg.ext$dy, leg.main.cex)
			
#		setHook("before.plot.new", 
#			function(...) {
#				m <- graphics::par()$mar
#				m[4] <- mars[4]
#				graphics::par(mar=m)
#				setHook("before.plot.new", NULL, action="replace")
#			}, action="replace")		
	}
)



setMethod("plot", signature(x="SpatRaster", y="missing"), 
	function(x, y, maxcell=100000, nc, nr, main, maxnl=16, add=FALSE, ...)  {
		nl <- min(nlyr(x), maxnl)
		if (nl == 0) {
			stop("SpatRaster has no cell values")
		}
		if (add) {
			plot(x, 1, maxcell=maxcell, add=TRUE, ...) 
			return(invisible(NULL))
		}

		if (nl==1) {
			if (missing(main)) {
				plot(x, 1, maxcell=maxcell, ...)
			} else {
				plot(x, 1, maxcell=maxcell, main=main, ...)
			}
			return(invisible(NULL))
		}
		if (missing(nc)) {
			nc <- ceiling(sqrt(nl))
		} else {
			nc <- max(1, min(nl, round(nc)))
		}
		if (missing(nr)) {
			nr <- ceiling(nl / nc)
		} else {
			nr <- max(1, min(nl, round(nr)))
			nc <- ceiling(nl / nr)
		}
		
		old.par <- graphics::par(no.readonly = TRUE) 
		on.exit(graphics::par(old.par))
		graphics::par(mfrow=c(nr, nc), mar=c(2, 2, 2, 4))

		maxcell=maxcell/(nl/2)
			
		if (missing("main")) {
			main <- names(x)
		} else {
			main <- rep_len(main, nl)	
		}
		x <- spatSample(x, maxcell, method="regular", as.raster=TRUE)
		if (onelegend) {
		
		} else {
			for (i in 1:nl) {
			#	image(x[[i]], main=main[i], ...)
				plot(x, i, main=main[i], ...)
			}
		}
	}
)



setMethod("lines", signature(x="SpatRaster"),
function(x, mx=50000, ...) {
	if(prod(dim(x)) > mx) {
		stop('too many lines')
	}
	v <- as.polygons(x)
	lines(v, ...)
}
)
