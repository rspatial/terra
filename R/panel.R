
setMethod("panel", signature(x="SpatRaster"),
	function(x, main, loc.main="topleft", nc, nr, maxnl=16, maxcell=500000, 
		box=FALSE, pax=list(), plg=list(), ...)  {

		dots <- list(...)
		if (!is.null(dots$type)) {
			error("panel", "you cannot set the plot type")
		}
		if (!is.null(dots$breaks)) {
			error("panel", "you cannot use breaks")
		}
		if (any(is.factor(x))) {
			error("panel", "cannot use categorical rasters")
		}

		nl <- max(1, min(nlyr(x), maxnl))

		if (nl==1) {
			out <- plot(x, 1, maxcell=maxcell, main=main[1], ...)
			return(invisible(out))
		}

		nrnc <- .get_nrnc(nr, nc, nl)
		old.par <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(old.par))
		mar <- rep(0.33, 4)
		m <- matrix((1:prod(nrnc))+1, nrow=nrnc[1], ncol=nrnc[2], byrow=TRUE)
		m <- rbind(1, cbind(1, m, max(m)+1), 1)
		w <- c(0.05, rep(0.8/nrnc[2], nrnc[2]), 0.15)
		h <- c(0.05, rep(0.94/nrnc[1], nrnc[1]), 0.05)
		graphics::layout(m, w, h)
		plot(0, type="n", axes=FALSE, xlab="", ylab="")

		maxcell <- maxcell/nl

		if (missing("main")) {
			tm <- time(x)
			if (!any(is.na(tm))) {
				main <- as.character(time(x))
			} else {
				main <- names(x)
			}
		} else {
			main <- rep_len(main, nl)
		}
		legend <- rep(FALSE , nl)
		legi <- max(1, ceiling(nrnc[1] / 2)) * nrnc[2]
		legend[legi] <- TRUE		
		left <- c(0,2)[(((1:nl)-1) %% nrnc[2] == 0)+1]
		b <- ((1:nl) > ((nrnc[1]-1) * nrnc[2]))+1
		i <- (prod(nrnc) - nl)
		if (i > 0) {
			b[((nrnc[1]-2) * nrnc[2]) + ((nrnc[2]-(i-1)):nrnc[2])] <- 2
		}
		bottom <- c(0,1)[b]

		rng <- range(minmax(x, FALSE))
		if (any(is.na(rng))) {
			xs <- spatSample(x, maxcell, ext=ext, method="regular", as.raster=TRUE, warn=FALSE)
			rng <- range(minmax(xs, TRUE))
		}
		if (is.null(plg$size)) plg$size <- max(1, nrnc[1] * 0.75)
		if (is.null(plg$cex))  plg$cex  <- 1.25 
		plg$yshift <- (nrnc[1] %% 2 == 0) 
		for (i in 1:nl) {
			pax$side <- c(bottom[i], left[i])
			plot(x, i, main=main[i], mar=mar, legend=legend[i], range=rng, pax=pax, box=box, 
				loc.main=loc.main, plg=plg, ...)
		}
	}
)
