
setMethod("panel", signature(x="SpatRaster"),
	function(x, main, loc.main="topleft", nc, nr, maxnl=16, maxcell=500000, 
		box=FALSE, pax=list(), plg=list(), range=NULL, ...)  {

		dots <- list(...)
		if (!is.null(dots$type)) {
			error("panel", "you cannot set the plot type")
		}
		if (!is.null(dots$breaks)) {
			error("panel", "you cannot use argument 'breaks'")
		}
		categorical <- FALSE
		if (any(is.factor(x))) {
		    lv <- levels(x)
			lv <- lv[sapply(lv, is.data.frame)]
			lv <- try(do.call(rbind, lv), silent=TRUE)
			if (inherits(lv, "try-error")) {
				error("panel", "cannot use non-matching categorical rasters")
			}
			lv <- unique(lv)
			if (length(unique(lv[,1])) < nrow(lv)) {
				error("panel", "cannot use rasters with conflicting categories")			
			}
			categorical <- TRUE
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

		if (!categorical) {
			if (is.null(range)) {
				if (all(hasMinMax(x))) {
					range <- range(minmax(x, FALSE))
				} else {
					x <- spatSample(x, maxcell, method="regular", as.raster=TRUE, warn=FALSE)
					range <- range(minmax(x, TRUE))
				}
			}
			if (diff(range) > 0) {
				ptype <- "continuous"
			} else {
				ptype <- "classes"
			}
		}
		if (is.null(plg$size)) plg$size <- max(1, nrnc[1] * 0.66)
		if (is.null(plg$cex))  plg$cex  <- 1.25 
		plg$yshift <- (nrnc[1] %% 2 == 0) 
		for (i in 1:nl) {
			pax$side <- c(bottom[i], left[i])
			if (categorical) {
				y <- x[[i]]
				levels(y) <- lv
				plot(y, 1, main=main[i], mar=mar, legend=legend[i], pax=pax, box=box, 
					loc.main=loc.main, plg=plg, type="classes", ...)
			} else {
				plot(x, i, main=main[i], mar=mar, legend=legend[i], range=range, pax=pax, box=box, 
					loc.main=loc.main, plg=plg, type=ptype, ...)
			}
		}
	}
)
