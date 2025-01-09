
setMethod("panel", signature(x="SpatRaster"),
	function(x, main, loc.main="topleft", nc, nr, maxnl=16, maxcell=500000, 
		box=FALSE, pax=list(), plg=list(), range=NULL, halo=TRUE, type=NULL, ...)  {

		if (!is.null(type)) {
			type <- match.arg(tolower(type), c("classes", "continuous", "interval"))
		}

#		dots <- list(...)
#		if (!is.null(dots$type)) {
#			error("panel", "you cannot set the plot type")
#		}
#		if (!is.null(dots$breaks)) {
#			error("panel", "you cannot use argument 'breaks'")
#		}
		categorical <- FALSE
		if (is.null(type) && any(is.factor(x))) {
			x <- combineLevels(x)
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
			x <- spatSample(x, maxcell, method="regular", as.raster=TRUE, warn=FALSE)
			r <- as.matrix(x)
			if (is.null(range)) {
				if (all(hasMinMax(x))) {
					range <- range(minmax(x, FALSE))
					if (any(is.nan(range) | is.infinite(range))) {
						r <- as.matrix(x)
						r[is.infinite(r)] <- NA
						range <- range(r, na.rm=TRUE)
					}
				} else {
					r[is.infinite(r)] <- NA
					range <- range(r, na.rm=TRUE)
				}
			}
			r <- unique(as.vector(r))
			if (is.null(type)) {
				if (length(r) > 10) {
					type <- "continuous"
				} else {
					type <- "classes"
				}
			}
			if (type == "classes") {
				levs <- data.frame(ID=1:length(r), sort(r))
				colnames(levs)[2] <- names(x)[1]
				x <- categories(x, 0, levs)
				categorical <- TRUE
			}
		}
		if (is.null(plg$size)) plg$size <- max(1, nrnc[1] * 0.66)
		if (is.null(plg$cex))  plg$cex  <- 1.25 
		plg$yshift <- (nrnc[1] %% 2 == 0) 

		for (i in 1:nl) {
			pax$side <- c(bottom[i], left[i])
			if (categorical) {
				plot(x[[i]], 1, main=main[i], mar=mar, legend=legend[i], pax=pax, box=box, 
					loc.main=loc.main, halo=halo, plg=plg, type="classes", all_levels=TRUE, maxcell=Inf, ...)
			} else {
				plot(x, i, main=main[i], mar=mar, legend=legend[i], range=range, pax=pax, box=box, 
					loc.main=loc.main, halo=halo, plg=plg, type=type, maxcell=Inf, ...)
			}
		}
	}
)

