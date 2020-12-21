# Author: Robert J. Hijmans
# Date :  June 2017
# Version 1.0
# License GPL v3


setMethod ("show" , "Rcpp_SpatCategories", 
	function(object) {
		print(data.frame(value=object$levels, label=object$labels))
	}
)


printDF <- function(object, n=6, first=FALSE) {
	n <- min(nrow(object), max(n, 0))
	old <- dim(object)
	if (old[2] == 0) { return() }
	h <- object
	if (old[1] > 0) {
		h <- h[1:n, ,drop=FALSE]
	}
	if (old[2] > 12) {
		h <- h[, 1:10]
	}
	d <- dim(object)

	cls <- sapply(object, class)
	cls <- gsub("integer", "int", cls)
	cls <- gsub("numeric", "num", cls)
	cls <- gsub("character", "chr", cls)
	cls <- paste0("<", cls, ">")
	cls <- data.frame(rbind(class=cls))
	names(cls) <- NULL

	nms <- colnames(h)
	nc <- nchar(nms)
	mx <- max(15, 100/d[2])
	nms[nc > mx] <- paste0(substr(nms[nc > mx], 1, (mx-1)), "~")
	if (d[1] > 0) {
		for (i in 1:ncol(h)) {
			if (is.character(h[[i]])) {
				n <- nchar(h[[i]])
				j <- (n > 11 & n > nc[i])
				h[[i]][j] <- paste0(substr(h[[i]][j], 1, (mx-1)), "~")
			} else if (is.numeric(h[[i]])) {
				h[[i]] <- formatC(h[[i]])		
			}
		}
	}
	
	h <- rbind(h[1,,drop=FALSE], h)
	h[1,] <- cls
	if (nrow(h) < d[1]) {
		h <- rbind(h, "...")
	}
	if (first) {
		h <- cbind("", h)
		colnames(h)[1] <- "names       :"
		h[1,1] <- "type        :"
		if (d[1] > 0) {
			h[2,1] <- "values      :"
		}
	}
	if (old[2] > d[2]) {
		name <- paste0("(and ", old[2] - d[2], " more)")
		h[[name]] <- ""
	}
	print(h, row.names = FALSE)
}


setMethod ("show" , "Rcpp_SpatDataFrame", 
	function(object) {
		cat("class       :" , class(object), "\n")
		object <- .getSpatDF(object)
		d <- dim(object)
		cat("dimensions  : ", d[1], ", ", d[2], "  (nrow, ncol)\n", sep="" ) 
		n <- 6
		if (d[1] > 6) {
			cat("values (head)\n") 
		} else {
			cat("values\n") 		
		}
		printDF(object)
	}
)


setMethod ("show" , "SpatExtent", 
	function(object) {
		e <- as.vector(object)
		e <- paste(e, collapse=", ")
		cat("SpatExtent :", e, "(xmin, xmax, ymin, ymax)\n")
	}
)	

setMethod ("show" , "SpatVectorCollection", 
	function(object) {
		cat(" class       :", class(object), "\n")
		cat(" size        :", length(object), "\n")
		for (i in 1:length(object)) {
			v <- object[i]
			cat("              ", geomtype(v), "\n")
		}
	}
)

setMethod ("show" , "SpatVector", 
	function(object) {
		e <- as.vector(ext(object))
		cat(" class       :", class(object), "\n")
		cat(" geometry    :", geomtype(object), "\n")
		d <- dim(object)
		cat(" dimensions  : ", d[1], ", ", d[2], "  (geometries, attributes)\n", sep="" ) 
		cat(" extent      : ", e[1], ", ", e[2], ", ", e[3], ", ", e[4], "  (xmin, xmax, ymin, ymax)\n", sep="")
		cat(" coord. ref. :", .proj4(object), "\n")
		if (all(d > 0)) {
			dd <- as.data.frame(object)
			printDF(dd, 3, TRUE)
		}
	}
)


setMethod ("show" , "SpatRaster", 
	function(object) {
		
		cat("class       :" , class(object), "\n")

		d <- dim(object)
		cat("dimensions  : ", d[1], ", ", d[2], ", ", d[3], "  (nrow, ncol, nlyr)\n", sep="" ) 
		#cat ("ncell       :" , ncell(object), "\n")

		xyres <- res(object)
		cat("resolution  : " , xyres[1], ", ", xyres[2], "  (x, y)\n", sep="")
		hw <- window(object)
		if (any(hw)) {
			w <- as.vector(ext(object))
			if (all(hw)) {
				txt <- "window      : "
			} else {
				txt <- "extent (win): "
			}
			cat(txt, w[1], ", ", w[2], ", ", w[3], ", ", w[4], "  (xmin, xmax, ymin, ymax)\n", sep="")		
			#e <- as.vector(object@ptr$source[[1]]$window$full_extent$vector)
			#cat("full extent : " , e[1], ", ", e[2], ", ", e[3], ", ", e[4], "  (xmin, xmax, ymin, ymax)\n", sep="")
		} else {
			e <- as.vector(ext(object))
			cat("extent      : " , e[1], ", ", e[2], ", ", e[3], ", ", e[4], "  (xmin, xmax, ymin, ymax)\n", sep="")
		}
		cat("coord. ref. :" , .proj4(object), "\n")

		mnr <- 6
		ln <- names(object)
		nl <- nlyr(object)
			
		if (nl > mnr) {
			ln <- c(ln[1:mnr], "...")
		}
		lnmx <- 60 / min(mnr, length(ln))
		b <- nchar(ln) > (lnmx+2)
		if (any(b)) {
			mid <- floor(lnmx/2)
			ln[b] <- paste(substr(ln[b], 1, mid), "~", substr(ln[b], nchar(ln[b])-mid+1, nchar(ln[b])), sep="")
		}


		if (hasValues(object)) {
			nsr <- nsrc(object)	
			m <- .inMemory(object)
			f <- .filenames(object)
			#f <- gsub("\\", "/", f, fixed=TRUE)
			f <- gsub("\"", "", basename(f))
			sources <- rep("memory", length(m))
			sources[!m] <- f[!m] 

			if (nsr > 1) {
				mxsrc <- 3
				lbs <- .nlyrBySource(object)
				lbsprint <- paste0(" (", lbs, " layers)")
				lbsprint[lbs == 1] <- ""
				cat("sources     :", sources[1], lbsprint[1], "\n")
				for (i in 2:(min(mxsrc, nsr))) {
					cat("             ", sources[i], lbsprint[i], "\n")
				}			
				
				if (nsr > mxsrc) {
					cat("             ", "... and", nsr-mxsrc, "more source(s)\n")				
				}
			} else {
				cat("source      :", sources[1], "\n")
			}
			varnms <- varnames(object)
			i <- varnms != ""
			if (any(i)) {
				longnms <- longnames(object)
				i <- longnms != ""
				if (any(i)) {
					varnms[i] <- paste0(varnms[i], " (", longnms[i], ")")
				}
				if (nsr == 1) {
					cat("varname     :", varnms[1], "\n")
				} else {
					cat("varnames    :", varnms[1], "\n")
					for (i in 2:(min(nsr, 3))) {
						cat("             ", varnms[i], "\n")
					}
				}
				if (nsr > 3) {
					cat("              ...\n")					
				}
			}
			uts <- units(object)
			hasunits <- !all(uts == "")
			if (nl > mnr) {
				uts <- c(uts[1:mnr], "...")
			}
			
			hMM <- .hasMinMax(object)
			if (any(hMM)) {
				r <- minmax(object)
				minv <- format(r[1,])
				maxv <- format(r[2,])
				minv <- gsub("Inf", " ? ", minv)
				maxv <- gsub("-Inf", "  ? ", maxv)
				minv[!hMM] <- gsub("NaN", " ? ", minv[!hMM])
				maxv[!hMM] <- gsub("NaN", " ? ", maxv[!hMM])
				minv[hw] <- paste0(">", minv[hw])
				maxv[hw] <- paste0(maxv[hw],"<")
				if (nl > mnr) {
					minv <- c(minv[1:mnr], "...")
					maxv <- c(maxv[1:mnr], "...")
				}
				isf <- is.factor(object)
				if (any(isf)) {
					lv <- levels(object)
					for (i in 1:length(isf)) {
						if (i > mnr) break
						if (isf[i]) {
							cats <- lv[[i]]
							minv[i] <- cats$labels[1]
							maxv[i] <- cats$labels[length(cats$labels)]					
						} 
					}
				}
				w <- pmax(nchar(ln), nchar(minv), nchar(maxv), nchar(uts))
				m <- rbind(ln, minv, maxv)
				if (hasunits) m <- rbind(m, uts)
				# a loop because "width" is not recycled by format
				for (i in 1:ncol(m)) {
					m[,i]   <- format(m[,i], width=w[i], justify="right")
				}
				if (ncol(m) == 1) {
					cat("name        :", paste(m[1,], collapse=", "), "\n")
					cat("min value   :", paste(m[2,], collapse=", "), "\n")
					cat("max value   :", paste(m[3,], collapse=", "), "\n")				
				} else {
					cat("names       :", paste(m[1,], collapse=", "), "\n")
					cat("min values  :", paste(m[2,], collapse=", "), "\n")
					cat("max values  :", paste(m[3,], collapse=", "), "\n")
				}
				if (hasunits) cat("unit        :", paste(m[4,], collapse=", "), "\n")

			} else {
				w <- pmax(nchar(ln), nchar(uts))
				m <- rbind(ln, uts)
				for (i in 1:ncol(m)) {
					m[,i]   <- format(m[,i], width=w[i], justify="right")
				}
				if (ncol(m) == 1) {
					cat("name        :", paste(m[1,], collapse=", "), "\n")
				} else {
					cat("names       :", paste(m[1,], collapse=", "), "\n")				
				}
				if (hasunits) cat("unit        :", paste(m[2,], collapse=", "), "\n")
			}			


			if (object@ptr$hasTime) {
				rtim <- range(time(object))
				utim <- unique(rtim)
				if (length(utim) > 1) {
					cat("time        :", paste(rtim, collapse=" to "), "\n")
				} else {
					cat("time        :", as.character(utim), "\n")
				}
			}		
		}
		# else {
		#	cat("data sources:", "no data\n")
		#	cat("names       :", paste(ln, collapse=", "), "\n")
		# }
	}
)



.sources <- function(x) {
	m <- .inMemory(x)
	f <- .filenames(x)
	f <- gsub("\"", "", basename(f))
	i <- grep(":", f)
	if (length(i) > 0) {
		for (j in i) {
			ff <- try(basename( strsplit(f[j], ':')[[1]][1]), silent=TRUE)
			if (!inherits(ff, "try-error")) {
				f[j] <- ff
			}
		}
	}
	sources <- rep("memory", length(m))
	sources[!m] <- f[!m] 
	unique(sources)
}


setMethod("show" , "SpatRasterDataset", 
	function(object) {
		
		cat("class       :" , class(object), "\n")
		ns <- length(object)
		cat("subdatasets :", ns, "\n") 
		if (ns == 0) return()
		
		d <- dim(object)
		cat("dimensions  :", paste(d, collapse=", "), "(nrow, ncol)\n") 
		nss <- nlyr(object)
		if (length(nss) > 10) {
			nss = c(as.character(nss[1:9], "..."))
		}
		cat("nlyr        :", paste(nss, collapse=", "), "\n") 

		
		xyres <- res(object)
		cat("resolution  : " , xyres[1], ", ", xyres[2], "  (x, y)\n", sep="")
		e <- as.vector(ext(object))
		cat("extent      : " , e[1], ", ", e[2], ", ", e[3], ", ", e[4], "  (xmin, xmax, ymin, ymax)\n", sep="")


		cat("coord. ref. :" , .proj4(object[1]), "\n")

		s <- unlist(lapply(1:ns, function(i) .sources(object[i])))
		s <- unique(s)
		cat("source(s)   :", paste(s, collapse=", "), "\n")

		ln <- names(object)
		if (any(ln != "")) {
			cat("names       :", paste(ln, collapse=", "), "\n")
		}
	}
)


setMethod ("head" , "SpatVector", 
	function(x, n=6L, ...) {
		utils::head(as.data.frame(x), n=n, ...)
	}
)


setMethod ("tail" , "SpatVector", 
	function(x, n=6L, ...) {
		utils::tail(as.data.frame(x), n=n,...)
	}
)


setMethod ("head" , "SpatRaster", 
	function(x, n=6L, ...) {
		utils::head(x[1:n], n=n, ...)
	}
)


setMethod ("tail" , "SpatRaster", 
	function(x, n=6L, ...) {
		nc = ncell(x)
		utils::tail(x[(nc-n+1):nc], n=n, ...)
	}
)


