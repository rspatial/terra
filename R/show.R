# Author: Robert J. Hijmans
# Date :  June 2017
# Version 1.0
# License GPL v3


printDF <- function(x, n=6, first=FALSE) {
	n <- min(nrow(x), max(n, 0))
	old <- dim(x)
	if (old[2] == 0) { return() }
	if (old[1] > 0) {
		x <- x[1:n, ,drop=FALSE]
	}
	if (old[2] > 12) {
		x <- x[, 1:10]
	}
	d <- dim(x)

	cls <- sapply(x, function(i){ a = class(i); a[length(a)]})
	cls <- gsub("integer", "int", cls)
	cls <- gsub("numeric", "num", cls)
	cls <- gsub("character", "chr", cls)
	cls <- paste0("<", cls, ">")
	cls <- data.frame(rbind(class=cls))
	names(cls) <- NULL

	nms <- colnames(x)
	nc <- nchar(nms)
	mx <- max(15, 100/d[2])
	i <- nc > (mx+2)
	nms[i] <- paste0(substr(nms[i], 1, (mx-1)), "~")
	if (d[1] > 0) {
		for (i in 1:ncol(x)) {
			if (is.character(x[[i]])) {
				x[[i]][is.na(x[[i]])] <- "NA"
				n <- nchar(x[[i]])
				j <- n > (mx+2)
				x[[i]][j] <- paste0(substr(x[[i]][j], 1, (mx-1)), "~")
			} else if (is.numeric(x[[i]])) {
				x[[i]] <- formatC(x[[i]])
			}
		}
	}
	x <- data.frame(lapply(x, as.character), check.names=FALSE)	
	x <- rbind(x[1,,drop=FALSE], x)
	x[1,] <- cls
	if (nrow(x) < d[1]) {
		x <- rbind(x, "...")
	}
	if (first) {
		x <- data.frame("", x, check.names=FALSE, stringsAsFactors=FALSE)
		colnames(x)[1] <- "names       :"
		x[1,1] <- "type        :"
		if (d[1] > 0) {
			x[2,1] <- "values      :"
		}
	}
	if (old[2] > d[2]) {
		name <- paste0("(and ", old[2] - d[2], " more)")
		x[[name]] <- ""
	}

	print(x, row.names = FALSE)
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

setMethod ("show" , "Rcpp_SpatCategories", 
	function(object) {
		show(object$df)
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
		cat(" length      :", length(object), "\n")
		for (i in 1:length(object)) {
			v <- object[i]
			if (i==1) {
				cat(" geometry    : ", geomtype(v), " (", nrow(v) , ")\n", sep="")
			} else {
				cat("               ", geomtype(v), " (", nrow(v) , ")\n", sep="")
			}
		}
	}
)

setMethod ("show" , "SpatVector", 
	function(object) {
		e <- as.vector(ext(object))
		d <- dim(object)
		cat(" class       :", class(object), "\n")
		cat(" geometry    :", geomtype(object), "\n")
		cat(" dimensions  : ", d[1], ", ", d[2], "  (geometries, attributes)\n", sep="" ) 
		cat(" extent      : ", e[1], ", ", e[2], ", ", e[3], ", ", e[4], "  (xmin, xmax, ymin, ymax)\n", sep="")
		if (object@ptr$source != "") {
			if (object@ptr$layer != tools::file_path_sans_ext(basename(object@ptr$source))) {
				cat(" source      : ", basename(object@ptr$source), " (", object@ptr$layer, ")\n", sep="")
			} else {
				cat(" source      : ", basename(object@ptr$source), "\n", sep="")
			}
		}
		cat(" coord. ref. :", .name_or_proj4(object), "\n")
		if (d[2] > 0) {
			nr <- min(d[1], 3)
			dd <- as.data.frame(object[1:nr,])
			printDF(dd, 3, TRUE)
		}
	}
)


setMethod ("show" , "SpatVectorProxy", 
	function(object) {
		e <- as.vector(ext(object))
		d <- dim(object)
		cat(" class       : SpatVectorProxy\n")
		cat(" geometry    :", geomtype(object), "\n")
		cat(" dimensions  : ", d[1], ", ", d[2], "  (geometries, attributes)\n", sep="" ) 
		cat(" extent      : ", e[1], ", ", e[2], ", ", e[3], ", ", e[4], "  (xmin, xmax, ymin, ymax)\n", sep="")
		if (object@ptr$v$layer != tools::file_path_sans_ext(basename(object@ptr$v$source))) {
			cat(" source      : ", basename(object@ptr$v$source), " (", object@ptr$v$layer, ")\n", sep="")
		} else {
			cat(" source      : ", basename(object@ptr$v$source), "\n", sep="")
		}
		cat(" coord. ref. :", .name_or_proj4(object), "\n")
		dd <- get.data.frame(object)
		printDF(dd, 0, TRUE)
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


		cat("coord. ref. :" , .name_or_proj4(object), "\n")


		if (hasValues(object)) {

			mnr <- 6
			ln <- names(object)
			nl <- d[3]

			if (nl > mnr) {
				ln <- c(ln[1:mnr], "...")
			}
			lnmx <- 60 / min(mnr, length(ln))
			b <- nchar(ln) > (lnmx+2)
			if (any(b)) {
				mid <- floor(lnmx/2)
				ln[b] <- paste(substr(ln[b], 1, mid), "~", substr(ln[b], nchar(ln[b])-mid+1, nchar(ln[b])), sep="")
			}

			nsr <- nsrc(object)
			m <- inMemory(object)
			f <- sources(object)
			hdf5 <- substr(f, 1, 5) == "HDF5:"
			f[!hdf5] <- basename(f[!hdf5])
			if (any(hdf5)) {
				ff <- strsplit(f[hdf5], "://")
				ff <- sapply(ff, function(i) paste(basename(i), collapse="://"))
				ff <- gsub('\"', "", ff)
				f[hdf5] <- ff
			}
			#f <- gsub("\\", "/", f, fixed=TRUE)

			f <- gsub("\"", "", f)
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
			rgb <- RGB(object)
			if (!is.null(rgb)) {				 
				cat(paste("colors", toupper(object@ptr$rgbtype), " :"), paste(rgb, collapse=", "), "\n")
			}
			hasct <- object@ptr$hasColors()
			if (any(hasct)) {
				cat("color table :", paste(which(hasct), collapse=", "), "\n")			
			}

		
			varnms <- varnames(object)
			fnms <- tools::file_path_sans_ext(f)
			if (any(fnms != varnms) && all(varnms != "")) {
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

			hMM <- hasMinMax(object)
			isB <- is.bool(object)
			if (any(hMM) || any(is.factor(object))) {
				#r <- minmax(object)
				r <- rbind(object@ptr$range_min, object@ptr$range_max)
				r[,!hMM] <- c(Inf, -Inf)
				minv <- format(r[1,])
				maxv <- format(r[2,])
				if (any(isB)) {
					minv[isB] <- ifelse(minv[isB]=="0", "FALSE", "TRUE")
					maxv[isB] <- ifelse(maxv[isB]=="0", "FALSE", "TRUE")
				}
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
					for (i in 1:length(isf)) {
						if (i > mnr) break
						if (isf[i]) {
							cats <- cats(object, i, TRUE)[[1]]
							j <- match(r[,i], cats[,1])
							cats <- cats[j, 2]
							if (length(cats) > 0) {
								minv[i] <- cats[1]
								maxv[i] <- cats[2]
							}
						} 
					}
				}
				u8 <- Encoding(ln) == "UTF-8"
				wln <- nchar(ln)
				if (any(u8)) {
					# for Chinese: wln <- wln + u8 * wln
					w <- pmax(wln, nchar(minv), nchar(maxv), nchar(uts), na.rm = TRUE)
					m <- rbind(paste0(rep(" ", max(wln)), collapse=""), minv, maxv)
					if (hasunits) m <- rbind(m, uts)
					# a loop because "width" is not recycled by format
					for (i in 1:ncol(m)) {
						m[,i] <- format(m[,i], width=w[i], justify="right")
						addsp <- w[i] - nchar(ln[i])
						m[1,i] <- paste0(paste0(rep(" ", addsp), collapse=""), ln[i])
					}
				} else {
					w <- pmax(wln, nchar(minv), nchar(maxv), nchar(uts), na.rm = TRUE)
					m <- rbind(ln, minv, maxv)
					if (hasunits) m <- rbind(m, uts)
					# a loop because "width" is not recycled by format
					for (i in 1:ncol(m)) {
						m[,i] <- format(m[,i], width=w[i], justify="right")
					}
				}
				if (ncol(m) == 1) {
					if (is.factor(object)) {
						g <- cats(object)[[1]]
						cat("categories  :", paste(colnames(g)[-1], collapse=", "), "\n")
					}
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

		}

		if (object@ptr$hasTime) {
			label <- "time        "
			rtim <- range(time(object))
			tims <- object@ptr$timestep
			if (tims == "yearmonths") {
				rtim <- format_ym(rtim)
				label <- "time (ymnts)"
			} else if (tims == "months") {
				rtim <- month.abb[rtim]
				label <- "time (mnts) "
			} else if (tims == "years") {
				label <- "time (years)"
			} else if (tims == "days") {
				label <- "time (days) "
			} else if (tims == "raw") {
				label <- "time (raw)  "
			}
			utim <- unique(rtim)
			if (length(utim) > 1) {
				ptim <- paste0(label, ": ", paste(rtim, collapse=" to "))
			} else {
				ptim <- paste0(label, ": ", as.character(utim))
			}
			if (tims == "seconds") {
				tz <- format(utim[1], format="%Z")
				ptim <- paste(ptim, tz)
			}
			cat(ptim, "\n")
		}

		# else {
		#	cat("data sources:", "no data\n")
		#	cat("names       :", paste(ln, collapse=", "), "\n")
		# }
	}
)



.sources <- function(x) {
	m <- inMemory(x)
	f <- sources(x)
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


		cat("coord. ref. :" , .name_or_proj4(object[1]), "\n")

		s <- unlist(lapply(1:ns, function(i) .sources(object[i])))
		s <- unique(s)
		cat("source(s)   :", paste(s, collapse=", "), "\n")

		ln <- names(object)
		if (any(ln != "")) {
			cat("names       :", paste(ln, collapse=", "), "\n")
		}
	}
)




setMethod("show" , "SpatRasterCollection", 
	function(object) {
		cat("class       :" , class(object), "\n")
		nr <- length(object)
		cat("length      :", nr, "\n") 
	}
)




setMethod ("head" , "SpatVector", 
	function(x, n=6L, ...) {
		nn <- min(n, nrow(x))
		if (nn > 0) { 
			x <- x[1:nn, ]
		} else {
			x <- x[0,]
		}
		as.data.frame(x)
	}
)


setMethod ("tail" , "SpatVector", 
	function(x, n=6L, ...) {
		nrx <- nrow(x)
		nn <- min(n, nrx)
		if (nn > 0) { 
			start <- nrx - nn + 1
			x <- x[start:nrx, ]
		} else {
			x <- x[0,]
		}
		x <- as.data.frame(x)
		if (nn > 0) { 
			rownames(x) <- start:nrx
		}
		x
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


