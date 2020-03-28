# Author: Robert J. Hijmans
# Date :  June 2017
# Version 1.0
# License GPL v3


setMethod ("show" , "Rcpp_SpatCategories", 
	function(object) {
		print(data.frame(value=object$levels, label=object$labels))
	}
)

setMethod ("show" , "Rcpp_SpatDataFrame", 
	function(object) {
		d <- .getSpatDF(object)
		head(d)
	}
)


setMethod ("show" , "SpatExtent", 
	function(object) {
		e <- as.vector(object)
		cat("class       :" , class(object), "\n")
		cat("xmin        :" , e[1], "\n")
		cat("xmax        :" , e[2], "\n")
		cat("ymin        :" , e[3], "\n")
		cat("ymax        :" , e[4], "\n")
	}
)	
	

setMethod ("show" , "SpatVector", 
	function(object) {
		e <- as.vector(ext(object))
		cat("class       :", class(object), "\n")
		cat("geometry    :", geomtype(object), "\n")
		cat("elements    : ", size(object), "\n", sep="" ) 
		cat("extent      : ", e[1], ", ", e[2], ", ", e[3], ", ", e[4], "  (xmin, xmax, ymin, ymax)\n", sep="")
		cat("coord. ref. :", crs(object), "\n")
		nms <- names(object)
		if (length(nms) > 0) {
			cat("names       :", paste(names(object), collapse=", "), "\n")		
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

		e <- as.vector(ext(object))
		cat("extent      : " , e[1], ", ", e[2], ", ", e[3], ", ", e[4], "  (xmin, xmax, ymin, ymax)\n", sep="")

		crs <- crs(object)
		cat("coord. ref. :" , crs(object), "\n")
		
		mnr <- 5

		ln <- names(object)
		nl <- nlyr(object)
			
		if (nl > mnr) {
			ln <- c(ln[1:mnr], "...")
		}

		if (.hasValues(object)) {
			nsr <- nsrc(object)	
			m <- .inMemory(object)
			f <- .filenames(object)
			f <- gsub("\\", "/", f, fixed=TRUE)
			sources <- rep("memory", length(m))
			sources[!m] <- f[!m] 

			if (nsr > 1) {
				lbs <- .nlyrBySource(object)
				lbsprint <- paste0(" (", lbs, " layers)")
				lbsprint[lbs == 1] <- ""
				cat("data sources:", sources[1], lbsprint[1], "\n")
				for (i in 2:(min(mnr, nsr))) {
					cat("             ", sources[i], lbsprint[i], "\n")
				}			
				
				if (nsr > mnr) {
					cat("             ", "... and", nsr-mnr, "more sources\n")				
				}
			} else {
				cat("data source :", sources[1], "\n")
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
				if (nl > mnr) {
					minv <- c(minv[1:mnr], "...")
					maxv <- c(maxv[1:mnr], "...")
				}
				
				
				n <- nchar(ln)
				if (nl > 5) {
					b <- n > 20
					if (any(b)) {
						mid <- floor(n/2)
						ln[b] <- paste(substr(ln[b], 1, 7), "~", substr(ln[b], nchar(ln[b])-7, nchar(ln[b])), sep="")
					}
				}
				
				w <- pmax(nchar(ln), nchar(minv), nchar(maxv))
				m <- rbind(ln, minv, maxv)
				# a loop because "width" is not recycled by format
				for (i in 1:ncol(m)) {
					m[,i]   <- format(m[,i], width=w[i], justify="right")
				}
				cat("names       :", paste(m[1,], collapse=", "), "\n")
				cat("min values  :", paste(m[2,], collapse=", "), "\n")
				cat("max values  :", paste(m[3,], collapse=", "), "\n")

			} else {
				cat("names       :", paste(ln, collapse=", "), "\n")
			}			
		} else {
			cat("data sources:", "no data\n")
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


