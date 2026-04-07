# Author: Robert J. Hijmans
# Date :  June 2017
# Version 1.0
# License GPL v3


setMethod("show", "SpatExtent", function(object) {
#	.show_ext(object)
	cat(object@pntr$show())
} )

setMethod("show", "SpatRaster", function(object) cat(object@pntr$show(TRUE)) )

setMethod("show", "SpatRasterCollection", function(object) cat(object@pntr$show()) )

setMethod("show", "SpatRasterDataset", function(object) cat(object@pntr$show()) )

setMethod("show", "SpatVector", function(object) cat(object@pntr$show()) )

setMethod("show", "SpatVectorCollection", function(object) cat(object@pntr$show()) )

setMethod("show", "SpatVectorProxy", function(object) cat(object@pntr$show()) )


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
	cls <- gsub("factor", "fact", cls)
	cls <- paste0("<", cls, ">")
	cls <- data.frame(rbind(class=cls), stringsAsFactors=FALSE)
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
	x <- data.frame(lapply(x, as.character), check.names=FALSE, stringsAsFactors=FALSE)
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


setMethod("show" , "SpatGraticule",
	function(object) {
		cat("class       :" , class(object), "\n")
		v <- vect()
		v@pntr <- object@pntr
		cat("lon         :" , stats::na.omit(v$lon), "\n")		
		cat("lat         :" , stats::na.omit(v$lat), "\n")		
		cat("coord. ref. :", .name_or_proj4(v), "\n")
		e <- as.vector(ext(v))
		cat("extent      : ", e[1], ", ", e[2], ", ", e[3], ", ", e[4], "  (xmin, xmax, ymin, ymax)\n", sep="")
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



str.SpatRaster <- function(object, ...) {
	cat("S4 class 'SpatRaster' [package \"terra\"]\n")
}
setMethod("str", signature(object="SpatRaster"), str.SpatRaster)


str.SpatVector <- function(object, ...) {
	cat("S4 class 'SpatVector' [package \"terra\"]\n")
}
setMethod("str", signature(object="SpatVector"), str.SpatVector)


str.SpatExtent <- function(object, ...) {
	cat("S4 class 'SpatExtent' [package \"terra\"]\n")
}
setMethod("str", signature(object="SpatExtent"), str.SpatExtent)

str.SpatGraticule <- function(object, ...) {
	cat("S4 class 'SpatGraticule' [package \"terra\"]\n")
}
setMethod("str", signature(object="SpatGraticule"), str.SpatGraticule)
