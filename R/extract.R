# Author: Robert J. Hijmans
# Date : November 2018
# Version 1.0
# License GPL v3


.big_number_warning <- function() {
# this warning should be given by C
	warn("big number", "cell numbers larger than ", 2^.Machine$double.digits, " are approximate")
}


ext_from_rc <- function(x, r1, r2, c1, c2){
	e <- as.vector(ext(x))
	r <- res(x)
	c1 <- min(max(c1, 1), ncol(x))
	c2 <- min(max(c2, 1), ncol(x))
	if (c1 > c2) {
		tmp <- c1
		c1 <- c2
		c2 <- tmp
	}
	r1 <- min(max(r1, 1), nrow(x))
	r2 <- min(max(r2, 1), nrow(x))
	if (r1 > r2) {
		tmp <- r1
		r1 <- r2
		r2 <- tmp
	}

	xn <- xFromCol(x, c1) - 0.5 * r[1]
	xx <- xFromCol(x, c2) + 0.5 * r[1]
	yx <- yFromRow(x, r1) + 0.5 * r[2]
	yn <- yFromRow(x, r2) - 0.5 * r[2]
	ext(c(sort(c(xn, xx))), sort(c(yn, yx)))
}


getLyrNrs <- function(layer, nms, n) { 
	nl <- length(nms)
	if (is.numeric(layer)) {
		layer <- round(layer)
		if (min(layer, na.rm=TRUE) < 1 || max(layer, na.rm=TRUE) > nl) {
			error("extract", "layer should be between 1 and nlyr(x)")
		}
	} else {
		layer <- match(layer, nms)
	}
	if (any(is.na(layer))) {
#		error("extract", "names in argument 'layer' do not match names(x)")
	}
	rep_len(layer, n)
}


extractCells <- function(x, y, raw=FALSE) {
	opt <- spatOptions()
	e <- x@pntr$extractCell(y-1, opt)
	e <- do.call(cbind, e)
	colnames(e) <- names(x)
	if (!raw) {
		e <- .makeDataFrame(x, e)
	}
	e
}


use_layer <- function(e, y, layer, nl, keepID) {
 	if (is.null(layer)) {
		return(e)
	}
	layer <- getLyrNrs(layer, colnames(e)[-1], nrow(y)) 
	idx <- cbind(1:nrow(e), layer[e[,1]] + 1)
	ee <- data.frame(e[,1,drop=FALSE], layer=colnames(e)[idx[,2]], value=e[idx])
	if (ncol(e) > (nl+1)) {
		ee <- cbind(ee, e[,(nl+2):ncol(e), drop=FALSE])
	} 
	if (!keepID) {
		ee <- ee[,-1]
	}
	ee
}


extract_table <- function(x, y, ID=FALSE, weights=FALSE, exact=FALSE, touches=FALSE, small=TRUE, na.rm=FALSE, ...) {
	if (weights && exact) {
		exact = FALSE
	}
	one_tab <- function(i) {
		e <- extract(x, y[i,], ID=FALSE, weights=weights, exact=exact, touches=touches, small=small, na.rm=na.rm)
		if (nrow(e) == 0) return(NULL)
		if (!(weights | exact)) {
			e$w_e_i_g_h_t123 <- 1
		}
		nms <- names(e)
		a <- aggregate(e[nms[length(nms)]], e[, c(nms[-length(nms)]), drop=FALSE], sum, na.rm=na.rm)
		colnames(a)[ncol(a)] <- "count"
		if (nrow(a) == 0) {
			return(NULL)
			#a <- a[1,]
			#a$count <- NULL 
			#a$count <- NA
		}
		cbind(ID=i, a)
	}
	# use a loop to keep mem consumption low
	out <- lapply(1:nrow(y), function(i) one_tab(i))
	out <- do.call(rbind, out)
	rownames(out) <- NULL
	if ((!ID) & (!is.null(out))) {
		split(out[,-1, drop=FALSE], out[,1])		
	} else {
		out
	}
}



extract_fun <- function(x, y, fun, ID=TRUE, weights=FALSE, exact=FALSE, touches=FALSE, small=TRUE, layer=NULL, bind=FALSE, na.rm=FALSE) {

	nl <- nlyr(x)
	nf <- length(fun)
	if ((nf > 1) & (!is.null(layer))) {
		error("extract", "you cannot use argument 'layer' when using multiple functions")
	}

	opt <- spatOptions()
	e <- x@pntr$extractVectorFlat(y@pntr, fun, na.rm, touches[1], small[1], "", FALSE, FALSE, weights, exact, opt)
	x <- messages(x, "extract")
	e <- data.frame(matrix(e, ncol=nl*nf, byrow=TRUE))
	if (nf == 1) {
		colnames(e) <- names(x)
	} else {
		colnames(e) <- apply(cbind(rep(fun, each=nl), names(x)), 1, function(i) paste(i, collapse="_"))
	}
	if (!is.null(layer)) {
		e <- cbind(ID=1:nrow(e), e)
		e <- use_layer(e, y, layer, nlyr(x), ID)
		if (!ID || bind) {
			e$ID <- NULL
		}
		ID <- FALSE
	} 
	if (bind) {
		if (nrow(e) == nrow(y)) {
			e <- data.frame(e)
			e <- cbind(y, e)
		} else {
			#? can this occur?
			warn("extract", "cannot return a SpatVector because the number of records extracted does not match he number of rows in y (perhaps you need to use a summarizing function")
		}
	} else if (ID) {
		e <- cbind(ID=1:nrow(e), e)
	}
	e
}


do_fun <- function(e, fun, ...) {		
	fun <- match.fun(fun)
	e <- aggregate(e[,-1,drop=FALSE], e[,1,drop=FALSE], fun, ...)
	m <- sapply(e, NCOL)
	if (any(m > 1)) {
		cn <- names(e)
		e  <- do.call(cbind, as.list(e))
		i <- rep(1, length(cn))
		i[m>1] <- m[m>1]
		cn <- rep(cn, i)
		cn <- make.names(cn, TRUE)
		if (length(cn) == ncol(e)) {
			colnames(e) <- cn
		}
	}
	e
}


setMethod("extract", signature(x="SpatRaster", y="SpatVector"),
function(x, y, fun=NULL, method="simple", cells=FALSE, xy=FALSE, ID=TRUE, weights=FALSE, exact=FALSE, touches=is.lines(y), small=TRUE, layer=NULL, bind=FALSE, raw=FALSE, search_radius=0, ...) {

	geo <- geomtype(y)
	if (!is.null(layer)) {
		if (length(layer) > nrow(y)) {
			error("extract", "length(layer) > nrow(y)")
		} else { # recycle
			if (is.numeric(layer)) {
				layer <- round(layer)
				if (min(layer, na.rm=TRUE) < 1 || max(layer, na.rm=TRUE) > nlyr(x)) {
					error("extract", "layer should be between 1 and nlyr(x)")
				}
			}
			#x <- x[[unique(layer)]]
			#if (is.numeric(layer)) { # Match new layer order
			#  layer <- match(layer, unique(layer))
			#}
			layer <- rep(layer, length.out=nrow(y))
		}
		#keepID <- ID
		#ID <- TRUE
	}

	opt <- terra:::spatOptions()

	if (geo == "points") {
		if (search_radius > 0) {
			pts <- crds(y)
			e <- x@pntr$extractBuffer(pts[,1], pts[,2], search_radius, opt)
			messages(x, "extract")
			e <- do.call(cbind, e)
			colnames(e) <- c(names(x)[1], "distance", "cell")		
			e[,3] <- e[,3] + 1
			if (xy) {
				e <- cbind(xyFromCell(x, e[,3]), e)
			}
			if (!raw) {
				e <- cbind(.makeDataFrame(x, e[,1,drop=FALSE]), e[,2:3])
			}
			if (bind) {
				e <- data.frame(e)
				e <- cbind(y, e)
			} else if (ID) {
				e <- cbind(ID=1:nrow(e), e) 
			}
			return(e)
		} else if (weights || exact) {
			method <- "bilinear"
			weights <- FALSE
			exact <- FALSE
		} 
		# method <- match.arg(tolower(method), c("simple", "bilinear"))
	} else if (!is.null(fun)) { # nothing to summarize for points
		txtfun <- .makeTextFun(fun)
		if (inherits(txtfun, "character")) {
			if (any(txtfun == "table")) {
				if (length(fun) > 1) {
					warn("extract", "'table' cannot be combined with other functions")
				}
				if (!is.null(layer)) {
					warn("extract", "argument 'layer' is ignored when 'fun=table'")
				}
				e <- extract_table(x, y, ID=ID, weights=weights, exact=exact, touches=touches, small=small, ...)
			} else {
				e <- extract_fun(x, y, txtfun, ID=ID, weights=weights, exact=exact, touches=touches, small=small, bind=bind, layer=layer, ...)
			}
			return(e)
		} else if (weights || exact) {
			error("extract", "if 'weights' or 'exact' is TRUE, you can only use functions mean, sum, min, max and table")
		}
		xy <- cells <- FALSE
		raw <- TRUE
	} 
	
	if (is.null(layer)) {
		e <- x@pntr$extractVectorFlat(y@pntr, "", FALSE, touches[1], small[1], method, isTRUE(cells[1]), isTRUE(xy[1]), isTRUE(weights[1]), isTRUE(exact[1]), opt)
		x <- messages(x, "extract")
		nc <- nl <- nlyr(x)
		cn <- c("ID", names(x))
	} else {

		lyrs <- unique(layer) 
		e <- rep(NA, length(layer))
		yy <- y[,0]
		for (lyr in lyrs) {
			xlyr <- x[[lyr]]
			i <- which(layer == lyr)
			ylyr <- yy[i]
			e[i] <- xlyr@pntr$extractVectorFlat(ylyr@pntr, "", FALSE, touches[1], small[1], method, isTRUE(cells[1]), isTRUE(xy[1]), isTRUE(weights[1]), isTRUE(exact[1]), opt)
		}
		nc <- nl <- 1
		cn <- c("ID", "value")
		x <- rast(xlyr) # for .makeDataFrame
	}

	if (cells) {
		cn <- c(cn, "cell")
		nc <- nc + 1
	}
	if (xy) {
		cn <- c(cn, "x", "y")
		nc <- nc + 2
	}
	if (weights) {
		cn <- c(cn, "weight")
		nc <- nc + 1
	} else if (exact) {
		cn <- c(cn, "fraction")
		nc <- nc + 1
	}
	if (geo == "points") {
		## this was? should be fixed upstream
		if (nc == nl) {
			e <- matrix(e, ncol=nc)
		} else {
			e <- matrix(e, ncol=nc, byrow=TRUE)
		}
		e <- cbind(1:nrow(e), e)
		if (nrow(e) > nrow(y)) { #multipoint
			g <- geom(y)
			e[,1] <- g[,1]
		}
	} else {
		e <- matrix(e, ncol=nc+1, byrow=TRUE)
	}
	colnames(e) <- cn
	if (!is.null(fun)) {
		e <- as.data.frame(e)
		e <- do_fun(e, fun, ...)
	}
	
	if (cells) {
		cncell <- cn =="cell"
		e[, cncell] <- e[, cncell] + 1
	}

	if (!raw) {
		if (method != "simple") {
			e <- as.data.frame(e)
		} else {
			id <- data.frame(e[,1,drop=FALSE])
			e <- cbind(id, .makeDataFrame(x, e[,-1,drop=FALSE]))
		}
	}
	
#	e <- use_layer(e, y, layer, nl, keepID)

	if (bind) {
		if (nrow(e) == nrow(y)) {
			e <- data.frame(e)
			e <- cbind(y, e[,-1,drop=FALSE])
		} else {
			warn("extract", "cannot return a SpatVector because the number of records extracted does not match the number of rows in y (perhaps you need to use a summarizing function")
		}
	} else if (!ID) {
		e <- e[,-1,drop=FALSE]
	}
	e
})



setMethod("extract", signature(x="SpatRaster", y="sf"),
	function(x, y, fun=NULL, method="simple", cells=FALSE, xy=FALSE, ID=TRUE, weights=FALSE, exact=FALSE, touches=is.lines(y), layer=NULL, bind=FALSE, ...) {
		y <- vect(y)
		extract(x, y, fun=fun, method=method, cells=cells, xy=xy, ID=ID, weights=weights, exact=exact, touches=touches, layer=layer, bind=bind, ...)
	}
)


setMethod("extract", signature(x="SpatRaster", y="data.frame"),
function(x, y, ...) {
	if (ncol(y) != 2) {
		error("extract", "extract expects a 2 column data.frame of x/y or lon/lat coordinates")
	}
	v <- vect(y, colnames(y), quiet=TRUE)
	if ((crs(v) != "") && is.lonlat(v, warn=FALSE) && is.lonlat(x, warn=FALSE)) {
		crs(v) <- NULL
	}
	extract(x, v, ...)
})


setMethod("extract", signature(x="SpatRaster", y="numeric"),
function(x, y, ...) {

	if (isTRUE((length(y)  == 2) && (any((y%%1)!=0)))) {
		warn("extract", "a vector of two decimal values is interpreted as referring to cell numbers, not to coordinates")
	}
	if (isTRUE(any((y < 1) | (y > ncell(x))))) {
		warn("extract", "out of range cell numbers detected")	
	}
	y <- xyFromCell(x, y)
	extract(x, y, ...)
})

setMethod("extract", signature(x="SpatRaster", y="matrix"),
function(x, y, ID=FALSE, ...) {
	extract(x, as.data.frame(y), ID=ID, ...)
})

setMethod("extract", signature(x="SpatRaster", y="SpatExtent"),
function(x, y, ...) {
	extract(x, xyFromCells(x, cells(x, y)), ID=FALSE, ...)
}
)

setMethod("extract", c("SpatVector", "SpatVector"),
function(x, y, count=FALSE) {

	e <- relate(y, x, "coveredby", pairs=TRUE, na.rm=FALSE)
	if (count) {
		if ((geomtype(x) == "polygons") && (geomtype(y) == "points")) {
			tab <- as.data.frame(table(e[,2]))
			i <- match(tab[,1], 1:nrow(x))
			count <- rep(NA, nrow(x))
			count[i] <- tab[,2]
			return(count)
		} else {
			error("extract", "count=TRUE is for point (y) in polygons (x) only")
		}
	}
	if (ncol(x) > 0) {
		d <- as.data.frame(x)
		e <- data.frame(id.y=e[,1], d[e[,2], ,drop=FALSE])
		rownames(e) <- NULL
	} else {
		colnames(e) <- c("id.y", "id.x")
	}
	e
})


setMethod("extract", signature(x="SpatVector", y="matrix"),
function(x, y) {
	stopifnot(ncol(y) == 2)
	.checkXYnames(colnames(y))
	y <- vect(y)
	extract(x, y)
})

setMethod("extract", signature(x="SpatVector", y="data.frame"),
function(x, y) {
	extract(x, as.matrix(y))
})


setMethod("extract", signature(x="SpatRasterCollection", y="ANY"),
function(x, y, ...) {
	lapply(x, function(r) extract(r, y, ...))
}
)

setMethod("extract", signature(x="SpatRasterDataset", y="ANY"),
function(x, y, ...) {
	lapply(x, function(r) extract(r, y, ...))
}
)




extractAlong <- function(x, y, ID=TRUE, cells=FALSE, xy=FALSE, online=FALSE, bilinear=TRUE) { 

	stopifnot(inherits(x, "SpatRaster"))
	if (inherits(y, "sf")) {
		y <- vect(y)
	} else {
		stopifnot(inherits(y, "SpatVector"))
	}
	stopifnot(geomtype(y) == "lines")
	
	spbb <- as.matrix(ext(y))
	rsbb <- as.matrix(ext(x))
	addres <- 2 * max(res(x))
	nlns <- nrow(y)
	res <- vector(mode = "list", length = nlns)

	if (spbb[1,1] > rsbb[1,2] | spbb[1,2] < rsbb[1,1] | spbb[2,1] > rsbb[2,2] | spbb[2,2] < rsbb[2,1]) {
		res <- data.frame(matrix(ncol=nlyr(x)+4, nrow=0))
		colnames(res) <- c("ID", "cell", "x", "y", names(x))
		if (!cells) res$cell <- NULL
		if (!xy) {
			res$x <- NULL
			res$y <- NULL
		}
		if (!ID) res$ID <- NULL
		return(res)
	}
	
	rr <- rast(x)
	g <- data.frame(geom(y))
	
	for (i in 1:nlns) {
		yp <- g[g$geom == i, ]
		nparts <- max(yp$part)
		vv <- NULL
		for (j in 1:nparts) {
			pp <- as.matrix(yp[yp$part==j, c("x", "y"), ])
			for (k in 1:(nrow(pp)-1)) {
				ppp <- pp[k:(k+1), ]
				spbb <- t(ppp)
				if (! (spbb[1,1] > rsbb[1,2] | spbb[1,2] < rsbb[1,1] | spbb[2,1] > rsbb[2,2] | spbb[2,2] < rsbb[2,1]) ) {
					lns <- vect(ppp, "lines")
					rc <- crop(rr, ext(lns) + addres)
					rc <- rasterize(lns, rc, touches=TRUE)
					cxy <- crds(rc)
					v <- cbind(row=rowFromY(rr, cxy[,2]), col=colFromX(rr, cxy[,1]))
					#up or down?
					updown <- c(1,-1)[(ppp[1,2] < ppp[2,2]) + 1]
					rightleft <- c(-1,1)[(ppp[1,1] < ppp[2,1]) + 1]

					v <- v[order(updown*v[,1], rightleft*v[,2]), ]
					vv <- rbind(vv, v)
				}
			} 
		}
		if (!is.null(vv)) {
			cell <- cellFromRowCol(rr, vv[,1], vv[,2])
			res[[i]] <- data.frame(i, cell)
		}
	}
	
	res <- do.call(rbind, res)
	if (is.null(res)) {
		if (xy) {
			res <- data.frame(matrix(ncol=nlyr(x)+4, nrow=0))		
			colnames(res) <- c("ID", "cell", "x", "y", names(x))
		} else {
			res <- data.frame(matrix(ncol=nlyr(x)+2, nrow=0))
			colnames(res) <- c("ID", "cell", names(x))
		}
	} else {
		colnames(res) <- c("ID", "cell")
		if (xy) {
			xycrd <- xyFromCell(x, res$cell)
			method <- "simple"
			if (online) {
				pts <- vect(xycrd, crs="local")
				crs(y) <- "local"
				n <- nearest(pts, y)
				xycrd <- crds(n)
				if (bilinear) method <- "bilinear"
				res <- data.frame(res, xycrd, extract(x, xycrd, method=method))
			} else {
				res <- data.frame(res, xycrd, extract(x, res$cell))
			}
		} else {
			res <- data.frame(res, extract(x, res$cell))		
		}
	}

	if (!cells) res$cell <- NULL
	if (!ID) res$ID <- NULL
	
	res

}


setMethod("extractRange", signature(x="SpatRaster", y="ANY"),
	function(x, y, first, last, lyr_fun=NULL, geom_fun=NULL, ID=FALSE, na.rm=TRUE, bind=FALSE, ...) {
		
		first <- getLyrNrs(first, names(x), NROW(y)) + 1 
		last  <- getLyrNrs(last,  names(x), NROW(y)) + 1
		if (inherits(y, "SpatVector")) {
			e <- extract(x, y, geom_fun, ID=TRUE, na.rm=na.rm, bind=FALSE, ...)
			if (nrow(e) != nrow(y)) {
				error("extractRange", "geom_fun must return a single value for each geometry/layer")
			}
		} else {
			e <- data.frame(extract(x, y, ID=TRUE, ...))
		}
		a <- lapply(1:nrow(e), function(i) e[i, c(first[i]:last[i]), drop=FALSE])
		if (!is.null(lyr_fun)) {
			a <- sapply(a, lyr_fun, na.rm=na.rm)
		}
		if (isTRUE(ID)) {
			if (is.list(a)) {
				a <- lapply(1:length(a), function(i) cbind(ID=i, a[[i]]))
			} else {
				a <- data.frame(ID=1:NROW(a), value=a)
			}
		}
		if (isTRUE(bind)) {
			if (is.list(a)) {
				warn("extractRange", "cannot bind these values")
			} else {
				if (is.vector(a)) {
					a <- data.frame(value=a)
				}
				a <- cbind(y, a)
			}
		}
		a
	}
)

