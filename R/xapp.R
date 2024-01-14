
setMethod("xapp", signature(x="SpatRaster", y="SpatRaster"),
function(x, y, fun, ..., filename="", overwrite=FALSE, wopt=list())  {

	compareGeom(x, y, crs=FALSE, warncrs=TRUE)
	if (!hasValues(x)) error("xapp", "x does not have values")
	if (!hasValues(y)) error("xapp", "y does not have values")

	fun <- match.fun(fun)
	out <- rast(x)
	nc <- ncol(x)
	readStart(x)
	readStart(y)
	on.exit(readStop(x))
	on.exit(readStop(y), add=TRUE)

	dots <- list(...)
	if (length(dots) > 0) {
		test <- any(sapply(dots, function(i) inherits(i, "SpatRaster")))
		if (test) {
			error("app", "additional arguments cannot be a SpatRaster")
		}
	}
	teststart <- max(1, 0.5 * nc - 6)
	testend <- min(teststart + 12, nc)
	ntest <- 1 + testend - teststart
	vx <- readValues(x, round(0.51*nrow(x)), 1, teststart, ntest, mat=TRUE)
	vy <- readValues(y, round(0.51*nrow(y)), 1, teststart, ntest, mat=TRUE)
	test <- sapply(1:nrow(vx), function(i) fun(vx[i, ], vy[i, ], ...))
	if (is.list(test)) {
		error("xapp", "'fun' returns a list (should be numeric or matrix)")
	}
	trans <- FALSE
	if (NCOL(test) > 1) {
		if (ncol(test) == ntest) {
			nlyr(out) <- nrow(test)
			trans <- TRUE
			nms <- rownames(test)
		} else if (nrow(test) == ntest) {
			nlyr(out) <- ncol(test)
			nms <- colnames(test)
		} else {
			error("xapp", "the number of values returned by 'fun' is not appropriate\n(it should be the product of the number of cells and and a positive integer)")
		}
		if (is.null(wopt$names)) {
			wopt$names <- nms
		}
	} else {
		if ((length(test) %% ntest) != 0) {
			error("xapp", "the number of values returned by 'fun' is not appropriate")
		} else {
			nlyr(out) <- length(test) / ntest
		}
	}

	ncops <- (nlyr(x)+nlyr(y)) / nlyr(out)
	ncops <- ifelse(ncops > 1, ceiling(ncops), 1) * 4
	b <- writeStart(out, filename, overwrite, wopt=wopt, n=ncops, sources=sources(x))
	for (i in 1:b$n) {
		vx <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
		vy <- readValues(y, b$row[i], b$nrows[i], 1, nc, TRUE)
		r <- sapply(1:nrow(vx), function(i) fun(vx[i, ], vy[i, ], ...))
		if (trans) {
			r <- t(r)
		}
		writeValues(out, r, b$row[i], b$nrows[i])
	}
	writeStop(out)
}
)


