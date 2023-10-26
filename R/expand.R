

setMethod("extend", signature(x="SpatExtent"),
function(x, y) {
	if (length(y) == 1) {
		y <- rep(y, 4)
	} else if (length(y) == 2) {
		y <- rep(y, each=2)
	} else if (! length(y) == 4 ) {
		error("extend", 'argument "y" should be a vector of 1, 2, or 4 elements')
	}
	y <- abs(y)
	e <- as.vector(x)
	e[1] <- e[1] - y[1]
	e[2] <- e[2] + y[2]
	e[3] <- e[3] - y[3]
	e[4] <- e[4] + y[4]
	ext(e)
}
)



setMethod("extend", signature(x="SpatRaster"),
function(x, y, snap="near", fill=NA, filename="", overwrite=FALSE, ...) {

	if (!inherits(y, "SpatExtent")) {

		if (is.vector(y)) {
			if (length(y) <= 2) {
				y <- rep_len(round(y), 2)
				stopifnot(all(y >= 0))
				adj <- rev(y) * res(x)
				y <- as.vector(ext(x))
				y[1] <- y[1] - adj[1]
				y[2] <- y[2] + adj[1]
				y[3] <- y[3] - adj[2]
				y[4] <- y[4] + adj[2]
				y <- ext(y)
			} else if (length(y) == 4) {
				y[1] <- y[1] - adj[1]
				y[2] <- y[2] + adj[2]
				y[3] <- y[3] - adj[3]
				y[4] <- y[4] + adj[4]
				y <- ext(y)
			} else {
				error("extend", "if 'y' is a vector it should have 1, 2, or four numbers")
			}
		} else {
			test <- try ( y <- ext(y), silent=TRUE )
			if (inherits(test, "try-error")) {
				error("extend", "cannot get a SpatExtent object from argument y")
			}
		}
	}

	opt <- spatOptions(filename, overwrite, ...)
	x@cpp <- x@cpp$expand(y@cpp, snap[1], fill[1], opt)
	messages(x, "extend")
}
)


