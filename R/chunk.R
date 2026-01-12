# Author: Robert J. Hijmans
# Date : December 2025
# Version 1.0
# License GPL v3


setMethod("chunk", signature(x="SpatRaster"),
	function(x, fun, ..., n=NULL, buffer=0, filename="", wopt=list()) {

		if (is.null(n)) {
			m <- mem_info(x, n=6, print=FALSE)
			m <- data.frame(t(m))
			m <- pmax(2, round(m$needed / m$available))
			n <- pmax(1, round(max(dim(x)) / m))
		}
		g <- getTileExtents(x, n, extend=FALSE, buffer=buffer)

		buf <- FALSE
		if (all(buffer > 0)) {
			be <- rep(res(x) * buffer, each=2)
			buf <- TRUE
		}

		x <- deepcopy(x)
		tmpdir <- tempfile()
		dir.create(tmpdir, FALSE, FALSE)
		tmpf <- file.path(tmpdir, paste0("tmp", 1:nrow(g), ".tif"))
		rlst <- vector(mode="list", nrow(g))

		nodots <- length(list(...)) == 0

		for (i in 1:nrow(g)) {
			window(x) <- g[i,]
			if (nodots) {
				r <- fun(x)
			} else {
				r <- fun(x, ...)			
			}
			if (buf) {
				rlst[[i]] <- crop(r, ext(r) - be, filename=tmpf[i])	
			} else {
				rlst[[i]] <- writeRaster(r, tmpf[i])
			}
		}
		s <- sprc(rlst)
		m <- merge(s, filename=filename, wopt=wopt)
		try(file.remove(tmpf), silent=TRUE)
		m
	}
)

