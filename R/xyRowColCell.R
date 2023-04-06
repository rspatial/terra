# Author: Robert J. Hijmans
# Date :  October 2017
# Version 1.0
# License GPL v3

point_on_border <- function(r, x, y, tolerance = sqrt(.Machine$double.eps)) {
    v <- h <- (x >= xmin(r)) & (x <= xmax(r)) & (y >= ymin(r)) & (y <= ymax(r))
    v[v] <- ((x[v] - xmin(r)) %% res(r)[1]) < tolerance
    h[h] <- ((y[h] - ymin(r)) %% res(r)[2]) < tolerance
    h | v
}


setMethod("yFromRow", signature(object="SpatRaster", row="numeric"),
	function(object, row) {
		object@ptr$yFromRow(row - 1)
	}
)

setMethod("yFromRow", signature(object="SpatRaster", row="missing"),
	function(object, row) {
	  row <- seq_len(object@ptr$nrow())
		object@ptr$yFromRow(row - 1)
	}
)

setMethod(xFromCol, signature(object="SpatRaster", col="numeric"),
	function(object, col) {
		object@ptr$xFromCol(col - 1)
	}
)
setMethod(xFromCol, signature(object="SpatRaster", col="missing"),
	function(object, col) {
	  col <- seq_len(object@ptr$ncol())
		object@ptr$xFromCol(col - 1)
	}
)

setMethod(colFromX, signature(object="SpatRaster", x="numeric"),
	function(object, x)	{
		cols <- object@ptr$colFromX(x) + 1
		cols[cols==0] <- NA
		cols
	}
)

setMethod(rowFromY, signature(object="SpatRaster", y="numeric"),
	function(object, y)	{
		rows <- object@ptr$rowFromY(y) + 1
		rows[rows==0] <- NA
		rows
	}
)

setMethod(cellFromXY, signature(object="SpatRaster", xy="matrix"),
	function(object, xy) {
		stopifnot(ncol(xy) == 2)
		#.checkXYnames(colnames(xy))
		object@ptr$cellFromXY(xy[,1], xy[,2], NA) + 1
	}
)

setMethod(cellFromXY, signature(object="SpatRaster", xy="data.frame"),
	function(object, xy) {
		stopifnot(ncol(xy) == 2)
		#.checkXYnames(colnames(xy))
		object@ptr$cellFromXY(xy[,1], xy[,2], NA) + 1
	}
)



setMethod(cellFromRowCol, signature(object="SpatRaster", row="numeric", col="numeric"),
	function(object, row, col) {
		object@ptr$cellFromRowCol(row-1, col-1) + 1
	}
)

setMethod(cellFromRowColCombine, signature(object="SpatRaster", row="numeric", col="numeric"),
	function(object, row, col) {
		object@ptr$cellFromRowColCombine(row-1, col-1) + 1
	}
)

setMethod(rowColCombine, signature(object="SpatRaster", row="numeric", col="numeric"),
	function(object, row, col) {
		cell <- object@ptr$cellFromRowColCombine(row-1, col-1)
		rc <- object@ptr$rowColFromCell(cell)
		rc <- do.call(cbind, rc)
		rc[rc < 0] <- NA
		rc+1
	}
)


setMethod(xyFromCell, signature(object="SpatRaster", cell="numeric"),
	function(object, cell) {
		xy <- object@ptr$xyFromCell(cell-1)
		xy <- do.call(cbind, xy)
		colnames(xy) <- c("x", "y")
		xy
	}
)

setMethod(yFromCell, signature(object="SpatRaster", cell="numeric"),
	function(object, cell) {
		xyFromCell(object, cell)[,2]
	}

)

setMethod(xFromCell, signature(object="SpatRaster", cell="numeric"),
	function(object, cell) {
		xyFromCell(object, cell)[,1]
	}
)

setMethod(rowColFromCell, signature(object="SpatRaster", cell="numeric"),
	function(object, cell) {
		rc <- object@ptr$rowColFromCell(cell-1)
		rc <- do.call(cbind, rc)
		rc[rc < 0] <- NA
		rc+1
	}
)

setMethod(rowFromCell, signature(object="SpatRaster", cell="numeric"),
	function(object, cell) {
		rowColFromCell(object, cell)[,1]
	}
)

setMethod(colFromCell, signature(object="SpatRaster", cell="numeric"),
	function(object, cell) {
		rowColFromCell(object, cell)[,2]
	}
)


.rep_fun <- function(v, n, N, m) {
	if (is.null(v)) {
		v
	} else if (n == 1) {
		replicate(N, v, FALSE)
	} else if (n == N) {
		as.list(v)
	} else if ((n == m) && all(v==(1:m))) {
		NULL
	} else {
		error("rcl", "if one argument is a list, the others should be a list,\n       or a vector of length 1, or have the length of the list argument")
	}
}


setMethod(rcl, signature(x="SpatRaster"),
	function(x, row=NULL, col=NULL, lyr=NULL) {

		lr <- is.list(row)
		lc <- is.list(col)
		ll <- is.list(lyr)
		ln <- sum(lr, lc, ll)
		if (ln > 0) {
			nr <- length(row)
			nc <- length(col)
			nl <- length(lyr)
			N <- unique(c(nr, nc, nl)[c(lr, lc, ll)])
			if (length(N) > 1) {
				N <- max(N)
				islst <- which(c(lr, lc, ll))
				for (i in islst) {
					if (i == 1) {
						row <- rep_len(row, N)
					} else if (i == 2) {
						col <- rep_len(col, N)
					} else {
						lyr <- rep_len(lyr, N)
					}
				}
			}
			notlst <- which(!c(lr, lc, ll))
			for (i in notlst) {
				if (i == 1) {
					row <- .rep_fun(row, nr, N, nrow(x))
				} else if (i == 2) {
					col <- .rep_fun(col, nc, N, ncol(x))
				} else {
					lyr <- .rep_fun(lyr, nl, N, nlyr(x))
				}
			}
			out <- lapply(1:N, function(i) { rcl(x, row[[i]], col[[i]], lyr[[i]]) })
			do.call(rbind, out)
		} else {
			hr <- !is.null(row)
			hc <- !is.null(col)
			hl <- !is.null(lyr)
			n <- sum(hr, hc, hl)

			if (hl) {
				lyr[!c(lyr %in% 1:nlyr(x))] <- NA 
			}
			if (n == 0) {
				out <- rowColFromCell(x, 1:ncell(x))
				out <- cbind(out[,1], out[,2], rep(1:nlyr(x), each=nrow(out)))
				colnames(out) <- c("row", "col", "lyr")
				return(out)
			} else if (hr & hc & hl) {
				out <- rowColCombine(x, row, col)
				out <- cbind(out[,1], out[,2], rep(lyr, each=nrow(out)))
			} else if (!hc) {
				out <- rowColCombine(x, row, 1:ncol(x))
				if (hl) {
					out <- cbind(out[,1], out[,2], rep(lyr, each=nrow(out)))
				} else {
					out <- cbind(out[,1], out[,2], rep(1:nlyr(x), each=nrow(out)))
				}
			} else if (!hr) {
				out <- rowColCombine(x, 1:nrow(x), col)
				if (hl) {
					out <- cbind(out[,1], out[,2], rep(lyr, each=nrow(out)))
				} else {
					out <- cbind(out[,1], out[,2], rep(1:nlyr(x), each=nrow(out)))
				}
			} else if (!hl) {
				out <- rowColCombine(x, row, col)
				out <- cbind(out[,1], out[,2], rep(1:nlyr(x), each=nrow(out)))
			}
			colnames(out) <- c("row", "col", "lyr")
			out
		}
	}
)



