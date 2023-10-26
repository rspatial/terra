
get_groups <- function(x, y) {
	j <- 1
	outx <- outy <- list()

	for (i in 1:length(x)) {
		if (is.na(x[i])) next
		gx <- stats::na.omit(x[x[i] == x] )
		gy <- y[x %in% gx]
		nx <- ny <- 0
		while(TRUE) {
			if (nx == length(gx)) break
			ny <- length(gy)
			nx <- length(gx)
			if ((ny == length(y) || (nx == length(x)))) break
			ux <- unique( x[y %in% gy] )
			gy <- y[x %in% ux]
			gx <- x[y %in% gy]
		}
		x[x %in% gx] <- NA
		y[y %in% gy] <- NA
		outx[[j]] <- gx
		outy[[j]] <- gy
		j <- j + 1
	}
	list(outx, outy)
}


connect_dateline <- function(x) {
	east <- west <- c()
	for (i in 1:nrow(x)) {
		e <- ext(x[i,])
		if (xmin(e) <= -180) {
			west <- c(west, i)
		} else if (xmax(e) >= 180) {
			east <- c(east, i)
		}
	}
	if ((length(east) == 0) || (length(west) == 0)) {
		return(x)
	}

	xx <- shift(x[west,], 360, 0)
	yy <- x[east, ]
	#zz <- x[-c(east, west), ]
	px <- py <- c()
	hasDF <- ncol(x) > 0
	for (i in 1:nrow(xx)) {
		for (j in 1:nrow(yy)) {
			if (hasDF) {
				if (all(as.data.frame(xx[i,]) != as.data.frame(yy[j,]))) {
					next
				}
			}
			if (is.related(xx[i,], yy[j,], "touches")) {
				px <- c(px, i)
				py <- c(py, j)
			}
		}
	}
	if ((length(px) == 0)) {
		return(x)
	}

	px <- west[px]
	py <- east[py]

	groups <- get_groups(px, py)
	xg <- groups[[1]]
	yg <- groups[[2]]
	vvv <- list()
	for (i in 1:length(xg)) {
		vvv[[i]] <- aggregate(x[unique(c(xg[[i]], yg[[i]])), ], dissolve=TRUE)
	}
	out <- x[-(unique(unlist(groups))), ]
	out <- c(vvv, out)
	do.call(rbind, out)
}


