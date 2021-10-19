
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
	px <- py <- c()
	hasDF <- ncol(x) > 0
	for (i in 1:nrow(xx)) {
		for (j in 1:nrow(yy)) {
			if (hasDF) {
				if (all(as.data.frame(xx[i,]) != as.data.frame(yy[j,]))) {
					next
				}
			}
			if (relate(xx[i,], yy[j,], "touches")) {
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

	upx <- unique(px)
	vvv <- list()
	hasDF <- ncol(x) > 0
	
	for (i in 1:length(upx)) {
		pyi <- py[px==upx[i]]
		vv <- x[px[i],]
		for (j in 1:length(pyi)) {

			vv <- aggregate(rbind(vv, x[pyi[j],]), dissolve=FALSE)
		}
		vvv[[i]] <- vv
	}

	out <- x[-c(px, py), ]
	out <- c(vvv, out)
	do.call(rbind, out)
}

