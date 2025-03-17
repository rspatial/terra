

hexcols <- function(out) {

	get_col <- function(cols, alpha) {
		if (isTRUE(alpha < 255)) {
			grDevices::rgb(t(grDevices::col2rgb(cols, alpha=TRUE)), alpha=alpha, maxColorValue=255)
		} else {
			i <- !grepl("^#", cols)
			cols[i] <- grDevices::rgb(t(grDevices::col2rgb(cols[i], alpha=FALSE)), maxColorValue=255)	
			cols
		}
	}

	if (NCOL(out$cols) == 1) {
		out$cols <- get_col(out$cols, out$alpha)
	} else if (NCOL(out$cols) == 2) {
		out$cols[,2] <- get_col(out$cols[,2], out$alpha)
	} else if (NCOL(out$cols) == 3) {
		out$cols[,3] <- get_col(out$cols[,3], out$alpha)
	}
	
	out

}


.default.pal <- function() {
	opt.pal <- options("terra.pal")[[1]]
	if (is.null(opt.pal))  {
		map.pal("viridis", 100)
	} else {
		opt.pal
	}
}

.get_nrnc <- function(nr, nc, nl) {
	if (missing(nc)) {
		nc <- ceiling(sqrt(nl))
	} else {
		nc <- max(1, min(nl, round(nc)))
	}
	if (missing(nr)) {
		nr <- ceiling(nl / nc)
	} else {
		nr <- max(1, min(nl, round(nr)))
		nc <- ceiling(nl / nr)
	}
	c(nr, nc)
}



.get_breaks <- function(x, n, method, r=NULL) {
	#x <- x[!is.na(x)]
	
	if (is.function(method)) {
		if (!is.null(r)) {
			if (!is.na(r[1])) { 
				x[ x < r[1] ] <- NA
			} 
			if (!is.na(r[2])) { 
				x[ x > r[2] ] <- NA
			} 
		}
		breaks <- method(x)
	} else if (method[1]=="cases") {
		if (!is.null(r)) {
			if (!is.na(r[1])) { 
				x[ x < r[1] ] <- NA
			} 
			if (!is.na(r[2])) { 
				x[ x > r[2] ] <- NA
			} 
		}
		n <- n+1
		i <- seq(0, 1, length.out=n)
		breaks <- quantile(x, i, na.rm=TRUE)
		breaks <- unique(breaks)
		if ((breaks[1] %% 1) != 0) {
			breaks[1] <- breaks[1] - 0.000001
		}
		if ((breaks[n] %% 1) != 0) {
			breaks[n] <- breaks[n] + 0.000001
		}
	} else { # if (method=="eqint") {
		if (is.null(r)) {
			r <- c(min(x, na.rm=TRUE), max(x, na.rm=TRUE))
		} else if (any(is.na(r))) {
			if (is.na(r[1])) r[1] <- min(x, na.rm=TRUE)
			if (is.na(r[2])) r[2] <- max(x, na.rm=TRUE)
		}
		small <- 1e-16
		if ((r[1] %% 1) != 0) { r[1] <- r[1] - small }
		if ((r[2] %% 1) != 0) { r[2] <- r[2] + small }
		breaks <- seq(r[1] , r[2], length.out=n+1)
	}
	unique(breaks)
}



