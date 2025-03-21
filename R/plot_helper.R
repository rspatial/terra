

# to be merged with the vector variant.
.generic.interval <- function(out, Z) {
	if (is.null(out$breaks)) {
		out$breaks <- 5
	}
	if (length(out$breaks) == 1) {
		out$breaks <- .get_breaks(Z, out$breaks, out$breakby, out$range)
	}

	if (!is.null(out$leg$digits)) {
#		out$leg$legend <- substr(formatC(levs, digits=digits, format = "f", flag="#"), 1, digits+1)
		fz <- cut(as.numeric(Z), out$breaks, include.lowest=TRUE, right=FALSE, dig.lab=out$leg$digits)
	} else {
		fz <- cut(as.numeric(Z), out$breaks, include.lowest=TRUE, right=FALSE)
	}


	out$vcut <- as.integer(fz)
	levs <- levels(fz)
	nlevs <- length(levs)

	cols <- out$cols
	ncols <- length(cols)
	if (nlevs < ncols) {
		i <- trunc((ncols / nlevs) * 1:nlevs)
		cols <- cols[i]
	} else {
		cols <- rep_len(cols, nlevs)
	}
	
	#out$cols <- cols
	out$leg$fill <- cols
	#out$leg$levels <- levels(fz)
	if (!is.null(out$leg$legend)) {
		stopifnot(length(out$leg$legend) == nlevs)
	} else {
		levs <- gsub("]", "", gsub(")", "", gsub("\\[", "", levs)))
		levs <- paste(levs, collapse=",")
		m <- matrix(as.numeric(unlist(strsplit(levs, ","))), ncol=2, byrow=TRUE)
		if (!is.null(out$leg$digits)) {
			m <- prettyNumbs(m, out$leg$digits)
		}
		m <- apply(m, 1, function(i) paste(i, collapse=" - "))
		m <- gsub("-Inf -", "<=", m)
		i <- grep("- Inf", m)
		if (length(i) == 1) {
			m[i] <- gsub("- Inf", "", m[i])
			m[i] <- paste(">", m[i])				
		}	
		out$leg$legend <- m
	}
	out$leg$digits <- NULL
	out
}


.get_nudge <- function(a) {
	if (is.null(a)) {
		a <- rep(0, 4)
	} else if (length(a) == 0) {
		a <- rep(0, 4)
	} else if (length(a) == 1) {
		a <- c(a, a, 0, 0)
	} else if (length(a) == 2) {
		a <- c(a[1], a[1], a[2], a[2])		
	} else if (length(a) == 3) {
		a <- c(a[1], a[2], a[3], a[3])		
	} 
	a
}

.nudge_ext <- function(x, e) {
	a <- .get_nudge(x$leg[["nudge"]]) 
	e <- x$leg[["ext"]]
	e$xmin <- e$xmin + a[1]
	e$xmax <- e$xmax + a[2]
	e$ymin <- e$ymin + a[3]
	e$ymax <- e$ymax + a[4]
	e$dy <- e$dy + a[4]-a[3]
	e$dx <- e$dx + a[2]-a[1]
	x$leg$ext <- e
	x
}

.nudge_xy <- function(xy, a) {
	a <- .get_nudge(a) 
	xy[1] <- xy[1] + a[1]
	xy[2] <- xy[2] + a[3]
	xy
}



prettyNumbs <- function(x, digits) {
	x <- formatC(x, digits=digits, format = "f", flag="#")
	x <- substr(x, 1, digits+1)
	gsub("\\.$", "", x)
}


add_more <- function(fun, i) {
	if (!is.null(fun) && is.function(fun)) {
		if (!is.null(formals(fun))) {
			fun(i)
		} else {
			fun()
		}
	}
}


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



