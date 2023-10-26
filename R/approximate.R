# Author: Robert J. Hijmans
# Date : February 2012
# Version 1.0
# Licence GPL v3


setMethod("approximate", signature(x="SpatRaster"),
function(x, method="linear", yleft, yright, rule=1, f=0, ties=mean, z=NULL, NArule=1, filename="", ...) {

	out <- rast(x, keeptime=TRUE)
	nl <- nlyr(out)
	if (nl < 2) {
		warning('cannot interpolate with a single layer')
		return(x)
	}

	if (is.null(z)) {
		xout <- time(x)
		if (any(is.na(xout))) {
			xout <- 1:nl
		}
	} else {
		if (length(z)!= nl) {
			error("approximate", "length of z does not match nlyr(x)")
		}
		xout <- z
	}

	ifelse((missing(yleft) & missing(yright)), ylr <- 0L, ifelse(missing(yleft), ylr <- 1L, ifelse(missing(yright), ylr <- 2L, ylr <- 3L)))

    nc <- ncol(out)

	readStart(x)
	on.exit(readStop(x))
	b <- writeStart(out, filename, sources=sources(x), ...)
	for (i in 1:b$n) {
		v <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
		s <- .rowSums(is.na(v), nrow(v), nl)

		if (isTRUE(NArule==1)) {
			j <- s == (nl-1) # one non-NA only
			if (length(j) > 0 ) {
				v[j, ] <- apply(v[j,,drop=FALSE ], 1, max, na.rm=TRUE)
			}
		}
		j <- (s < nl-1) # need at least two
		if (length(j) > 0 ) {
			if (ylr==0) {
				v[j,] <- t( apply(v[j,,drop=FALSE], 1, function(x) stats::approx(x=xout, y=x, xout=xout, method=method, rule=rule, f=f, ties=ties)$y ) )
			} else if (ylr==1) {
				v[j,] <- t( apply(v[j,,drop=FALSE], 1, function(x) stats::approx(x=xout, y=x, xout=xout, method=method, yright=yright, rule=rule, f=f, ties=ties)$y ) )
			} else if (ylr==2) {
				v[j,] <- t( apply(v[j,,drop=FALSE], 1, function(x) stats::approx(x=xout, y=x, xout=xout, method=method, yleft=yleft, rule=rule, f=f, ties=ties)$y ) )
			} else {
				v[j,] <- t( apply(v[j,,drop=FALSE], 1, function(x) stats::approx(x=xout, y=x, xout=xout, method=method, yright=yright, yleft=yleft, rule=rule, f=f, ties=ties)$y ) )
			}
		}
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	writeStop(out)
}
)

