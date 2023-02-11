# Author: Robert J. Hijmans
# Date : April 2011 / Jan 2021
# Version 1.0
# License GPL v3


.checkngb <- function(ngb, mustBeOdd=FALSE) {
	ngb <- as.integer(round(ngb))
	if (length(ngb) == 1) {
		ngb <- c(ngb, ngb)
	} else if (length(ngb) > 2) {
		error("autocor", "ngb should be a single value or two values")
	}
	if (min(ngb) < 1) { stop("ngb should be larger than 1") }
	if (mustBeOdd) {
		if (any(ngb %% 2 == 0)) {
			error("autocor", "neighborhood size must be an odd number")
		}
	}
	return(ngb)
}


.getFilter <- function(w, warn=TRUE) {
	if (!is.matrix(w)) {
		w <- .checkngb(w)
		w <- matrix(1, nrow=w[1], ncol=(w[2]))
		w[ceiling(dim(w)[1]/2), ceiling(dim(w)[2]/2)] <- 0
	} else {
		if (w[ceiling(dim(w)[1]/2), ceiling(dim(w)[2]/2)] != 0) {
			if (warn) {
				warning('central cell of weights matrix (filter) was set to zero')
			}
			w[ceiling(dim(w)[1]/2), ceiling(dim(w)[2]/2)] <- 0
		}
		stopifnot(all(w >= 0))
	}
	if (min(dim(w) %% 2)==0) {
		stop('dimensions of weights matrix (filter) must be uneven')
	}
	w
}



setMethod("autocor", signature(x="numeric"),
	function(x, w, method="moran") {
		method <- match.arg(tolower(method), c("moran", "geary", "gi", "gi*", "mean", "locmor"))

		if (all(is.na(x))) {
			error("autocor", "all values are NA")
		}
		if (any(is.na(w))) {
			error("autocor", "NA value(s) in the weight matrix")
		}

		d <- dim(w)
		n <- length(x)
		if ((d[1] != d[2]) || (d[1] != n)) {
			error("autocor", "w must be a square matrix with sides the size of x")
		}

		if (method %in% c("moran", "geary", "locmor", "gi")) {
			if (any(as.numeric(diag(w)) != 0)) {
				warn("autocor", paste("it is unexpected that a weight matrix for", method, "has diagonal values that are not zero"))
			}
		} else if (method %in% c("gi*")) {
			if (any(as.numeric(diag(w)) == 0)) {
				warn("autocor", paste("it is unexpected that a weight matrix for", method, "has diagonal values that are zero"))
			}
		}


		if (method == "moran") {
			dx <- x - mean(x, na.rm=TRUE)
			pm <- matrix(rep(dx, each=n) * dx, ncol=n)
			(n / sum(dx^2)) * sum(pm * w) / sum(w)
		} else if (method == "geary") { # geary
			dx <- x - mean(x, na.rm=TRUE)
			pm <- matrix(rep(dx, each=n) - dx, ncol=n)^2
			((n-1)/sum((dx)^2)) * sum(w * pm) / (2 * sum(w))
		} else if (method == "gi") {
			if (any(as.numeric(diag(w)) != 0)) {
				warn("autocor", "it is unexpected that a weight matrix for Gi has diagonal values that are not zero")
			}
			diag(w) <- 0
			sumxminx <- sum(x, na.rm=TRUE) - x
			Gi <- colSums(x * w) / sumxminx
			Ei <- rowSums(w) / (n-1)

			# variance following spdep::localG
			xibar <- sumxminx/(n - 1)
		    si2 <- (sum(x^2) - x^2)/(n - 1) - xibar^2
			VG <- si2 * (((n - 1) * rowSums(w^2) - rowSums(w)^2)/(n - 2))
			VG <- VG/((sumxminx)^2)

			(Gi-Ei)/sqrt(VG)

		} else if (method == "gi*") {

			if (any(as.numeric(diag(w)) == 0)) {
				warn("autocor", "it is unexpected that a weight matrix for Gi* has diagonal values that are zero")
			}
			Gi <- colSums(x * w) / sum(x)
			Ei <- rowSums(w) / n
			# variance following spdep::localG
			si2 <- sum(scale(x, scale = FALSE)^2)/n
			VG <- (si2 * ((n * rowSums(w^2) - rowSums(w)^2)/(n - 1))) / (sum(x)^2)

			(Gi-Ei)/sqrt(VG)

		} else if (method == "locmor") {
			if (any(as.numeric(diag(w)) != 0)) {
				warn("autocor", "it is unexpected that a weight matrix for local Moran has diagonal values that are not zero")
			}
			z <- x - mean(x, na.rm=TRUE)
			mp <- z / ( (sum(z^2, na.rm=TRUE) / length(x)) )
			mp * apply(w, 1, function(i) {
					sum(z * i, na.rm=TRUE)
				} )

		} else if (method == "mean") {
			j <- is.na(x)
			x[j] <- 0
			w[j,j] <- 0
			m <- apply(w, 1, function(i) {
					sum(x * i) / sum(i)}
				)
			m[j] <- NA
			m
		}
 	}
)



setMethod("autocor", signature(x="SpatRaster"),
	function(x, w=matrix(c(1,1,1,1,0,1,1,1,1),3), method="moran", global=TRUE) {

		method <- match.arg(tolower(method), c("moran", "geary"))

		if (nlyr(x) > 1) {
			warn("autocor", "only the first layer of x is used")
			x <- x[[1]]
		}

		if (global) {
			if (method == "moran") {
				z <- x - unlist(global(x, "mean", na.rm=TRUE))
				wZiZj <- focal(z, w=w, fun="sum", na.rm=TRUE)
				wZiZj <- wZiZj * z
				wZiZj <- unlist(global(wZiZj, "sum", na.rm=TRUE))
				z2 <- unlist(global(z*z, "sum", na.rm=TRUE))
				n <- ncell(z) - unlist(global(is.na(z), "sum"))
				zz <- ifel(is.na(x), NA, 1)
				W <- focal( zz, w=w, fun="sum")
				NS0 <- n / unlist(global(W, "sum", na.rm=TRUE))
				m <- NS0 * wZiZj / z2
				names(m) <- names(x)
				m
			} else { # geary
				w <- .getFilter(w, warn=FALSE)
				i <- trunc(length(w)/2)+1
				n <- ncell(x) - unlist(global(is.na(x), "sum"))
				fun <- function(x,...) sum((x-x[i])^2, ...)
				f <- focal(x, w=dim(w), fun=fun, na.rm=TRUE)
				Eij <- unlist(global(f, "sum", na.rm=TRUE))
				xx <- ifel(is.na(x), NA ,1)
				W <- focal(xx, w=w, na.rm=TRUE )
				z <- 2 * unlist(global(W, "sum", na.rm=TRUE)) *
					unlist(global((x - unlist(global(x, "mean", na.rm=TRUE)))^2, "sum", na.rm=TRUE))
				g <- (n-1)*Eij/z
				names(g) <- names(x)
				g
			}
		} else { # local
			if (method == "moran") {
				z <- x - unlist(global(x, "mean", na.rm=TRUE))
				zz <- ifel(is.na(x), NA, 1)
				W  <- focal(zz, w=w, na.rm=TRUE)
				lz <- focal(z, w=w, na.rm=TRUE) / W

				##n <- ncell(x) - unlist(global(is.na(x), "sum"))
				s2 <- unlist(global(x, "sd", na.rm=TRUE)^2)
				m <- (z / s2) * lz
				names(m) <- names(x)
				m
			} else {
				w <- .getFilter(w)
				i <- trunc(length(w)/2)+1
				fun <- function(x,...) sum((x-x[i])^2, ...)
				Eij <- focal(x, w=dim(w), fun=fun, na.rm=TRUE)
				s2 <- unlist(global(x, "sd", na.rm=TRUE))^2
				##n <- ncell(x) - unlist(global(is.na(x), "sum"))
				g <- Eij / s2
				names(g) <- names(x)
				g
			}
		}
	}
)



