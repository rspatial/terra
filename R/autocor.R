# Author: Robert J. Hijmans
# Date : April 2011 / Jan 2021
# Version 1.0
# Licence GPL v3


.checkngb <- function(ngb, mustBeOdd=FALSE) {
	ngb <- as.integer(round(ngb))
	if (length(ngb) == 1) {
		ngb <- c(ngb, ngb)
	} else if (length(ngb) > 2) {
		stop('ngb should be a single value or two values')
	}
	if (min(ngb) < 1) { stop("ngb should be larger than 1") } 
	if (mustBeOdd) {
		if (any(ngb %% 2 == 0)) {
			stop('neighborhood size must be an odd number')
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
		method <- match.arg(tolower(method), c("moran", "geary", "gi", "gi*"))
		d <- dim(w)
		n <- length(x)
		if ((d[1] != d[2]) || (d[1] != n)) {
			stop("w must be a square matrix with sides the size of x")
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
			diag(w) <- 0
			Gi <- colSums(x * w) / (sum(x)-x)
			Ei <- rowSums(w) / (n-1) 

			# variance following spdep::localG
			xibar <- (sum(x) - x)/(n - 1)
		    si2 <- (sum(x^2) - x^2)/(n - 1) - xibar^2
			VG <- si2 * (((n - 1) * rowSums(w^2) - rowSums(w)^2)/(n - 2))
			VG <- VG/((sum(x) - x)^2)

			(Gi-Ei)/sqrt(VG)

		} else if (method == "gi*") {
			if (any(diag(w) == 0)) {
				warn("autocor", "it is unexpected that a weight matrix for Gi* has diagonal values that are zero")
			}
			Gi <- colSums(x * w) / sum(x)
			Ei <- rowSums(w) / n
			# variance following spdep::localG
			si2 <- sum(scale(x, scale = FALSE)^2)/n
			VG <- (si2 * ((n * rowSums(w^2) - rowSums(w)^2)/(n - 1))) / (sum(x)^2)

			(Gi-Ei)/sqrt(VG)

			
		} else {
			stop()
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
				wZiZj <- unlist(global(wZiZj, "sum"))
				z2 <- unlist(global(z*z, "sum"))
				n <- ncell(z) - unlist(global(is.na(z), "sum", na.rm=TRUE))
				zz <- ifel(is.na(x), NA, 1)
				W <- focal( zz, w=w, fun="sum") 
				NS0 <- n / unlist(global(W, "sum", na.rm=TRUE))
				NS0 * wZiZj / z2
			} else { # geary
				w <- .getFilter(w, warn=FALSE)
				i <- trunc(length(w)/2)+1 
				n <- ncell(x) - unlist(global(is.na(x), "sum"))
				fun <- function(x,...) sum((x-x[i])^2, ...)
				f <- focal(x, w=dim(w), fun=fun, na.rm=TRUE)
				Eij <- unlist(global(f, "sum"))
				xx <- ifel(is.na(x), NA ,1)
				W <- focal(xx, w=w, na.rm=TRUE ) 
				z <- 2 * unlist(global(W, "sum")) * unlist(global((x - unlist(global(x, "mean")))^2, "sum"))
				(n-1)*Eij/z
			}
		} else { # local
			if (method == "moran") {	
				z <- x - unlist(global(x, "mean", na.rm=TRUE))
				zz <- ifel(is.na(x), NA, 1)
				W  <- focal(zz, w=w, na.rm=TRUE)		
				lz <- focal(z, w=w, na.rm=TRUE) / W
					
				n <- ncell(x) - unlist(global(is.na(x), "sum", na.rm=TRUE))
				s2 <- unlist(global(x, "sd")^2 )
				(z / s2) * lz	
			} else {
				w <- .getFilter(w)
				i <- trunc(length(w)/2)+1 
				fun <- function(x,...) sum((x-x[i])^2, ...)
				Eij <- focal(x, w=dim(w), fun=fun, na.rm=TRUE)
				s2 <- unlist(global(x, "sd"))^2 
				n <- ncell(x) - unlist(global(is.na(x), "sum"))	
				Eij / s2
			}
		}
	}
)



