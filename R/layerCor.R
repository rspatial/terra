# Robert Hijmans
# Date : Nov 2021
# Version 1.0
# Licence GPL v3

# Computation of the weighted covariance and (optionally) weighted means of bands in an Raster.

# based on code for "raster" by Jonathan Greenberg and Robert Hijmans partly based on code by Mort Canty


old_pearson <- function(x, asSample, na.rm, nl, n, mat) {
	if (na.rm) {
		means <- matrix(NA, nrow=2, ncol=nlyr(x))
		for(i in 1:(nl-1)) {
			for(j in (i+1):nl) {
				m <- anyNA(x[[c(i,j)]])
				a <- mask(x[[i]], m, maskvalue=TRUE)
				b <- mask(x[[j]], m, maskvalue=TRUE)
				xx <- c(a, b)
				mns <- unlist(global(xx, fun="mean", na.rm=na.rm) )
				means[2,i] <- mns[1]
				means[1,j] <- mns[2]						
				sds <- unlist(global(xx, fun="sd", na.rm=na.rm) )
				r <- prod(xx - mns)
				nas <- unlist(global(is.na(r), fun="sum"))
				v <- unlist(global(r, fun="sum", na.rm=na.rm))
				v <- v / ((n - nas - asSample) * sds[1] * sds[2])
				mat[j,i] <- mat[i,j] <- v
			}
		}
		colnames(means) <- names(x)
	} else {
		means <- unlist(global(x, fun="mean", na.rm=na.rm) )
		sds <- unlist(global(x, fun="sd", na.rm=na.rm) )
		x <- (x - means)
		
		for(i in 1:(nl-1)) {
			for(j in i:nl) {
				r <- x[[i]] * x[[j]]
				v <- unlist(global(r, fun="sum", na.rm=na.rm))
				v <- v / ((n - asSample) * sds[i] * sds[j])
				mat[j,i] <- mat[i,j] <- v
			}
		}
		means <- matrix(means, nrow=1)
		colnames(means) <- names(x)
	}
	diag(mat) <- 1
	covar <- list(mat, means)
	names(covar) <- c("pearson", "mean")
	return(covar)
}
	


if (!isGeneric("layerCor")) {setGeneric("layerCor", function(x, ...) standardGeneric("layerCor"))}

setMethod("layerCor", signature(x="SpatRaster"),
	function(x, fun, w, asSample=TRUE, na.rm=FALSE, maxcell=Inf, ...) {

		stopifnot(is.logical(asSample) & !is.na(asSample))
		nl <- nlyr(x)
		if (nl < 2) {
			error("layerCor", "x must have at least 2 layers")
		}
		
		n <- ncell(x)
		mat <- matrix(NA, nrow=nl, ncol=nl)
		colnames(mat) <- rownames(mat) <- names(x)


		if (inherits(fun, "character")) {
			fun <- tolower(fun)
			stopifnot(fun %in% c("cov", "weighted.cov", "pearson"))
		} else {
			FUN <- fun
			fun <- ""
		}

		if (maxcell < Inf) {
			x <- spatSample(x, size=maxcell, "regular", as.raster=TRUE)
		}
		
		if (fun == "weighted.cov") {
			if (missing(w))	{
				stop("to compute weighted covariance a weights layer should be provided")
			}
			stopifnot( nlyr(w) == 1 )

			sumw <- unlist(global(w, fun="sum", na.rm=na.rm) )
			means <- unlist(global(x * w, fun="sum", na.rm=na.rm)) / sumw
			sumw <- sumw - asSample
			x <- (x - means) * sqrt(w)

			for(i in 1:nl) {
				for(j in i:nl) {
					if (na.rm) {
						m <- anyNA(x[[c(i,j)]])
						a <- mask(x[[i]], m, maskvalue=TRUE)
						b <- mask(x[[j]], m, maskvalue=TRUE)
						r <- a * b
					} else {
						r <- x[[i]] * x[[j]]
					}
					v <- unlist(global(r, fun="sum", na.rm=na.rm)) / sumw
					mat[j,i] <- mat[i,j] <- v
				}
			}
			names(means) <- names(x)
			cov.w <- list(mat, means)
			names(cov.w) <- c("weighted_covariance", "weighted_mean")
			return(cov.w)

		} else if (fun == "cov") {

			means <- unlist(global(x, fun="mean", na.rm=na.rm) )
			x <- (x - means)

			for(i in 1:nl) {
				for(j in i:nl) {
					if (na.rm) {
						m <- anyNA(x[[c(i,j)]])
						a <- mask(x[[i]], m, maskvalue=TRUE)
						b <- mask(x[[j]], m, maskvalue=TRUE)
						r <- x[[i]] * x[[j]]
						v <- unlist(global(r, fun="sum", na.rm=na.rm)) / (n - unlist(global(r, fun="isNA")) - asSample)
					} else {
						r <- x[[i]] * x[[j]]
						v <- unlist(global(r, fun="sum", na.rm=na.rm)) / (n - asSample)
					}
					mat[j,i] <- mat[i,j] <- v
				}
			}

			names(means) <- names(x)
			covar <- list(mat, means)
			names(covar) <- c("covariance", "mean")
			return(covar)

		} else if (fun == "pearson") {
			if (isTRUE(list(...)$old)) {
				old_pearson(x, asSample=asSample, na.rm=na.rm, nl=nl, n=n, mat=mat)
			} else {
				opt <- spatOptions()
				m <- x@ptr$layerCor("pearson", na.rm, asSample, opt)
				x <- messages(x)
				mat <- matrix(m[[1]], nrow=nl, byrow=TRUE)
				means <- matrix(m[[2]], nrow=nl, byrow=TRUE)
				return( list(pearson=mat, mean=means) )
			}
		} else {

			v <- spatSample(x, size=maxcell, "regular", na.rm=na.rm, warn=FALSE)
			for(i in 1:nl) {
				for(j in i:nl) {
					mat[j,i] <- mat[i,j] <- FUN(v[,i], v[,j], ...)
				}
			}
			mat
		}
	}
)

