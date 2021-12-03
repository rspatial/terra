# Robert Hijmans
# Date : Nov 2021
# Version 1.0
# Licence GPL v3

# Computation of the weighted covariance and (optionally) weighted means of bands in an Raster.

# based on code for "raster" by Jonathan Greenberg and Robert Hijmans partly based on code by Mort Canty

if (!isGeneric("layerCor")) {setGeneric("layerCor", function(x, ...) standardGeneric("layerCor"))}

setMethod("layerCor", signature(x="SpatRaster"), 
	function(x, fun, w, asSample=TRUE, na.rm=FALSE, maxcell=Inf, ...) {

		stopifnot(is.logical(asSample) & !is.na(asSample))
		nl <- nlyr(x)
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
		
		if (fun == "weighted.cov") {
			if (missing(w))	{
				stop("to compute weighted covariance a weights layer should be provided")
			}
			stopifnot( nlyr(w) == 1 )

			if (na.rm) {
			# a cell is set to NA if it is NA in any layer. That is not ideal, but easier and quicker
				nas <- sum(x) * w
				x <- mask(x, nas)
				w <- mask(w, nas)
			}

			sumw <- unlist(global(w, fun="sum", na.rm=na.rm) )
			means <- unlist(global(x * w, fun="sum", na.rm=na.rm)) / sumw
			sumw <- sumw - asSample
			
			x <- (x - means) * sqrt(w)
			

			for(i in 1:nl) {
				for(j in i:nl) {
					r <- x[[i]] * x[[j]]
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
					r <- x[[i]] * x[[j]]
					if (na.rm) {
						v <- unlist(global(r, fun="sum", na.rm=na.rm)) / (n - unlist(global(is.na(r), fun="sum")) - asSample)
					} else {
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

			means <- unlist(global(x, fun="mean", na.rm=na.rm) )
			sds <- unlist(global(x, fun="sd", na.rm=na.rm) )
			x <- (x - means)
			
			for(i in 1:nl) {
				for(j in i:nl) {
					r <- x[[i]] * x[[j]]
					v <- unlist(global(r, fun="sum", na.rm=na.rm))
					if (na.rm) {
						nas <- unlist(global(is.na(r), fun="sum"))
						v <- v / ((n - nas - asSample) * sds[i] * sds[j])
					} else {
						v <- v / ((n - asSample) * sds[i] * sds[j])
					}
					mat[j,i] <- mat[i,j] <- v
				}
			}
			names(means) <- names(x)
			covar <- list(mat, means)
			names(covar) <- c("pearson", "mean")
			return(covar)
			
		} else {

			v <- spatSample(x, size=maxcell, "regular", na.rm=na.rm)
			for(i in 1:nl) {
				for(j in i:nl) {
					mat[j,i] <- mat[i,j] <- FUN(v[,i], v[,j])
				}
			}
			mat
			
		}
	}
)

