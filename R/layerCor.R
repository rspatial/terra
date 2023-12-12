# Robert Hijmans
# Date : Nov 2021
# Version 1.0
# Licence GPL v3

# Computation of the weighted covariance and (optionally) weighted means of bands in an Raster.
# based on code for "raster" by Jonathan Greenberg and Robert Hijmans
# partly based on code by Mort Canty



setMethod("layerCor", signature(x="SpatRaster"),
	function(x, fun, w, asSample=TRUE, use="everything", maxcell=Inf, ...) {

		ops <- c("everything", "complete.obs", "pairwise.complete.obs", "masked.complete")

		# backwards compatibility 
		na.rm <- list(...)$na.rm
		if (isTRUE(na.rm) && (use == "everything")) {
			use <- "pairwise.complete.obs"
		} 

		use <- match.arg(use, ops)
		if (use != "everything") {
			na.rm <- TRUE
		} else {
			na.rm <- FALSE			
		}

		stopifnot(is.logical(asSample) & !is.na(asSample))
		
		nl <- nlyr(x)
		if (nl < 2) {
			error("layerCor", "x must have at least 2 layers")
		}
		
		if (inherits(fun, "character")) {
			fun <- tolower(fun)
			if (!(fun %in% c("cov", "weighted.cov", "pearson", "cor"))) {
				error("layerCor", "character function names must be one of: 'cor', 'cov', or 'weighted.cov'")
			}
			# backwards compatibility
			if (fun == "pearson") fun = "cor"
			if (fun == "weighted.cov") {
				if (missing(w))	{
					error("layerCor", "to compute weighted covariance a weights layer should be provided")
				}
				stopifnot( nlyr(w) == 1 )
				x <- c(w, x)
			}
		} else {
			FUN <- fun
			fun <- ""
		}
		
		if (maxcell < ncell(x)) {
			x <- spatSample(x, size=maxcell, "regular", as.raster=TRUE)
		}
		n <- ncell(x)

		# for "cor" masking is done in cpp code
		if ((use == "complete.obs") && (fun != "cor")) { 
			x <- mask(x, anyNA(x), maskvalue=TRUE)
		}
		
		if (fun == "weighted.cov") {
			w <- x[[1]]
			x <- x[[-1]]

			means <- mat <- matrix(NA, nrow=nl, ncol=nl)
			colnames(means) <- rownames(means) <- colnames(mat) <- rownames(mat) <- names(x)
			sqrtw <- sqrt(w)
			for(i in 1:nl) {
				for(j in i:nl) {
					s <- c(x[[c(i,j)]])
					if (use == "pairwise.complete.obs") {
						s <- mask(c(s, w), anyNA(s), maskvalue=TRUE)
						ww <- s[[3]]
						s <- s[[1:2]]
						sumw <- unlist(global(ww, fun="sum", na.rm=TRUE) )
						avg <- unlist(global(s * ww, fun="sum", na.rm=TRUE)) / sumw
					} else {
						sumw <- unlist(global(w, fun="sum", na.rm=na.rm) )
						avg <- unlist(global(s * w, fun="sum", na.rm=na.rm)) / sumw
					}					
					sumw <- sumw - asSample
					s <- prod( (s - avg) * sqrtw )
					v <- unlist(global(s, fun="sum", na.rm=na.rm)) / sumw
					mat[j,i] <- mat[i,j] <- v
					means[i,j] <- avg[1]
					means[j,i] <- avg[2]
				}
			}
			return( list(weighted_covariance=mat, weighted_mean=means) )
		} else if (fun == "cov") {
			means <- mat <- nn <- matrix(NA, nrow=nl, ncol=nl)
			colnames(means) <- rownames(means) <- colnames(mat) <- rownames(mat) <- names(x)
			n_ij <- n
			for(i in 1:nl) {
				for(j in i:nl) {
					s <- x[[c(i,j)]]
					if (use == "pairwise.complete.obs") {
						m <- anyNA(s)
						s <- mask(s, m, maskvalue=TRUE)
						n_ij <- n - global(m, fun="sum")$sum
					}
					avg <- unlist(global(s, fun="mean", na.rm=na.rm) )
					r <- prod(s - avg)
					v <- unlist(global(r, fun="sum", na.rm=na.rm)) / (n_ij - asSample)
					mat[j,i] <- mat[i,j] <- v
					means[i,j] <- avg[1]
					means[j,i] <- avg[2]
					nn[i,j] <- nn[j,i] <- n_ij
				}
			}

			return( list(covariance=mat, mean=means, n=nn) )

		} else if (fun == "cor") {
			if (isTRUE(list(...)$old)) {
				old_pearson(x, asSample=asSample, na.rm=na.rm, nl=nl, n=n, mat=mat)
			} else {			
				opt <- spatOptions()
				m <- x@cpp$layerCor("cor", use, asSample, opt)
				x <- messages(x)
				m <- lapply(m, function(i) {
					matrix(i, nrow=nl, byrow=TRUE, dimnames=list(names(x), names(x)))
				})
				names(m) <- c("correlation", "mean", "n")
				return(m)
			}
		} else {
			v <- spatSample(x, size=maxcell, "regular", na.rm=na.rm, warn=FALSE)
			if (use %in% c("complete.obs", "complete.masked")) {
				v <- na.omit(v)
			}
			for(i in 1:nl) {
				for(j in i:nl) {
					if (use == "pairwise.complete.obs") {
						vij <- na.omit(v[,c(i,j)])
						mat[j,i] <- mat[i,j] <- FUN(vij[,1], vij[,2], ...)
					} else {
						mat[j,i] <- mat[i,j] <- FUN(v[,i], v[,j], ...)
					}
				}
			}
			mat
		}
	}
)




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
	
