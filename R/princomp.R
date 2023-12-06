

setMethod("princomp", signature(x="SpatRaster"),
	function(x, cor=FALSE, fix_sign=TRUE) {
		stopifnot(nlyr(x) > 1)
		xcov <- layerCor(x, fun="cov", na.rm=TRUE, asSample=FALSE)
		model <- princomp(covmat=xcov$covariance, cor=cor, fix_sign=fix_sign)
		model$center <- diag(xcov$mean)
		n <- diag(xcov$n)
		if (cor) {
		## scale as population sd
			S <- diag(xcov$covariance)
			model$scale <- sqrt(S)
		}
		model$n.obs <- as.dist(xcov$n)
		model
	}
)

