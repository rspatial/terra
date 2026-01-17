# Author: Robert J. Hijmans
# Date : March 2009
# Version 1.0
# License GPL v3
# revised April 2011
# adapted November 2020


setMethod("crosstab", signature(x="SpatRaster", y="missing"),
	function(x, digits=0, long=FALSE, useNA=FALSE) {

		nl <- nlyr(x)
		if (nl < 2) {
			error("crosstab", "needs at least 2 layers")
		}
		nms <- names(x)
		opt <- spatOptions()

		res <- x@pntr$crosstab(digits, !useNA, opt)
		if (length(res) == 0) {
			res <- data.frame(matrix(nrow=0, ncol=nl+1))
		} else {
			res <- matrix(res, ncol=nl+1, byrow=TRUE)
			res <- as.data.frame(res)
		}
		colnames(res) <- c(nms, "Freq")

		ff <- is.factor(x)
		if (any(ff)) {
			ff <- which(ff)
			v <- levels(x)
			for (i in ff) {
				j <- match(res[,i], v[[i]][,1])
				res[,i] <- v[[i]][j,2]
			}
		}

		if (!long) {
			nms <- make.names(nms, unique=TRUE)
			colnames(res) <- c(nms, "Freq")
			f <- eval(parse(text=paste("Freq ~ ", paste(nms , collapse="+"))))
			res <- stats::xtabs(f, data=res, addNA=useNA)
		} else {
			res <- res[res$Freq > 0,  ,drop=FALSE]
			# sorting is already done in C++
			rownames(res) <- NULL
			colnames(res)[ncol(res)] <- "n"
		}
		return(res)
	}
)