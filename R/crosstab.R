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

		b <- blocks(x, 4)
		readStart(x)
		on.exit(readStop(x))

		res <- NULL
		nc <- ncol(x)
		for (i in 1:b$n) {
			d <- readValues(x, b$row[i], b$nrows[i], 1, nc, TRUE)
			d <- lapply(1:nl, function(i) round(d[, i], digits=digits))
			d <- do.call(table, c(d, useNA="ifany"))
			d <- as.data.frame(d)
			res <- rbind(res, d)
		}

		res <- res[res$Freq > 0,  ,drop=FALSE]

		# some complexity to aggregate keeping
		# variables that are NA
		if (useNA) {
			for (i in 1:(ncol(res)-1)) {
				if (any(is.na(res[,i]))) {
					res[,i] <- factor(res[,i], levels=c(levels(res[,i]), NA), exclude=NULL)
				}
			}
		}
		res <- aggregate(res[, ncol(res), drop=FALSE], res[, 1:(ncol(res)-1), drop=FALSE], sum)

		for (i in 1:(ncol(res)-1)) {
		# get rid of factors
			res[,i] <- as.numeric(as.character(res[,i]))
		}
		if (nrow(res) == 0) {
			res <- data.frame(matrix(nrow=0, ncol=length(nms)+1))
		}
		nms <- make.names(nms, unique=TRUE)
		colnames(res) <- c(nms, "Freq")

		if (! useNA ) {
			i <- apply(res, 1, function(x) any(is.na(x)))
			res <- res[!i,  ,drop=FALSE]
		}

		ff <- is.factor(x)
		if (any(ff) && (digits >= 0)) {
			ff <- which(ff)
			v <- levels(x)
			for (i in ff) {
				j <- match(res[,i], v[[i]][,1])
				res[,i] <- v[[i]][j,2]
			}
		}

		if (!long) {
			f <- eval(parse(text=paste("Freq ~ ", paste(nms , collapse="+"))))
			res <- stats::xtabs(f, data=res, addNA=useNA)
		} else {
			res <- res[res$Freq > 0,  ,drop=FALSE]
			res <- res[order(res[,1], res[,2]), ]
			rownames(res) <- NULL
			colnames(res)[ncol(res)] <- "n"
		}
		return(res)
	}
)


