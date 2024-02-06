# Author: Robert Hijmans
# November 2009
# License GPL3


movingFun <- function(x, n, fun=mean, type="around", circular=FALSE, na.rm=FALSE, ...)  { 
	n <- round(abs(n))
    if (n == 0) { stop("n == 0")  }
    x = as.vector(x)
	lng <- length(x)
	if (type == "around") {
		hn <- floor(n/2)
		if (circular) {
			x <- c(x[(lng-hn+1):lng], x, x[1:hn])
		} else { 
			x <- c(rep(NA, hn), x, rep(NA, hn)) 
		}
	} else if (type == "to") {
		if (circular) { 
			x <- c(x[(lng-n+2):lng], x)
		} else { 
			x <- c(rep(NA, n-1), x) 
		}
	} else if (type == "from") {
		if (circular) { 
			x <- c(x,  x[1:n])
		} else { 
			x <- c(x, rep(NA, n))	
		}
	} else {
		stop('unknown type; should be "around", "to", or "from"')
	}
	m <- matrix(ncol=n, nrow=lng)
    for (i in 1:n) { m[,i] <- x[i:(lng+i-1)] }
	if (na.rm) {
		apply(m, MARGIN=1, FUN=fun, na.rm=na.rm, ...)
	} else {
		apply(m, MARGIN=1, FUN=fun, ...)	
	}
}


.roll <- function(x, n) {
# by Josh O'Brien
    x[(seq_along(x) - (n+1)) %% length(x) + 1]
}


setMethod("roll", signature(x="numeric"),
	function(x, n, fun=mean, type="around", circular=FALSE, na.rm=FALSE, ...) {
		movingFun(x, n, fun, type=type, circular=circular, na.rm=na.rm, ...)
	}
)



setMethod("roll", signature(x="SpatRaster"),
	function(x, n, fun="mean", type="around", circular=FALSE, na.rm=FALSE, filename="", ..., wopt=list()) {
		txtfun <- .makeTextFun(match.fun(fun))
		if (inherits(txtfun, "character")) {
			if (txtfun %in% .cpp_funs) {
				opt <- spatOptions(filename, ...)
				x@ptr <- x@ptr$roll(n, txtfun, type, circular, na.rm, opt)
				return (messages(x, "roll")	)
			}
		} else {
			f <- function(i) {
				movingFun(i, n, fun, type=type, circular=circular, na.rm=na.rm, ...)
			}
			app(x, f, filename=filename, ..., wopt=wopt)
		}
	}
)

