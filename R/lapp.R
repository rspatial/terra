# Author: Robert J. Hijmans
# Date : May 2020
# Version 1.0
# License GPL v3


.getNL <- function(v, fun, usenames, ...) {
# figure out the shape of the output
	nr <- nrow(v)
	if (!usenames) colnames(v) <- NULL
	vtst <- try(do.call(fun, c(v, list(...))), silent=TRUE)
	if (inherits(vtst, "try-error")) {
		return(-1)
	}
	if (length(vtst) >= nr) {
		if ((length(vtst) %% nr) == 0) {
			nl <- length(vtst) / nr
		} else {	
			return(-1)
		}
	} else {
		return(-1)
	}
	nl
}


setMethod("lapp", signature(x="SpatRaster"), 
function(x, fun, ..., usenames=FALSE, filename="", overwrite=FALSE, wopt=list())  {
	
	stopifnot(!missing(fun))
	if (usenames) {
		fnames <- names(formals(fun))
		x <- x[[names(x) %in% fnames]]
	}
	readStart(x)
	ncx <- ncol(x)
	v <- readValues(x, round(0.5*nrow(x)), 1, 1, ncx, dataframe=TRUE)
	nl <- .getNL(v, fun, usenames, ...)
	if (nl < 1) stop("lapp does not like 'fun'")
	out <- rast(x, nlyr=nl)
	b <- writeStart(out, filename, overwrite, wopt)
	for (i in 1:b$n) {
		v <- readValues(x, b$row[i], b$nrows[i], 1, ncx, dataframe=TRUE)
		if (!usenames) colnames(v) <- NULL
		v <- do.call(fun, c(v, list(...)))
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	readStop(x)
	out <- writeStop(out)
	return(out)
}
)


.getNLstack <- function(v, fun, recycle, ...) {
# figure out the shape of the output
	nr <- nrow(v[[1]])
	if (recycle) {
		v <- lapply(v, as.vector)
	}
	vtst <- try(do.call(fun, c(v, list(...))), silent=TRUE)
	if (inherits(vtst, "try-error")) {
		return(-1)
	}
	if (length(vtst) >= nr) {
		if ((length(vtst) %% nr) == 0) {
			nl <- length(vtst) / nr
		} else {	
			return(-1)
		}
	} else {
		return(-1)
	}
	nl
}



setMethod("lapp", signature(x="SpatStack"), 
function(x, fun, ..., recycle=FALSE, filename="", overwrite=FALSE, wopt=list())  {
	
	stopifnot(!missing(fun))
	
	ncx <- ncol(x[1])
	nrx <- nrow(x[1])
	readStart(x)
	v <- lapply(1:nsds(x), function(i) readValues(x[i], round(0.5*nrx), 1, 1, ncx, mat=TRUE))
	nl <- .getNLstack(v, fun, recycle, ...)
	if (nl < 1) stop("lapp does not like 'fun'")
	out <- rast(x[1], nlyr=nl)
	b <- writeStart(out, filename, overwrite, wopt)
	for (i in 1:b$n) {
		v <- lapply(1:nsds(x), function(s) readValues(x[s], b$row[i], b$nrows[i], 1, ncx, mat=TRUE))
		if (recycle) {
			v <- lapply(v, as.vector)
		}
		v <- do.call(fun, c(v, list(...)))
		writeValues(out, v, b$row[i], b$nrows[i])
	}
	readStop(x)
	out <- writeStop(out)
	return(out)
}
)




# setMethod("overlay", signature(x="SpatRaster", y="SpatRaster"), 
# function(x, y, fun, ..., filename="", overwrite=FALSE, wopt=list())  {
	
	# stopifnot(!missing(fun))
	
	# compareGeom(x, y, lyrs=FALSE, crs=TRUE, warncrs=TRUE, ext=TRUE, rowcol=TRUE, res=FALSE)
	
	# nc <- ncol(x)
	# readStart(x)
	# readStart(y)
	
# # figure out the shape of the output by testing with one row
	# vx <- readValues(x, round(0.5*nrow(x)), 1, 1, nc, TRUE)
	# vy <- readValues(y, round(0.5*nrow(y)), 1, 1, nc, TRUE)
	# vxy <- cbind(vx, vy) # recycling
	# test <- TRUE
	# ncvx <- ncol(vx)
	# vtst <- try(fun(vxy[,1:ncvx], vxy[,((ncvx+1):ncol(vxy))], ...), silent=TRUE)
	# nl <- 1
	
	# if (length(vtst) >= nc) {
		# if ((length(vtst) %% nc) == 0) {
			# dofun <- TRUE
			# nl <- length(vtst) / nc
			# test <- FALSE
		# }
	# }
	# if (test) {
		# vtst <- try(apply(vx, 1, fun, ...), silent=TRUE)
		# if (length(vtst) >= nc) { # should be nrow?		
			# if ((length(vtst) %% nc) == 0) {
				# dofun <- FALSE
				# nl <- length(vtst) / nc
				# # check for tranpose
				# if (NROW(vtst) > 1) {
					# transpose <- TRUE
				# } else {
					# transpose <- FALSE
				# }
			# } else {
				# stop('overlay cannot use this function "fun"')
			# }
		# }
	# }


	# out <- rast(x, nlyr=nl)
	# b <- writeStart(out, filename, overwrite, wopt)
	# nc <- ncol(x)
# #	nl <- nlyr(x)
# #	fnames <- names(formals(fun))
# #	if (length(fnames) != nl) {	dnames <- NULL 	} else { dnames <- list(list(), fnames)	}
	# for (i in 1:b$n) {
		# vx <- readValues(x, b$row[i], b$nrows[i], 1, nc)
		# vy <- readValues(y, b$row[i], b$nrows[i], 1, nc)
		# vx <- cbind(vx, vy) #recycling
		# if (dofun) {
			# r <- fun(vx[,1], vx[,2], ...)
		# } else {
			# r <- apply(vx, 1, fun, ...)
			# if (transpose) {
				# r <- as.vector(t(r))
			# } 
		# }
		# writeValues(out, r, b$row[i], b$nrows[i])
	# }
	# readStop(x)
	# out <- writeStop(out)
	# return(out)
# }
# )

