# Author: Robert J. Hijmans
# Date : December 2017
# Version 1.0
# License GPL v3

setMethod("as.list", signature(x="SpatRaster"), 
function(x, ...)  {
	out <- list()
	for (i in 1:nlyr(x)) {
		out[[i]] <- x[[i]]
	}
	out
}
)


setMethod("overlay", signature(x="SpatRaster", y="SpatRaster"), 
function(x, y, fun, ..., filename="", overwrite=FALSE, wopt=list())  {
	
	stopifnot(!missing(fun))
	
	compareGeom(x, y, lyrs=FALSE, crs=TRUE, warncrs=TRUE, ext=TRUE, rowcol=TRUE, res=FALSE)
	
	nc <- ncol(x)
	readStart(x)
	readStart(y)
	
# figure out the shape of the output by testing with one row
	vx <- readValues(x, round(0.5*nrow(x)), 1, 1, nc, TRUE)
	vy <- readValues(y, round(0.5*nrow(y)), 1, 1, nc, TRUE)
	vxy <- cbind(vx, vy) # recycling
	test <- TRUE
	ncvx <- ncol(vx)
	vtst <- try(fun(vxy[,1:ncvx], vxy[,((ncvx+1):ncol(vxy))], ...), silent=TRUE)
	nl <- 1
	
	if (length(vtst) >= nc) {
		if ((length(vtst) %% nc) == 0) {
			dofun <- TRUE
			nl <- length(vtst) / nc
			test <- FALSE
		}
	}
	if (test) {
		vtst <- try(apply(vx, 1, fun, ...), silent=TRUE)
		if (length(vtst) >= nc) {		
			if ((length(vtst) %% nc) == 0) {
				dofun <- FALSE
				nl <- length(vtst) / nc
				# check for tranpose
				if (NROW(vtst) > 1) {
					transpose <- TRUE
				} else {
					transpose <- FALSE
				}
			} else {
				stop('overlay cannot use this function "fun"')
			}
		}
	}


	out <- rast(x, nlyr=nl)
	b <- writeStart(out, filename, overwrite, wopt)
	nc <- ncol(x)
#	nl <- nlyr(x)
#	fnames <- names(formals(fun))
#	if (length(fnames) != nl) {	dnames <- NULL 	} else { dnames <- list(list(), fnames)	}
	for (i in 1:b$n) {
		vx <- readValues(x, b$row[i], b$nrows[i], 1, nc)
		vy <- readValues(y, b$row[i], b$nrows[i], 1, nc)
		vx <- cbind(vx, vy) #recycling
		if (dofun) {
			r <- fun(vx[,1], vx[,2], ...)
		} else {
			r <- apply(vx, 1, fun, ...)
			if (transpose) {
				r <- as.vector(t(r))
			} 
		}
		writeVals(out, r, b$row[i], b$nrows[i])
	}
	readStop(x)
	out <- writeStop(out)
	return(out)
}
)

