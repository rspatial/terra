# Author: Robert J. Hijmans
# Date : April 2015
# Version 1.0
# Licence GPL v3


# setMethod("rowSums", signature(x="SpatRaster"), 
	# function(x, na.rm = FALSE, dims = 1L, ...) {
	
		# nl <- nlyr(x)
		# nc <- ncol(x)
		# readStart(x)
		# on.exit(readStop(x))
		# b <- blocks(x, n=4)
		# s <- list()
		# for (i in 1:b$n) {
			# v <- readValues(x, row=b$row[i], nrows=b$nrows[i])
			# s[[i]] <- .colSums(v, nc, b$nrows[i]*nl, na.rm=na.rm, ...)
		# }
		# s <- t(matrix(unlist(s), nrow=nl, byrow=TRUE))
		# colnames(s) <- names(x)
		# s
	# }
# )


# setMethod("colSums", signature(x="SpatRaster"), 
	# function(x, na.rm = FALSE, dims = 1L, ...) {
		# nl <- nlyr(x)
		# nc <- ncol(x)
		# readStart(x)
		# on.exit(readStop(x))
		# b <- blocks(x, n=4)
		# s <- matrix(nrow=b$n, ncol=nc*nl) 
		# for (i in 1:b$n) {
			# v <- readValues(x, row=b$row[i], nrows=b$nrows[i], mat=TRUE)
			# for (j in 1:nl) {
				# k <- (j-1) * nc + 1
				# k <- k:(k+nc-1)
				# s[i, k] <- .colSums(matrix(v[,j], nrow=b$nrows[i], byrow=TRUE), b$nrows[i], nc, na.rm=na.rm, ...)
			# }
		# }
		# s <- matrix(.colSums(s, nrow(s), ncol(s), na.rm=na.rm), ncol=nl)
		# colnames(s) <- names(x)
		# return(s)
	# }
# )



# setMethod("rowMeans", signature(x="SpatRaster"), 
	# function(x, na.rm = FALSE, dims = 1L, ...) {
	
		# nl <- nlyr(x)
		# nc <- ncol(x)
		# readStart(x)
		# on.exit(readStop(x))
		# b <- blocks(x, n=4)
		# s <- list()
		# for (i in 1:b$n) {
			# v <- readValues(x, row=b$row[i], nrows=b$nrows[i])
			# s[[i]] <- .colMeans(v, nc, b$nrows[i]*nl, na.rm=na.rm, ...)
		# }
		# s <- t(matrix(unlist(s), nrow=nl, byrow=TRUE))
		# colnames(s) <- names(x)
		# s
	# }
# )


# setMethod("colMeans", signature(x="SpatRaster"), 
	# function(x, na.rm = FALSE, dims = 1L, ...) {
		# nl <- nlyr(x)
		# nc <- ncol(x)
		# readStart(x)
		# on.exit(readStop(x))
		# b <- blocks(x, n=4)
		# s <- matrix(nrow=b$n, ncol=nc*nl) 
		# for (i in 1:b$n) {
			# v <- readValues(x, row=b$row[i], nrows=b$nrows[i], mat=TRUE)
			# for (j in 1:nl) {
				# k <- (j-1) * nc + 1
				# k <- k:(k+nc-1)
				# s[i, k] <- .colMeans(matrix(v[,j], nrow=b$nrows[i], byrow=TRUE), b$nrows[i], nc, na.rm=na.rm, ...)
			# }
		# }
		# s <- matrix(.colSums(s, nrow(s), ncol(s), na.rm=na.rm), ncol=nl)
		# colnames(s) <- names(x)
		# return(s)
	# }
# )


