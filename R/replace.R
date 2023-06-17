# # Author: Robert J. Hijmans
# # Date:  October 2018
# # Version 1.0
# # License GPL v3


# .rast_replace <- function(x, name, value, caller="$<-") {
	# if (inherits(value, "SpatRaster")) {
		# value <- value[[1]]
		# names(value) <- name
	# } else if (!is.null(value)) {
		# y <- rast(x, nlyrs=1)
		# test <- try(values(y) <- value, silent=TRUE)
		# if (inherits(test, "try-error")) {
			# error(caller, "the replacement value is not valid")
		# }
		# value <- y
		# names(value) <- name
	# }
	# i <- which(name == names(x))[1]
	# if (is.null(value)) {
		# if (is.na(i)) {
				# return(x)
		# } else {
			# return(subset(x, -i, NSE=FALSE))
		# }
	# }

	# if (is.na(i)) {
		# c(x, value)
	# } else if (nlyr(x) == 1) {
		# value$deepcopy()
	# } else if (i == 1) {
		# c(value, x[[2:nlyr(x)]])
	# } else if (i == nlyr(x)) {
		# c(x[[1:(nlyr(x)-1)]], value)
	# } else {
		# c(x[[1:(i-1)]], value, x[[(i+1):nlyr(x)]])
	# }
# }

# setMethod("$<-", "SpatRaster",
	# function(x, name, value) {
		# .rast_replace(x, name, value, "$<-")
	# }
# )


# setReplaceMethod("[[", c("SpatRaster", "character", "missing"),
	# function(x, i, j, value) {
		# if (inherits(value, "SpatRaster")) {
			# if (nlyr(value) != length(i)) {
				# error("`[[`", "length of names must be equal to the number of layers")
			# }
			# names(value) <- i
		# } else if (length(i) > 1) {
			# if (NCOL(value) > 1) {
				# value <- as.list(data.frame(value))
			# } else {
				# stopifnot(length(i) == length(value))
			# }
		# } else if (!is.list(value)) {
			# value <- list(value)
		# }
		# for (k in 1:length(i)) {
			# x <- .rast_replace(x, i[k], value[[k]], " [[<- ")
# #			eval(parse(text = paste0("x$", i[k], " <- value[[k]]")))
		# }
		# x
	# }
# )

# setReplaceMethod("[[", c("SpatRaster", "numeric", "missing"),
	# function(x, i, j, value) {
		# if (!inherits(value, "SpatRaster")) {
			# error("`[[<-`", "Expected a SpatRaster as replacement value")
		# }
		# if (nlyr(value) != length(i)) {
			# error("`[[`", "length of indices must be equal to the number of layers")
		# }
		# if (any(i<1) | any(i > nlyr(x))) {
			# error("`[[`", "indices must be between 1 and the number of layers")
		# }
		# if (nlyr(x) == 1) {
			# compareGeom(x, value, crs=FALSE, warncrs=TRUE)
			# return(value)
		# }
		# for (k in 1:length(i)) {
			# if (i[k] == 1) {
				# x <- c(value[[k]], x[[2:nlyr(x)]])
			# } else if (i[k] == nlyr(x)) {
				# x <- c(x[[1:(nlyr(x)-1)]], value[[k]])
			# } else {
				# x <- c(x[[1:(i[k]-1)]], value[[k]], x[[(i[k]+1):nlyr(x)]])
			# }
		# }
		# g <- gc()
		# x
	# }
# )



# setReplaceMethod("[", c("SpatRaster", "missing", "missing"),
	# function(x, i, j, value) {

		# nl <- nlyr(x)
		# if (is.matrix(value)) {
			# d <- dim(value)
			# if (!all(d == c(ncell(x), nl))) {
				# if ((d[2] == nl) && (d[1] < ncell(x))) {
					# value <- apply(value, 2, function(i) rep_len(i, ncell(x)))
				# } else {
					# error("`[`","dimensions of the matrix do not match the SpatRaster")
				# }
			# }
			# x <- try( setValues(x, value, TRUE, TRUE) )
		# } else {
			# v <- try( matrix(nrow=ncell(x), ncol=nl) )
			# if (! inherits(x, "try-error")) {
				# v[] <- value
				# x <- try( setValues(x, v, TRUE, TRUE) )
			# }
		# }
		# if (inherits(x, "try-error")) {
			# error("`[`", "cannot set values")
		# }
		# return(x)
	# }
# )



# setReplaceMethod("[", c("SpatRaster","numeric", "missing"),
	# function(x, i, j, value) {
		# theCall <- sys.call(-1)
		# narg <- length(theCall)-length(match.call(call=sys.call(-1)))
		# if (narg > 0) { # row
			# i <- cellFromRowColCombine(x, i, 1:ncol(x))
		# }
		# #if (any(is.na(i))) {
		# #	warn("`[`", "indices should not be NA")
		# #}
		# bylyr = FALSE
		# if (!is.null(dim(value))) {
			# #x@ptr <- x@ptr$replaceValues(i, value, ncol(value))
			# stopifnot(ncol(value) == nlyr(x))
			# bylyr <- TRUE
			# if (inherits(value, "data.frame")) {
				# value <- as.matrix(value)
			# }
			# value <- as.vector(value)
		# }

		# x@ptr <- x@ptr$deepcopy()
		# opt <- spatOptions()
		# if (!x@ptr$replaceCellValues(i-1, value, bylyr, opt)) {
			# messages(x, "`[<-`")
		# } else {
			# x
		# }
	# }
# )


# setMethod("set.values", signature(x="SpatRaster"),
	# function(x, cells, values, layer=0)  {

		# #if (any(is.na(cells))) {
		# #	warn("set.values", "cells should not be NA")
		# #}

		# if (is.character(layer)) {
			# layer <- match(layer, names(x))
			# if (any(is.na(layer))) {
				# error("set.values", "invalid layer")
			# }
		# }
		# layer <- round(layer)

		# if (all(layer > 0)) {
			# if (missing(cells) && missing(values)) {
				# return(invisible(TRUE));
			# }
			# if (any(is.na(layer))) { error("set.values", "layers cannot be NA")}
			# if (inherits(layer, "character")) {
				# layer <- match(layer, names(x))
				# if (any(is.na(layer))) { error("set.values", "invalid layer names")}
			# }
			# if (any((layer < 1) | (layer > nlyr(x))))  { error("set.values", "invalid layer numbers") }
			# n <- length(layer)
			# if (n > length(unique(layer)))  { error("set.values", "duplicated layers") }

			# bylyr <- FALSE
			# if (!is.null(dim(values))) {
				# if (ncol(values) != n) {
					# error("set.values", "ncol(values) does not match the `length(layer)`")
				# }
				# bylyr <- TRUE
				# #if (inherits(values, "data.frame")) {
				# #	values <- as.matrix(values)
				# #}
				# values <- as.vector(values)
			# }
			# ok <- x@ptr$replaceCellValuesLayer(layer-1, cells-1, values, bylyr, spatOptions())
			# messages(x)
			# invisible(TRUE)
		# } else {
			# if (any(layer > 0)) {
				# error("set.values", "some (but not all) layer numbers are < 1")
			# }
			# if (missing(cells) && missing(values)) {
				# x@ptr$readAll()
				# return(invisible(TRUE));
			# }
			# bylyr <- FALSE
			# if (!is.null(dim(values))) {
				# if (ncol(values) != nlyr(x)) {
					# error("set.values", "ncol(values) does not match the nlyr(x)")
				# }
				# bylyr <- TRUE
				# #if (inherits(values, "data.frame")) {
				# #	values <- as.matrix(values)
				# #}
				# values <- as.vector(values)
			# }
			# ok <- x@ptr$replaceCellValues(cells-1, values, bylyr, spatOptions())
			# messages(x)
		# }
		# invisible(TRUE)
	# }
# )

# setReplaceMethod("[", c("SpatRaster", "numeric", "numeric"),
	# function(x, i, j, value) {
		# i <- cellFromRowColCombine(x, i, j)
		# x[i] <- value
		# x
	# }
# )



# setReplaceMethod("[", c("SpatRaster","missing", "numeric"),
	# function(x, i, j, value) {
		# i <- cellFromRowColCombine(x, 1:nrow(x), j)
		# x[i] <- value
		# x
	# }
# )



# setReplaceMethod("[", c("SpatRaster", "logical", "missing"),
	# function(x, i, j, value) {
		# i <- which(rep_len(i, ncell(x)))
		# x[i] <- value
		# x
	# }
# )


# setReplaceMethod("[", c("SpatRaster", "SpatRaster", "ANY"),
	# function(x, i, j, value) {
		# theCall <- sys.call(-1)
		# narg <- length(theCall)-length(match.call(call=sys.call(-1)))
		# if (narg > 0) { # row
			# error("`[`", "you cannot use a SpatRaster as a row index")
		# }
		# if (inherits(value, "SpatRaster")) {
			# x <- mask(x, i, maskvalues=TRUE)
			# cover(x, value)
		# } else {
			# if (NCOL(value) > 1) {
				# error(" [", "cannot use a data.frame with multiple columns")
			# }
			# value <- unlist(value)
			# if (length(value) == 1) {
				# mask(x, i, maskvalues=TRUE, updatevalue=value[1])
			# } else {
				# i <- as.logical(values(i))
				# i[is.na(i)] <- TRUE
				# i <- which(i)
				# x[i] <- value
				# x
			# }
		# }
	# }
# )


# setReplaceMethod("[", c("SpatRaster", "SpatVector", "missing"),
	# function(x, i, j, value) {
		# theCall <- sys.call(-1)
		# narg <- length(theCall)-length(match.call(call=sys.call(-1)))
		# if (narg > 0) { # row
			# error("`[`", "you cannot use a SpatVector as a row index")
		# }
		# if (length(value) > 1) {
			# value <- rep_len(value, length.out=length(x))
		# }
		# rasterize(i, x, field=value, update=TRUE)
	# }
# )