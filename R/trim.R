# Author: Robert J. Hijmans
# Date : December 2009
# Version 1.0
# Licence GPL v3

.mytrimws <- function(x, internal=FALSE, ...) {
	if (internal) {
		gsub("^ *|(?<= ) | *$", "", x, perl=TRUE)
	} else {
		gsub("^\\s+|\\s+$", "", x)
	}
}



# setMethod('trim', signature(x='data.frame'), 
	# function(x, ...) {
		# for (i in 1:ncol(x)) {
			# if (class(x[,i]) == 'character') {
				# x[,i] <- trim(x[,i])
			# } else if (class(x[,i]) == 'factor') {
				# x[,i] <- as.factor(trim(as.character(x[,i])))
			# }	
		# }
		# return(x)
	# }
# )


# setMethod('trim', signature(x='matrix'), 
	# function(x, ...) {
		# if (is.character(x)) {
			# x[] = trim(as.vector(x))
		# } else {
			# rows <- rowSums(is.na(x))
			# cols <- colSums(is.na(x))
			# rows <- which(rows != ncol(x))
			# cols <- which(cols != nrow(x))
			# if (length(rows)==0) {
				# x <- matrix(ncol=0, nrow=0)
			# } else {
				# x <- x[min(rows):max(rows), min(cols):max(cols), drop=FALSE]
			# }
		# }
		# return(x)
	# }
# )



setMethod('trim', signature(x='SpatRaster'), 
function(x, padding=0, filename='', ...) {

	overwrite <- .overwrite(...)
	r <- methods::new('SpatRaster')
	ptr <- try(x@ptr$trim(padding, filename, overwrite))
	if (class(ptr) == 'try-error') { stop("trim error") } else { r@ptr <- ptr }
	r
}
)

