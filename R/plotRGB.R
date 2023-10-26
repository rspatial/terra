# Author: Robert J. Hijmans
# Date :  April 2010
# Version 0.9
# License GPL v3

# ..linStretch <- function (x) {
    # v <- stats::quantile(x, c(0.02, 0.98), na.rm = TRUE)
    # temp <- (255 * (x - v[1]))/(v[2] - v[1])
    # temp[temp < 0] <- 0
    # temp[temp > 255] <- 255
    # return(temp)
# }

# # Histogram equalization stretch
# ..eqStretch <- function(x){
	# ecdfun <- stats::ecdf(x)
	# ecdfun(x)*255
# }

# ..rgbstretch <- function(RGB, stretch, caller="") {
	# stretch = tolower(stretch)
	# if (stretch == 'lin') {
		# RGB[,1] <- ..linStretch(RGB[,1])
		# RGB[,2] <- ..linStretch(RGB[,2])
		# RGB[,3] <- ..linStretch(RGB[,3])
	# } else if (stretch == 'hist') {
		# RGB[,1] <- ..eqStretch(RGB[,1])
		# RGB[,2] <- ..eqStretch(RGB[,2])
		# RGB[,3] <- ..eqStretch(RGB[,3])
	# } else if (stretch != '') {
		# warn(caller, "invalid stretch value")
	# }
	# RGB
# }


