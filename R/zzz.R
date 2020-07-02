

gdal_version <- function() {
	.gdalversion()
}

.init <- function() {
	path = ""
	if (file.exists(system.file("proj/nad.lst", package = "terra")[1])) {
		path <- system.file("proj", package="terra")
	} 
	s <- SpatRaster$new()
	s$spatinit(path)

	.create_options()

}

loadModule("spat", TRUE)

## we need the below, but this gives an error in R CMD check
#setLoadActions(
#	function(ns) {
#		try(.init(), silent=FALSE)
#	}
#) 


.onAttach <- function(libname, pkgname) {
	tv <- utils::packageVersion("terra")
	m <- paste("terra version", tv, "(alpha-release)")
	gdv <- gdal_version()
	if (gdv < "3.0.0") {
		a <- paste("\n\nNOTE: using GDAL version", gdv, "\nFor full functionality you need at least version 3.0.4\n")
		m <- c(m, a)
	}
	packageStartupMessage(m)

	.init()
}

