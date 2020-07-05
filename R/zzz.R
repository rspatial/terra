

.init <- function() {
	path = ""
	if (file.exists(system.file("proj/nad.lst", package = "terra")[1])) {
		path <- system.file("proj", package="terra")
	} 
	s <- SpatRaster$new()
	s$spatinit(path)

	.create_options()

}

.init2 <- function() {
	path = ""
	if (file.exists(system.file("proj/nad.lst", package = "terra")[1])) {
		path <- system.file("proj", package="terra")
	} 
	.gdalinit(path)
}



loadModule("spat", TRUE)

## we need the below, but this gives an error in R CMD check
setLoadActions(
	function(ns) {
#		print(.gdalversion())
		.init2()
#		try(.init(), silent=FALSE)
	}
) 


.onAttach <- function(libname, pkgname) {
	tv <- utils::packageVersion("terra")
	m <- paste("terra version", tv, "(beta-release)")
	gdv <- gdal_version()
	if (gdv < "3.0.0") {
		a <- paste("\n\nNOTE: using GDAL version", gdv, "\nFor full functionality you need at least version 3.0.4\n")
		m <- c(m, a)
	}
	packageStartupMessage(m)

	.init()
}

