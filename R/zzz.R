
.gdinit <- function() {
	path = ""
	if (file.exists(system.file("proj/nad.lst", package = "terra")[1])) {
		path <- system.file("proj", package="terra")
	} 
	.gdalinit(path)
}

loadModule("spat", TRUE)

#setLoadActions(
#	function(ns) {
#		.gdinit()
#	}
#)
 
.onLoad <- function(libname, pkgname) {
	.gdinit()
}


.onAttach <- function(libname, pkgname) {
	tv <- utils::packageVersion("terra")
	#m <- paste("terra version", tv, "(beta-release)")
	m <- paste("terra version", tv)
	gdv <- gdal_version()
	if (gdv < "3.0.0") {
		a <- paste("\n\nNOTE: using GDAL version", gdv, "\nFor full functionality you need at least version 3.0.4\n")
		m <- c(m, a)
	}
	packageStartupMessage(m)
	.create_options()
}

