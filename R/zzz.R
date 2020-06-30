
loadModule("spat", TRUE)

gdal_version <- function() {
	.gdalversion()
}


.onAttach <- function(libname, pkgname) {
	tv <- utils::packageVersion("terra")
	m <- paste("This is terra version", tv, "(alpha-release)")
	gdv <- gdal_version()
	if (gdv < "3.0.0") {
		a <- paste("\n\nNOTE: You are using GDAL version", gdv, "\nFor full functionality you need at least version 3.0.4\n")
		m <- c(m, a)
	}
	packageStartupMessage(m)

	.create_options()
	
	path = ""
	if (file.exists(system.file("proj/nad.lst", package = "terra")[1])) {
		path <- system.file("proj", package="terra")
	} 
	s <- SpatRaster$new()
	s$spatinit(path)

}

