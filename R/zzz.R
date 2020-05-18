
loadModule("spat", TRUE)

gdal_version <- function() {
	.gdalversion()
}

.onAttach <- function(libname, pkgname) {
	#tv <- utils::packageVersion("terra")
	gdv <- gdal_version()
	if (gdv < "3.0.0") {
		m <- paste0("You are using GDAL", gdv, "\n For full functionality you need at least version 3.0.0")
		packageStartupMessage(m)
	}

##############################
	.create_options()
	
	path = ""
	if (file.exists(system.file("proj/nad.lst", package = "terra")[1])) {
		path <- system.file("proj", package="terra")
	} 
	SpatRaster$new()$spatinit(path)
}

