
loadModule("spat", TRUE)

.onAttach <- function(libname, pkgname) {
	tv <- utils::packageVersion("terra")
	m <- paste0("This is version ", tv, " of the \"terra\" package, for evaluation only\n")
	packageStartupMessage(m)
	
##############################
	.create_options()
	
	path = ""
	if (file.exists(system.file("proj/nad.lst", package = "terra")[1])) {
		path <- system.file("proj", package="terra")
	} 
	SpatRaster$new()$spatinit(path)
}

