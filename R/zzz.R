
loadModule("spat", TRUE)


.onAttach <- function(libname, pkgname) {
	tv <- utils::packageVersion("terra")
	m <- paste0("This is version ", tv, " of the \"terra\" package, for evaluation only\n")
	
	#min_raster <- "3.0-5"
	#e <- try(vv <- utils::packageVersion("raster"), silent=TRUE)
	#if (class(e)[1] == "try-error") {
	#	m <- paste0(m,"For this package to work you need to install version ", min_raster, 
	#	" or higher of the \"raster\" package\n" )	
	#} else if (vv < min_raster) {
	#	m <- paste0(m,"For this package to work you need to install version ", min_raster, 
	#	" or higher of the \"raster\" package\nYou have version: ", vv,"\n" )
	#} 

	packageStartupMessage(m)
	
##############################
	.create_options()
	x <- SpatRaster$new()
	x$spatinit()
}


