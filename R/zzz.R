
loadModule("spat", TRUE)


.onAttach <- function(libname, pkgname) {
	minv <- "2.8.7"
	m <- "\n\nThis is an early version of the \"terra\" package, for evaluation only\n"
	e <- try(vv <- utils::packageVersion("raster"), silent=T )
	if (class(e)[1] == "try-error") {
		m <- paste0(m,"For this package to work you need to install version ", minv, 
		" or higher of the \"raster\" package\n" )	
	} else if (vv < minv ) {
		m <- paste0(m,"For this package to work you need to install version ", minv, 
		" or higher of the \"raster\" package\nYou have version: ", vv,"\n" )
	} 
	packageStartupMessage(m)
}


