
loadModule("spat", TRUE)


.onAttach <- function(libname, pkgname) {
	tv <- utils::packageVersion("terra")
	m <- paste0("This is version ", tv, " of the \"terra\" package, for evaluation only\n")
	packageStartupMessage(m)
	
##############################
	.create_options()
	SpatRaster$new()$spatinit()
}

