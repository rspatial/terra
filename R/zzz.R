
.gdinit <- function() {
	path = ""
	if (file.exists(system.file("proj/nad.lst", package = "terra")[1])) {
		path <- system.file("proj", package="terra")
	} 
	.gdalinit(path)
}

loadModule("spat", TRUE)


.onLoad <- function(libname, pkgname) {
	.gdinit()
}


.onAttach <- function(libname, pkgname) {
	packageStartupMessage("terra ", utils::packageVersion("terra"))
	.create_options()
#	terraOptions(todisk=TRUE)
}

