
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

	if (length(grep(.geos_version(FALSE, TRUE), .geos_version(TRUE))) != 1) {
		packageStartupMessage("WARNING: different compile-time and run-time versions of GEOS")
		packageStartupMessage("Compiled with:", .geos_version(FALSE, TRUE))
		packageStartupMessage(" Running with:", .geos_version(TRUE, TRUE))
		packageStartupMessage("\nYou should reinstall package 'terra'\n")
	}
#	terraOptions(todisk=TRUE, steps=2)
#	terraOptions(memfrac=0)
}

