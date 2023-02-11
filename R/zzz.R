
## from htmltools for future use with knitr
## method = c(package, genname, class).
#registerMethods <- function(methods) {
#	lapply(methods, function(method) {
#		pkg <- method[[1]]
#		generic <- method[[2]]
#		class <- method[[3]]
#		func <- get(paste(generic, class, sep="."))
#		if (pkg %in% loadedNamespaces()) {
#			registerS3method(generic, class, func, envir = asNamespace(pkg))
#		}
#		setHook(
#			packageEvent(pkg, "onLoad"),
#				function(...) {
#					registerS3method(generic, class, func, envir = asNamespace(pkg))
#				}
#			)
#		}
#	)
#}



.gdinit <- function() {
	path = ""
	sf <- system.file("", package="terra")
	if (file.exists(file.path(sf, "proj/nad.lst"))) {
		path <- system.file("proj", package="terra")
	}
	.gdalinit(path, file.path(sf, "gdal"))
	if (gdal() == "3.6.0") {
		message("Using GDAL version 3.6.0 which was retracted because it cannot write large GPKG files")
	}
}

loadModule("spat", TRUE)

.onLoad <- function(libname, pkgname) {
	.gdinit()
	#registerMethods(list(c(package, genname, class)))
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

