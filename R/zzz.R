
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
	proj_path <- system.file("proj", package="terra")
	gdal_path <- system.file("gdal", package="terra")
	
	has_proj <- dir.exists(proj_path) && (file.exists(file.path(proj_path, "proj.db")) || file.exists(file.path(proj_path, "nad.lst")))
	has_gdal <- dir.exists(gdal_path)
	
	if (has_proj) {
		projPaths(proj_path, TRUE)
	}
	
	.gdalinit(path, ifelse(has_gdal, gdal_path, ""))
	
	if (libVersion("gdal") == "3.6.0") {
		message("Using GDAL version 3.6.0 which was retracted because it cannot write large GPKG files")
	}
	if (!proj_ok()) {
		message("There appears to be a problem with the PROJ installation")
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

