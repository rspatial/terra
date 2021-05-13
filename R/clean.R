

simplify <- function(x, tolerance=0, preserveTopology=TRUE) {
	x@ptr <- x@ptr$simplify(tolerance, preserveTopology)
	messages(x)	
}


make_nodes <- function(x) {
	x@ptr <- x@ptr$make_nodes()
	messages(x)
}

line_merge <- function(x) {
	x@ptr <- x@ptr$line_merge()
	messages(x)
}

clean_further <- function(x) {
	out <- as.lines(x)
	out <- make_nodes(out)
	out <- line_merge(out)
	as.polygons(out)
}


clean <- function(x) {
	x <- clean_further(x)

	out = x[1]
	for (i in 2:nrow(x)) {
		out <- erase(out, x[i])
		out <- rbind(out, x[i])
	}
	out <- clean_further(out)
	p <- as.polygons(floor(ext(out)+1), crs=crs(out))
	e <- disaggregate(erase(p, out))
	if (nrow(e) > 1) {
		xmin = ext(p)[1]
		i <- sapply(1:nrow(e), function(i) ext(e[i])[1] > xmin)
		e <- e[i]
		out <- rbind(out, e)
	}
	out
}


mergebyborder <- function(x, field) {
	i <- is.na(x[[field, drop=TRUE]])
	if (!any(i)) return(x)
	s <- sharedPaths(x)
	s$length <- perim(s)
	from <- x[i]
	x <- x[!i]

	for (i in 1:nrow(from)) {

	}
} 


#p <- vect(system.file("ex/lux.shp", package="terra"))
#h = convHull(p[-12], "NAME_1")
#x <- clean(h)
#y <- clean_further(x)


