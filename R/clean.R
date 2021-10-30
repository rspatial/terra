
..gaps <- function(x) {
	p <- as.polygons(floor(ext(x)+1), crs=crs(x))
	e <- disagg(erase(p, x))
	if (nrow(e) > 1) {
		xmin = ext(p)[1]
		i <- sapply(1:nrow(e), function(i) ext(e[i])[1] > xmin)
		e <- e[i]
		x <- rbind(x, e)
	}
	x
}


width <- function(x) {
	x@ptr$width()
}

clearance <- function(x) {
	x@ptr$clearance()
}


removeDupNodes <- function(x, digits=-1) {
	x@ptr <- x@ptr$remove_duplicate_nodes(digits)
	messages(x, "removeDupNodes")	
}


simplify <- function(x, tolerance=0, preserveTopology=TRUE) {
	x@ptr <- x@ptr$simplify(tolerance, preserveTopology)
	messages(x, "simplify")	
}


make_nodes <- function(x) {
	x@ptr <- x@ptr$make_nodes()
	messages(x, "make_nodes")
}

line_merge <- function(x) {
	x@ptr <- x@ptr$line_merge()
	messages(x, "line_merge")
}

clean_further <- function(x, tolerance=0.0001) {
	out <- as.lines(x)
	out <- snap(out, tolerance)
	out <- make_nodes(out)
	out <- line_merge(out)
	as.polygons(out)
}

clean <- function(x) {
	g <- gaps(x)
	out <- erase(x)
	out <- rbind(out, g)
	g <- gaps(out)
	rbind(out, g)
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

snap <- function(x, tolerance) {
	x@ptr <- x@ptr$snap(tolerance)
	messages(x, "snap")
}


centerline <- function(p) {
	v <- as.points(voronoi(p, tolerance=0))
	v <- intersect(v, p)
	v
}



#library(terra); messages = terra:::messages
#p <- vect(system.file("ex/lux.shp", package="terra"))
#h <- convHull(p[-12], "NAME_1")
#x <- clean(h)
#y <- clean_further(x)

#hh <- rbind(h, h)
#e <- erase(hh)
#g <- gaps(e)

#v1 = as.polygons(ext(0,1,0,1))
#v2 = as.polygons(ext(1.01,2,0,1))
#v <- rbind(v1, v2)
#s = snap(v)

