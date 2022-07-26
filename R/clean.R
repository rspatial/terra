
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


clean_further <- function(x, tolerance=0.0001) {
	out <- as.lines(x)
	out <- snap(out, tolerance)
	out <- makeNodes(out)
	out <- mergeLines(out)
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

