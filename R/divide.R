# Author: Robert Hijmans
# December 2024
# License GPL3


setMethod("divide", signature(x="SpatRaster"),
	function(x, n=2, start="ns", as.raster=FALSE, na.rm=TRUE) {

		if (nlyr(x) > 1) {
			warn("divide", "only the first layer is used")
			x <- x[[1]]
		}

		n <- round(n)
		north <- TRUE
		if (length(n) > 1) {
			if (any(n == 0)) {
				error("divide", "n values cannot be 0")
			}
		} else {
			if (!isTRUE(n > 0)) {
				error("divide", "n must be > 0")
			}
			start <- match.arg(tolower(start), c("ns", "ew"))
			north <- start == "ns"
		}

		if (!hasValues(x)) {
			out <- x
			if (length(n) > 1) {
				nrow(out) <- sum(abs(n[n<0]) + 1)
				ncol(out) <- sum(abs(n[n<0]) + 1)
			} else {
				if (north) {
					rows <- ceiling(n / 2)
					cols <- n - rows
				} else {
					cols <- ceiling(n / 2)
					rows <- n - cols
				}
				nrow(out) <- rows * 2
				ncol(out) <- cols * 2
			}
			out <- as.polygons(out)
			out$zones <- 1:nrow(out)
			if (isTRUE(as.raster) || is.na(as.raster)) {
				r <- rasterize(out, x, "zones")
				if (is.na(as.raster)) {
					return(list(r=r, v=out))
				} else {
					return(r)
				}
			}
			return(out)
		}

		opt <- spatOptions()
		nn <- as.integer(n)
		if (isTRUE(as.raster) || is.na(as.raster)) {
			r <- x
			r@pntr <- x@pntr$divide(nn, north, opt)
			r <- messages(r, "divide")
			if (na.rm) {
				r <- mask(r, x)
			}
			if (isTRUE(as.raster)) {
				return(r)
			}
		}
		out <- methods::new("SpatVector")
		out@pntr <- x@pntr$divide_vector(nn, north, opt)
		out <- messages(out, "divide")
		if (na.rm) {
			border <- as.polygons(not.na(x, TRUE))
			out <- crop(out, border)
		}
		if (is.na(as.raster)) {
			return(list(r=r, v=out))
		}
		out
	}
)



check_frac <- function(f) {
	if (is.null(f)) return(f)
	if ((any(f <= 0)) || (any(f >= 1))) {
		stop("f values must be > 0 and < 1")
	}
	if (length(f) == 1) {
		f <- seq(f, 1, f)
		f[length(f)] <- 1
	} else {
		if (length(unique(f)) < length(f)) {
			stop("f values must be unique")	
		}
		ord <- order(f)
		if (!all(ord == 1:length(ord))) {
			stop("f values must be in ascending order")		
		}
		if (!isTRUE(all.equal(sum(f), 1, tolerance=0.0001))) {
			stop("f values must sum to 1")		
		}
	}
	f
}



strip_polygon <- function(x, vertical, horizontal) {

## based on a function by Barry Rowlinson 

	totalArea <- expanse(x, transform=FALSE, unit="km")
	e <- ext(x)
	ex <- data.frame(t(as.vector(e + 1)))
	e <- data.frame(t(as.vector(e)))

	if (!is.null(vertical)) {
		edges <- sapply(vertical, function(fraction){
			target <- totalArea * fraction
			target_fun <- function(xm){
				expanse(crop(x, ext(ex$xmin, xm, ex$ymin, ex$ymax)), transform=FALSE, unit="km") - target
			}
			stats::uniroot(target_fun, lower=e$xmin+0.0000001, upper=e$xmax)$root
		})
		xbnds <- matrix(c(ex$xmin, rep(edges,rep(2,length(edges))), ex$xmax), ncol=2, byrow=TRUE)
		xbnds <- cbind(xbnds, ex$ymin, ex$ymax)
		xbnds <- do.call(rbind, apply(xbnds, 1, function(i) as.polygons(ext(i))))
		xbnds$vid <- 1:nrow(xbnds)
	}
	if (!is.null(horizontal)) {
		edges <- sapply(horizontal, function(fraction){
			target <- totalArea * fraction
			target_fun <- function(ym){
				expanse(crop(x, ext(ex$xmin, ex$xmax, ex$ymin, ym)), transform=FALSE, unit="km") - target
			}
			stats::uniroot(target_fun, lower=e$ymin+0.0000001, upper=e$ymax)$root
		})
		ybnds <- matrix(c(ex$ymin, rep(edges,rep(2,length(edges))), ex$ymax), ncol=2, byrow=TRUE)
		ybnds <- cbind(ex$xmin, ex$xmax, ybnds)
		ybnds <- do.call(rbind, apply(ybnds, 1, function(i) as.polygons(ext(i))))
		ybnds$hid <- 1:nrow(ybnds)

		if (!is.null(vertical)) {
			bnds <- union(xbnds, ybnds)
		} else {
			bnds <- ybnds
		}	
	}
	intersect(x, bnds)
}	


divide_polygon <- function(x, n, w, alpha, ...) {
	xcrs <- crs(x) 
	crs(x) <- "+proj=utm +zone=1"
	s <- terra::spatSample(x, max(n*4, 1000, log(n) * 100), "regular")
	xy <- terra::crds(s)
	if (!is.null(w)) {
		e <- extract(w, s, ID=FALSE)
		alpha <- rep_len(alpha, 2)
		xy[,1] <- xy[,1] * alpha[1]
		xy[,2] <- xy[,2] * alpha[2]
		xy <- na.omit(cbind(xy, e))
	} 
	k <- stats::kmeans(xy, centers = n, ...)
	ctrs <- k$centers[, 1:2]
	if (!is.null(w)) {
		ctrs[,1] <- ctrs[,1] / alpha[1]
		ctrs[,2] <- ctrs[,2] / alpha[2]
	}
	v <- terra::voronoi(vect(ctrs, crs=xcrs), bnd=x)
	terra::crop(v, x)
}


setMethod("divide", signature(x="SpatVector"),
	function(x, n=5, w=NULL, alpha=1, ...) {
		if (geomtype(x) != "polygons") {
			error("divide", "the geometry type must be polgyons")
		}
		if (is.list(n)) {
			vertical <- check_frac(n$v)
			horizontal <- check_frac(n$h)
			if (is.null(vertical) && is.null(horizontal)) return(x)
			out <- lapply(1:nrow(x), function(i) strip_polygon(x[i], vertical, horizontal))
		} else {
			n <- round(n)
			stopifnot(n > 0)
			if (n == 1) return(deepcopy(x))
			out <- lapply(1:nrow(x), function(i) divide_polygon(x[i], n, w, alpha, ...))
		}
		do.call(rbind, out)
	}
)

