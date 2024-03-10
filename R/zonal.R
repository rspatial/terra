
replace_with_label <- function(x, v, colnr) {
	ff <- is.factor(x)
	if (any(ff)) {
		cgs <- cats(x)
		for (f in which(ff)) {
			cg <- cgs[[f]]
			if (length(ff) == 1) {
				r <- 1:nrow(v)
			} else {
				r <- which(v[,1] == f)
			}
			i <- match(v[r,colnr], cg[,1])
			act <- activeCat(x, f) + 1
			if (!inherits(cg[[act]], "numeric")) {
				v[r, colnr] <- as.character(factor(cg[i, act], levels=unique(cg[[act]])))
			} else {
				v[r, colnr] <- cg[i, act]
			}
		}
	}
	v
}


setMethod("zonal", signature(x="SpatRaster", z="SpatRaster"),
	function(x, z, fun="mean", ..., w=NULL, wide=TRUE, as.raster=FALSE, filename="", overwrite=FALSE, wopt=list())  {

		group <- FALSE
		made_unique <- FALSE
		grast <- rast()
		nlz <- nlyr(z)
		znms <- names(z)
		if (nlz == 1) {
			group <- FALSE
		} else if ((!as.raster) && (nlz == 2)) {
			grast <- z[[2]]
			z <- z[[1]]
			group <- TRUE
		} else {
			ff <- is.factor(z)
			if (any(ff)) {
				levs <- levels(z)
				levels(z) <- NULL
			}
			z <- unique(z, as.raster=TRUE)
			made_unique <- TRUE
		}
		txtfun <- .makeTextFun(fun)
		if (inherits(txtfun, "character") && 
			(txtfun %in% c("max", "min", "mean", "sum", "notNA", "isNA"))) {

			if ((nlyr(z) > 1) && (nlyr(x) > 1)) {
				error("zonal", "x and z cannot both have more than one layer")
			}
			na.rm <- isTRUE(list(...)$na.rm)
			opt <- spatOptions()
			if (!is.null(w)) {
				if (txtfun != "mean") {
					error("zonal", "fun must be 'mean' when using weights")
				}
				sdf <- x@ptr$zonal_weighted(z@ptr, w@ptr, na.rm, opt)			
			} else {
				sdf <- x@ptr$zonal(z@ptr, grast@ptr, txtfun, na.rm, opt)
			}
			sdf <- messages(sdf, "zonal")
			out <- .getSpatDF(sdf)
			nz <- 1
			if (group) {
				out$layer <- out$layer + 1
				out <- replace_with_label(z, out, 2)
				out <- replace_with_label(grast, out, 3)
				if (nlyr(x) > 1) {
					out <- split(out[-1], out$layer)
					out <- Reduce(function(x, y) merge(x=x, y=y, by=1:2, all=TRUE), out)
					out <- out[order(out[,1], out[,2]), ]
				} else {
					out <- out[,-1]
				}
				colnames(out) <- c(znms, names(x))
				nz <- 2
			} else {
				if (made_unique && (!as.raster)) {
					ulevs <- cats(z)[[1]][, -c(1:2)]
					if (any(ff)) {
						for (f in which(ff)) {
							i <- match(ulevs[,f], levs[[f]][,1])
							ulevs[,f] <- levs[[f]][i,2]
						}
					}
					out <- cbind(ulevs, out[,-1,drop=FALSE])
					nz <- ncol(ulevs)
				} else {
					out <- replace_with_label(z, out, 1)
					colnames(out)[1] <- znms
				}
			}
			colnames(out) <- make.unique(colnames(out))
		} else {
			if (!is.null(w)) {
				error("zonal", "can only use weights when fun=mean")
			}
			compareGeom(x, z, lyrs=FALSE, crs=FALSE, ext=TRUE, rowcol=TRUE)
			#if (nlyr(z) > 1) {
			#	warn("zonal", "z can only have one layer with this function")
			#	z <- z[[1]]
			#}

			fun <- match.fun(fun)
			nl <- nlyr(x)
			nms <- names(x)
			if (group) {
				gzx <- c(grast, z, x)
				v <- as.data.frame(gzx, na.rm=FALSE)
				out <- stats::aggregate(v[,-c(1:2)], v[,1:2], fun, ...)
				colnames(out)[-c(1:2)] <- nms			
			} else {
				for (i in 1:nl) {
					xz <- c(x[[i]], z)
					v <- as.data.frame(xz, na.rm=FALSE)
					d <- stats::aggregate(v[,1], v[,2,drop=FALSE], fun, ...)
					colnames(d)[2] <- nms[i]
					if (i == 1) {
						out <- d
					} else {
						out <- merge(out, d, by=1)				
					}
				}
			}
		}		
		if (as.raster) {
			if (is.null(wopt$names) && (nlyr(x) == 1)) {
				wopt$names <- names(x)
			}
			levels(z) <- NULL
			out <- subst(z, out[,1], out[,-1], filename=filename, wopt=wopt)
		}		
		if (wide) {
			if (group) {
				nms <- names(out)
				isch <- inherits(out[,2], "character")
				#out <- stats::reshape(out, direction="wide", idvar=nms[c(1,3)], timevar=nms[2])
				out <- stats::reshape(out, direction="wide", idvar=nms[1], timevar=nms[2])
				if (isch) {
					colnames(out) <- gsub(paste0("^", nms[3], "."), "", colnames(out))
				}
				if (inherits(txtfun, "character") && (txtfun == "sum")) {
					out[is.na(out)] <- 0
				}
			}
		} else if (nz == 1){
			nls <- as.character(1:nlyr(x))
			colnames(out)[-1] <- nls
			if (colnames(out)[1] == "layer") colnames(out)[1] <- "zone"
			out <- stats::reshape(out, direction="long", varying=nls, timevar="layer",v.names="value")
			out <- out[, c(2,1,3)]
			rownames(out) <- NULL
		}
		out
	}
)


setMethod("zonal", signature(x="SpatRaster", z="SpatVector"),
	function(x, z, fun="mean", na.rm=FALSE, w=NULL, weights=FALSE, exact=FALSE, touches=FALSE, small=TRUE, as.raster=FALSE, as.polygons=FALSE, wide=TRUE, filename="", wopt=list())  {
		opt <- spatOptions()
		txtfun <- .makeTextFun(fun)
		if (!inherits(txtfun, "character")) {
			error("zonal", "this 'fun' is not supported. You can use extract instead")
		} else {
			if (txtfun == "table") {
				if (!is.null(w)) {
					error("cannot use 'w' when 'fun=table'")
				}
				v <- x@ptr$zonal_poly_table(z@ptr, weights[1], exact[1], touches[1], small[1], na.rm, opt)
				messages(x, "zonal")

				v <- lapply(v, function(i) if (length(i) == 0) NA else i)
				v <- lapply(1:length(v), function(i) cbind(i, matrix(v[[i]], ncol=2)))
				v <- do.call(rbind, v)
				v <- as.data.frame(v)
				colnames(v) <- c("zone", "value", "count")
				ff <- is.factor(x)[1]
				if (ff) {
					cg <- cats(x)[[1]]
					i <- match(v$value, cg[,1])
					act <- activeCat(x, 1) + 1
					v$value <- cg[i, act]
				}
				if (as.polygons | wide) {
					nms <- names(v)
					v <- stats::reshape(v, direction="wide", idvar=nms[1], timevar=nms[2])
					names(v) <- gsub("count.", "", names(v))
					v[is.na(v)] <- 0
					rownames(v) <- NULL
				}
				if (as.polygons) {
					values(z) <- v
					return(z)
				}
				return(v)
			} else {
				if (is.null(w)) {
					out <- x@ptr$zonal_poly(z@ptr, txtfun, weights[1], exact[1], touches[1], small[1], na.rm, opt)
				} else {
					if (txtfun != "mean") {
						error("zonal", "fun must be 'mean' when using weights")
					}
					out <- x@ptr$zonal_poly_weighted(z@ptr, w@ptr, weights[1], exact[1], touches[1], small[1], na.rm, opt)
				}
				messages(out, "zonal")
				out <- .getSpatDF(out)
			}
		}
		if (as.raster) {
			if (is.null(wopt$names)) {
				wopt$names <- names(x)
			}
			x <- rasterize(z, x, 1:nrow(z))
			subst(x, 1:nrow(out), out, filename=filename, wopt=wopt)
		} else if (as.polygons) {
			cbind(z, out)	
		} else {
			out
		}
	}
)


setMethod("zonal", signature(x="SpatVector", z="SpatVector"),
	function(x, z, fun=mean, ..., weighted=FALSE, as.polygons=FALSE)  {
		if (geomtype(z) != "polygons") {
			error("zonal", "x must be points, and z must be polygons")
		}
		if (nrow(x) == 0) {
			error("zonal", "x is empty")
		}
		isn <- which(sapply(values(x[1,]), is.numeric))
		if (!any(isn)) {
			error("zonal", "x has no numeric variables (attributes) to aggregate")
		}
		x <- x[,isn]
		if (geomtype(x) == "points") {
			r <- !relate(x, z, "disjoint", pairs=FALSE)
			i <- apply(r, 1, function(i) if(any(i)) which(i) else (NA))
			if (length(i) == 0) {
				error("zonal", "there are no points in x that overlap with the polygons in z")
			}
			a <- aggregate(values(x), data.frame(zone=i), fun, ...)
		} else {
			if (as.polygons) {
				zz <- z
				values(zz) <- data.frame(zone = 1:nrow(zz))
				i <- intersect(zz, x)
			} else {
				values(z) <- data.frame(zone = 1:nrow(z))
				i <- intersect(z, x)
			}
			if (nrow(i) == 0) {
				error("zonal", "the intersection of x and z is empty")
			}
			v <- values(i)
			if (weighted) {
				if (geomtype(i) == "lines") {
					v$w <- perim(i)
				} else {
					v$w <- expanse(i)
				}
				s <- split(v, v$zone)
				n <- ncol(v)-2
				s <- lapply(s, function(d) {
						out <- rep(NA, n)
						for (i in 2:n) {
							out[i-1] <- weighted.mean(d[[i]], w = d$w)
						}
						out
					})
				a <- data.frame(as.integer(names(s)), do.call(rbind, s))
				colnames(a) <- names(v)[-ncol(v)]
			} else {
				a <- aggregate(v[,-1,drop=FALSE], v[,1,drop=FALSE], fun, ...)
			}
		}
		if (as.polygons) {
			f <- basename(tempfile())
			z[[f]] <- 1:nrow(z)
			names(a)[1] = f
			a <- merge(z, a, by=f, all.x=TRUE)
			a[[f]] <- NULL
		}
		a
	}
)


setMethod("global", signature(x="SpatRaster"),
	function(x, fun="mean", weights=NULL, maxcell=Inf, ...)  {

		nms <- names(x)
		nms <- make.unique(nms)
		txtfun <- .makeTextFun(fun)

		opt <- spatOptions()
		if (!is.null(weights)) {
			stopifnot(inherits(weights, "SpatRaster"))
			stopifnot(txtfun %in% c("mean", "sum"))
			na.rm <- isTRUE(list(...)$na.rm)
			ptr <- x@ptr$global_weighted_mean(weights@ptr, txtfun, na.rm, opt)
			messages(ptr, "global")
			res <- .getSpatDF(ptr)
			rownames(res) <- nms
			return(res)
		}

		if (inherits(txtfun, "character")) {
			if (any(is.na(txtfun))) error("global", "fun cannot be NA")
			if (all(txtfun %in% c("prod", "max", "min", "mean", "sum", "range", "rms", "sd", "std", "sdpop", "notNA", "isNA"))) {
				txtfun[txtfun == "sdpop"] <- "std"
				i <- grep("range", txtfun)
				if (length(i) > 0) {
					txtfun <- txtfun[-i]
					txtfun <- c(txtfun, "min", "max")
				}
				txtfun <- unique(txtfun)
				na.rm <- isTRUE(list(...)$na.rm)
				#if (isTRUE(list(...)$old)) {
				#	ptr <- x@ptr$global(txtfun, na.rm, opt)			
				#} else {
				ptr <- x@ptr$mglobal(txtfun, na.rm, opt)
				#}
				messages(ptr, "global")
				res <- .getSpatDF(ptr)
				rownames(res) <- nms
				return(res)
			}
		}

		nl <- nlyr(x)
		res <- list()
		if (is.finite(maxcell)) {
			maxcell <- round(maxcell)
			if (maxcell < 1) error("global", "maxcell should be positive")
			x <- spatSample(x, maxcell, "regular", as.raster=TRUE)
		}
		for (i in 1:nl) {
			res[[i]] <- fun(values(x[[i]]), ...)
		}
		res <- do.call(rbind, res)
		res <- data.frame(res)

		# more efficient but more risky:
		#apply(data.frame(x), 2, fun, ...)

		if ((ncol(res) == 1) && (colnames(res) == "res")) {
			colnames(res) <- "global"
		}

		rownames(res) <- nms
		res
	}
)



setMethod("freq", signature(x="SpatRaster"),
	function(x, digits=0, value=NULL, bylayer=TRUE, usenames=FALSE, zones=NULL, wide=FALSE) {

		if (!is.null(zones)) {
			vna <- (!is.null(value) && is.na(value[1]))
#			if (vna) levels(x) <- NULL
			if (inherits(zones, "SpatVector")) {
				out <- vector("list", nrow(zones))
				for (i in 1:nrow(zones)) {
					z <- zones[i,]
					e <- align(ext(z), x, snap="near")
					if (!is.null(intersect(e, ext(x)))) {
						r <- crop(x, zones[i,], mask=TRUE, touches=FALSE)
						if (vna) {
							ra <- rasterize(zones[i,], r, NA, background=0, touches=FALSE)
							r <- cover(ra, r)
						}
						out[[i]] <- freq(r, digits=digits, value=value, bylayer=bylayer, usenames=usenames, zones=NULL)
						out[[i]]$zone <- i
					}
				}
			} else if (inherits(zones, "SpatRaster")) {
				compareGeom(x, zones, crs=FALSE)
				if (nlyr(zones) > 1) zones <- zones[[1]]
				u <- unlist(unique(zones))
				out <- vector("list", length(u))
				for (i in 1:length(u)) {
					r <- mask(x, zones, maskvalues=u[i], inverse=TRUE)
					out[[i]] <- freq(r, digits=digits, value=value, bylayer=bylayer, usenames=usenames, zones=NULL, wide=FALSE)
					out[[i]]$zone <- i
				}
			} else {
				error("freq", "zones must be a SpatVector or a SpatRaster")
			}
			out <- do.call(rbind, out)
			if (is.null(out)) return(out)
			out <- out[!is.na(out$count), ]
			if (nrow(out) == 0) return(out)
			out <- out[order(out$layer), ]
			if (wide) {
				out$count[is.na(out$count)] <- 0
				if (vna) {
					out$value <- "NA"
				}
				out <- stats::reshape(out, idvar=c("layer", "zone"), timevar="value", direction="wide")
				colnames(out) <- gsub("count.", "", colnames(out))
				out[is.na(out)] <- 0
			}
			return(out)
		}

		opt <- spatOptions()
		if (!bylayer) usenames <- FALSE

		if (!is.null(value)) {
			value <- unique(value)
			if (length(value) > 1) {
				error("freq", "value must have a length of one")
			}
			if (is.character(value)) {
				value <- value[value != ""]
				if (length(value) == 0) {
					error("freq", "no valid value")
				}
				ff <- is.factor(x)
				if (!any(ff)) {
					error("freq", "a character value is only meaningful for categorical rasters")
				}
				f <- freq(x[[ff]])
				if (usenames) {
					f$layer <- names(x)[f$layer]
				}
				f <- f[f$label == value,]
				return(f)
			}

			if (is.na(digits)) {
				v <- x@ptr$count(value, bylayer[1], FALSE, 0, opt)
			} else {
				v <- x@ptr$count(value, bylayer[1], TRUE, digits, opt)
				value <- round(value, digits)
			}
			if (bylayer) {
				v <- data.frame(layer=1:nlyr(x), value=value, count=v)
			} else {
				v <- data.frame(value=value, count=v)
			}

		} else {
			if (is.na(digits)) {
				v <- x@ptr$freq(bylayer[1], FALSE, 0, opt)
			} else {
				v <- x@ptr$freq(bylayer[1], TRUE, digits, opt)
			}
			v <- lapply(v, function(i) if (length(i) == 0) NA else i)

			v <- lapply(1:length(v), function(i) cbind(i, matrix(v[[i]], ncol=2)))
			v <- do.call(rbind, v)
			v <- as.data.frame(v)
			colnames(v) <- c("layer", "value", "count")
			ff <- is.factor(x)
			if (any(ff)) {
				cgs <- cats(x)
				v <- data.frame(v)
				for (f in which(ff)) {
					cg <- cgs[[f]]
					j <- which(v[,1] == f)
					i <- match(v[j,2], cg[,1])
					act <- activeCat(x, f) + 1
					if (!inherits(cg[[act]], "numeric")) {
						v[j, 2] <- as.character(factor(cg[i, act], levels=unique(cg[[act]])))
					} else {
						v[j, 2] <- cg[i, act]
					}
				}
			}
			if (!bylayer) {
#				if (nlyr(x) > 1)
#					v <- aggregate(v[,"count",drop=FALSE], v[,"value", drop=FALSE], sum)
#				} 
				v <- v[,-1]
			}
		}
		if (usenames) {
			v$layer <- names(x)[v$layer]
		}
		if (wide) {
			v$count[is.na(v$count)] <- 0
			if ((!is.null(value)) && is.na(value)) {
				v$value <- "NA"
			}
			v <- stats::reshape(v, idvar="layer", timevar="value", direction="wide")
			colnames(v) <- gsub("count.", "", colnames(v))
			v[is.na(v)] <- 0
		}
	
		v
	}
)





setMethod ("expanse", "SpatRaster",
	function(x, unit="m", transform=TRUE, byValue=FALSE, zones=NULL, wide=FALSE, usenames=FALSE) {
		opt <- spatOptions()
		if (!is.null(zones)) {
			if (!inherits(zones, "SpatRaster")) {
				error("expanse", "zones must be a SpatRaster")
			}
			compareGeom(x, zones, lyrs=FALSE, crs=FALSE, ext=TRUE, rowcol=TRUE)
			v <- x@ptr$sum_area_group(zones@ptr, unit[1], transform[1], byValue[1], opt)
			messages(x)
			v <- lapply(v, function(i) matrix(i, ncol=4, byrow=TRUE))
			v <- data.frame(do.call(rbind, v))
			colnames(v) <- c("layer", "value", "zone", "area")
			v[,1] <- v[,1] + 1
			if (byValue) {
				v <- replace_with_label(x, v, 2)	
				v <- replace_with_label(zones, v, 3)	
			} else {
				v <- replace_with_label(zones, v, 3)	
				v$value <- NULL
			}
			if (wide) {
				if (byValue) {
					v <- stats::reshape(v, idvar=c("layer", "zone"), timevar="value", direction="wide")
					colnames(v) <- gsub("area.", "", colnames(v))
				} else {
					v <- stats::reshape(v, idvar=c("layer"), timevar="zone", direction="wide")
					colnames(v) <- gsub("area.", "", colnames(v))
				}
				v[is.na(v)] <- 0
			}
		} else {
			v <- x@ptr$sum_area(unit, isTRUE(transform[1]), isTRUE(byValue[1]), opt)
			x <- messages(x, "expanse")
			if (byValue) {
				v <- lapply(1:length(v), function(i) cbind(i, matrix(v[[i]], ncol=2, byrow=TRUE)))
				v <- data.frame(do.call(rbind, v))
				colnames(v) <- c("layer", "value", "area")
				v <- replace_with_label(x, v, 2)	
			} else {
				v <- v[[1]]
				v <- data.frame(layer=1:length(v), area=v)
			}
			if (wide) {
				if (byValue) {
					v <- stats::reshape(v, idvar="layer", timevar="value", direction="wide")
					colnames(v) <- gsub("area.", "", colnames(v))
				}
				v[is.na(v)] <- 0
			}
		}
		if (usenames) {
			v$layer <- names(x)[v$layer]
		}
		v
	}
)


