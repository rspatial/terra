# these methods require the dev version of leaflet

checkLeafLetVersion <- function() {
	v <- utils::packageVersion("leaflet")
	if (v < "2.1.2.9000") {
		error("plet", "plet needs the development version of leaflet")
	}
}

popUp <- function(x) {
	nms <- names(x)
	if (length(nms) > 0) {
		s <- sapply(1:length(nms), function(i) paste0(nms[i], ": ", x[[i, drop=TRUE]]))
		if (is.null(dim(s))) {
			paste(s, collapse="<br>")
		} else {
			apply(s, 1, function(i) paste(i, collapse="<br>"))		
		}
	} else {
		paste("geom", 1:nrow(x), sep="_")
	}
}

makelonlat <- function(x) {
	geo <- is.lonlat(x)
	if (is.na(geo)) {
		geo <- is.lonlat(x, TRUE, TRUE)
		if (geo) return(x)
		error("plet", "coordinate reference system is unknown and does not look like lon/lat")
	} else if (!geo) {
		project(x, "+proj=longlat")
	} else {
		x
	}
}

setMethod("plet", signature(x="missing"),
	function(x) {
		leaflet::leaflet()
	}
)

#.leaflet-container {
#    background: #FFF;
#}


baselayers <- function(tiles, wrap=TRUE) {
	map <- leaflet::leaflet()
	if ((!is.null(tiles)) && (length(tiles) > 0)) {
		if ("Streets" %in% tiles) {
			map <- leaflet::addTiles(map, group="Streets", 
			       options=leaflet::tileOptions(noWrap = !wrap))
		} 
		tiles <- tiles[tiles != "Streets"]
		if (length(tiles) > 0) {
			for (i in 1:length(tiles)) {
				map <- leaflet::addProviderTiles(map, tiles[i], group=tiles[i], 
					       options=leaflet::tileOptions(noWrap = !wrap))
			}
		}
	}
	map
}


.get_leg <- function(v, type="", dig.lab=3, cols, breaks=NULL, breakby="eqint", sort=TRUE, decreasing=FALSE,  ...) {
	out <- list(v=v, leg=list())
	
	if (is.null(type)) type <- ""
	if (type == "continuous") type <- "interval"	
	if ((!is.numeric(v)) || (length(unique(v)) < 11)) {
		type <- "classes"
	} else if (type == "") {
		type <- "interval"
	} else {
		type <- match.arg(type, c("interval", "classes"))
	}
	out$legend_type <- type

	out$uv <- unique(v)
	
	if (inherits(cols, "function")) {
		cols <- cols(100)
	}
	out$cols <- cols
		
	if (out$legend_type == "classes") {
		out$legend_sort <- sort[1]
		out$legend_sort_decreasing <- decreasing[1]
		out <- .vect.legend.classes(out)
	} else if (out$legend_type == "interval") {
		out$breaks <- breaks
		out$breakby <- breakby
		out$range <- range(v, na.rm=TRUE)
		out <- .vect.legend.interval(out, dig.lab=dig.lab)
	}

	out
}


setMethod("plet", signature(x="SpatVector"),
	function(x, y="", col, fill=0.2, main=y, cex=1, lwd=2, border="black", alpha=1, popup=TRUE, label=FALSE, split=FALSE, tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), wrap=TRUE, legend="bottomright", collapse=FALSE, type=NULL, breaks=NULL, breakby="eqint", sort=TRUE, decreasing=FALSE, map=NULL, ...)  {

		#checkLeafLetVersion()
		y <- unique(y)
		if (length(y) > 1) {
			y = y[1]
#			xvc <- svc(lapply(y, function(i) x[,i]))
#			if (is.numeric(y)) {
#				names(xvc) <- names(x)[y]
#			} else {
#				names(xvc) <- y			
#			}
#			plet(xvc, col=col, fill=0.2, alpha=alpha, cex=cex, lwd=lwd, border=border, popup=popup, label=label, split=split, tiles=tiles, wrap=wrap, legend=legend, collapse=collapse,  map=map, ...)	
#type=type, breaks=breaks, breakby=breakby, sort=sort, decreasing=decreasing,
		}
		
		if (missing(col)) col <- grDevices::rainbow
		alpha <- max(0, min(1, alpha))
		fill <- max(0, min(1, fill))
		x <- makelonlat(x)

		if (is.null(map)) {
			tiles <- unique(as.character(tiles))
			tiles <- tiles[tiles!=""]
			map <- baselayers(tiles, wrap)
		} else {
			tiles <- NULL
		}
		g <- geomtype(x)
		leg <- NULL
		if (y == "") { # no legend
			group <- x@cpp$layer
			if (group == "") group = g
			cols <- .getCols(nrow(x), col)
			pop  <- lab <- NULL
			if (isTRUE(popup[1])) pop <- popUp(x)
			if (isTRUE(label[1])) lab <- 1:nrow(x)
			if (g == "polygons") {
				map <- leaflet::addPolygons(map, data=x, label=lab, group=group, 
							fillColor=cols, fillOpacity=fill, opacity=alpha, 
							popup=pop, color=border, weight=lwd, ...)
			} else if (g == "lines") {
				map <- leaflet::addPolylines(map, data=x, label=lab, group=group,  
							color=cols, opacity=alpha,  popup=pop, weight=lwd, ...)
			} else {
				map <- leaflet::addCircleMarkers(map, data=x, radius=cex, popup=pop, group=group, 
							label=lab, opacity=alpha, color=cols, ...)
			}

			if (length(tiles) > 1) {
				map <- leaflet::addLayersControl(map, baseGroups = tiles, overlayGroups=group,
						options = leaflet::layersControlOptions(collapsed=collapse))
			} else {
				map <- leaflet::addLayersControl(map, overlayGroups = group, 
						options = leaflet::layersControlOptions(collapsed=collapse))
			}			
			map
		} else { # legend
			y <- y[1]
			if (is.numeric(y)) {
				y <- round(y)
				stopifnot((y > 0) && (y <= nlyr(x)))
				y <- names(x)[y]
			}
			stopifnot(y %in% names(x))
			x <- x[, y]
			v <- values(x)[,1]

			if (split) { 
				u <- unique(v)
				cols <- .getCols(length(u), col)
				for (i in seq_along(u)) {
					s <- x[v == u[i], ]
					pop <- lab <- NULL
					if (isTRUE(popup[1])) pop <- popUp(s)
					if (isTRUE(label[1])) lab <- u
					if (g == "polygons") {
						map <- leaflet::addPolygons(map, data=s, label=lab[i], group=u[i], 
							fillColor=cols[i],  fillOpacity=fill, opacity=alpha, popup=pop, 
							col=border, ...)
					} else if (g == "lines") {
						map <- leaflet::addPolylines(map, data=s, label=lab[i], group=u[i], 
							col=cols[i], opacity=alpha, popup=pop, ...)
					} else {
						map <- leaflet::addCircleMarkers(map, data=s, label=lab[i], group=u[i], 
							col=cols[i], fillOpacity=fill, opacity=alpha, popup=pop, radius=cex, ...)
					}
				}
				if (length(tiles) > 1) {
					map <- leaflet::addLayersControl(map, baseGroups = tiles, overlayGroups = u, 
						options = leaflet::layersControlOptions(collapsed=collapse))
				} else {
					map <- leaflet::addLayersControl(map, overlayGroups = u, 
						options = leaflet::layersControlOptions(collapsed=collapse))
				}
			} else { # do not split
				#vcols <- cols[1:length(v)]
				leg <- .get_leg(v, type=type, dig.lab=3, cols=col, breaks=breaks, breakby=breakby, sort=sort, decreasing=decreasing, ...)
				pop <- lab <- NULL
				if (isTRUE(popup[1])) pop <- popUp(x)
				if (isTRUE(label[1])) lab <- v
				if (g == "polygons") {
					map <- leaflet::addPolygons(map, data=x, label=lab, group=y, 
						fillColor=leg$main_cols, opacity=alpha, fillOpacity=fill, 
						col = border, popup=pop, ...)
				} else if (g == "lines") {
					map <- leaflet::addPolylines(map, data=x, label=lab, group=y, 
						col=leg$main_cols, popup=pop, opacity=alpha, ...)
				} else {
					map <- leaflet::addCircleMarkers(map, data=x, label=lab, group=y, 
						col=leg$main_cols, radius=cex, popup=pop, fillOpacity=fill, opacity=alpha, ...)
				}
				if (length(tiles) > 1) {
					map <- leaflet::addLayersControl(map, baseGroups = tiles, overlayGroups=y,
						options = leaflet::layersControlOptions(collapsed=collapse))
				} else {
					map <- leaflet::addLayersControl(map, overlayGroups = y, 
						options = leaflet::layersControlOptions(collapsed=collapse))
				}

			}
			if ((!is.null(legend)) && (!is.null(leg))) {
				if (leg$legend_type != "") {
					main <- gsub("\n", "</br>", main[1])
					op = ifelse(g == "polygons", fill, 1)
					map <- leaflet::addLegend(map, position=legend, colors=leg$leg$fill, 
							labels=as.character(leg$leg$legend), opacity=op, title=main)
				}
			}
			map
		}
	}
)



setMethod("plet", signature(x="SpatVectorCollection"),
	function(x, col, fill=0, cex=1, lwd=2, border="black", alpha=1, popup=TRUE, label=FALSE, tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), wrap=TRUE, legend="bottomright", collapse=FALSE, map=NULL)  {

		#checkLeafLetVersion()

		if (is.null(map)) {
			tiles <- unique(as.character(tiles))
			tiles <- tiles[tiles!=""]
			map <- baselayers(tiles, wrap)
		} else {
			tiles <- NULL
		}
		nms <- names(x)
		nms[nchar(nms) == 0] <- "X"
		nms <- make.unique(nms)

		if (missing(col)) col <- grDevices::rainbow
		n <- length(x)
		if (is.function(col)) {
			cols <- col(n)
		} else {
			cols <- rep_len(col, n) 
		}

		lwd <- rep_len(lwd, n) 

		alpha <- rep_len(alpha, n) 
		alpha <- pmax(0, min(1, alpha))
		fill <- rep_len(fill, n) 
		fill <- pmax(0, min(1, fill))
		popup <- rep_len(popup, n) 
		label <- rep_len(label, n) 
		border <- rep_len(border, n) 

		for (i in 1:n) {
			v <- x[i]
			v <- makelonlat(v)
			g <- geomtype(v)
			pop <- NULL
			lab <- NULL
			if (popup[i]) {
				pop <- popUp(v)
			} 
			if (label[i]) {
				lab <- 1:nrow(v)
			}
			if (g == "polygons") {
				map <- leaflet::addPolygons(map, data=v, weight=lwd[i], fillColor=cols[i], 
				fillOpacity=fill[i], col=border[i], opacity=alpha[i], popup=pop, 
				label=lab, group=nms[i])
			} else if (g == "lines") {
				map <- leaflet::addPolylines(map, data=v, weight=lwd[i], opacity=alpha[i], 
				col=cols[i], group=nms[i], popup=pop, label=lab)
			} else {
				map <- leaflet::addCircleMarkers(map, data=v, radius=cex[i], popup=pop, label=lab, opacity=alpha[i], col=cols[i], group=nms[i])
			}
		}
		if (length(tiles) > 1) {
			map <- leaflet::addLayersControl(map, baseGroups = tiles, overlayGroups = nms, 
				options = leaflet::layersControlOptions(collapsed=collapse))
		} else {
			map <- leaflet::addLayersControl(map, overlayGroups = nms, 
				options = leaflet::layersControlOptions(collapsed=collapse))
		}
		map
	}
)


setMethod("lines", signature(x="leaflet"),
	function(x, y, col, lwd=2, alpha=1, ...)  {
		if (inherits(y, "SpatVector")) {
			if (nrow(y) == 0) return(x)
			y <- makelonlat(y)
			if (missing(col)) col <- "black"
			if (!(geomtype(y) %in% c("lines", "polygons"))) {
				error("lines", "SpatVector y must have either lines or polygons geometry")
			}
			leaflet::addPolylines(x, data=y, weight=lwd, opacity=alpha, col=col, ...)
		} else if (inherits(y, "SpatVectorCollection")) {
			nms <- names(y)
			n <- length(y)
			nms[nchar(nms) == 0] <- "X"
			nms <- make.unique(nms)
			if (is.function(col)) {
				cols <- col(n)
			} else {
				cols <- rep_len(col, n) 
			}
			lwd <- rep_len(lwd, n) 
			alpha <- rep_len(alpha, n) 
			for (i in 1:length(nms)) {
				x <- leaflet::addPolylines(x, data=y[i], weight=lwd[i], opacity=alpha[i], col=cols[i], group=nms[i], ...)
			}
			collapse=FALSE
			leaflet::addLayersControl(x, overlayGroups = nms, options = leaflet::layersControlOptions(collapsed=collapse))
		}
	}
)

setMethod("points", signature(x="leaflet"),
	function(x, y, col, cex=1, alpha=1, popup=FALSE, ...)  {
		stopifnot(inherits(y, "SpatVector"))
		if (nrow(y) == 0) return(x)
		y <- makelonlat(y)
		if (missing(col)) col <- "black"
		if (!(geomtype(y) == "points")) {
			if (geomtype(y) == "polygons") {
				y <- centroids(y)
			} else {
				y <- as.points(y)
			}
		}
		if (popup) {
			popup=popUp(y)
		} else {
			popup <- NULL
		}
		leaflet::addCircleMarkers(x, data=y, radius=cex, popup=popup, label=1:nrow(y), opacity=alpha, col=col, ...)
	}
)




make.panel <- function(x, maxcell) {

	nl <- nlyr(x)
	x <- spatSample(x, maxcell/nl, "regular", as.raster=TRUE, warn=FALSE)
	if (is.lonlat(x)) {
		asp <- 1/cos((mean(ext(x)[3:4]) * pi)/180)
	} else {
		asp <- 1
	} 
	crs(x, warn=FALSE) <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
	ext(x) = c(0,1,0,asp)
	#if (!is.null(add)) {
	#	e <- as.vector(ext(x))
	#	for (i in 1:length(add)) {
	#		v <- project(v, crs(x))
	#		add[i] <- rescale(v, fx=1/diff(e[1:2]), fy=1/diff(e[3:4]))
	#	}
	#}

	asp <- asp * nrow(x) / ncol(x)
	nc <- ceiling(2*sqrt(nl/2) * asp)
	nr <- ceiling(nl / nc)
	nc <- ceiling(nl / nr)

	e <- as.vector(ext(x))
	r <- res(x)
	skiprow <- -asp - max(1, min(10, trunc(nrow(x)/20))) *r[2]
	skipcol <- 1 + max(1, min(10, trunc(ncol(x)/20))) *r[1]
#	skipcol <- 1 + max(r[1], min(10, r[1] * trunc(nrow(x)/20)))
	labs <- data.frame(x=0, y=0, label=names(x))
	rw = cl = 0
	y <- vector(mode="list", length=nl)
	off <- 0 #.4 / nr
	for (i in 1:nl) {
		y[[i]] <- shift(x[[i]], cl * skipcol, rw * skiprow)
		e <- as.vector(ext(y[[i]]))
		labs[i,1:2] <- c(mean(e[1:2]), e[4]-off)
		cl <- cl + 1
		if (cl == nc) {
			cl <- 0
			rw <- rw + 1
		}
	}
	labs <- vect(labs, geom=c("x", "y"), crs=crs(x))
	labs <- project(labs, "+proj=longlat")
	x <- merge(sprc(y))
	list(x, labs)
}


setMethod("plet", signature(x="SpatRaster"),
	function(x, y=1, col, alpha=0.8, main=names(x), tiles=NULL, wrap=TRUE, maxcell=500000, legend="bottomright", shared=FALSE, panel=FALSE, collapse=TRUE, map=NULL)  {

		#checkLeafLetVersion()

		if (is.na(crs(x))) {
			error("plet", "x must have a crs")
		}


#		if (!is.null(add)) {
#			if (inherits(add, "SpatVector")) {
#				add <- svc(makelonlat(add))
#			} else if (inherits(add, "SpatVectorCollection")) {
#				for (i in 1:length(add)) {
#					add[i] <- makelonlat(add[i])
#				}
#			} else {
#				error("plet", "add should be a SpatVector or a SpatVectorCollection")
#			}
#		}

		alpha <- max(0, min(1, alpha))
		
		e <- ext(x)
		if (is.lonlat(x) && ((e$ymin < -85) || (e$ymax > 85))) {
			yr1 <- e$ymax - e$ymin
			e$ymin <- max(e$ymin, -85)
			e$ymax <- min(e$ymax, 85)
			yr2 <- e$ymax - e$ymin
			x <- spatSample(x[[y]], (yr1/yr2) * maxcell, "regular", as.raster=TRUE, warn=FALSE)
			x <- crop(x, e) 
		}
		
		if (panel) {
			tiles <- NULL
			p <- make.panel(x, maxcell) #, add)
			x <- p[[1]]
			p <- p[[2]]
			add <- p[[3]]
			main <- ""
		} else {
			x <- spatSample(x[[y]], maxcell, "regular", as.raster=TRUE, warn=FALSE)
		}

		if (is.null(map)) {
			tiles <- unique(as.character(tiles))
			tiles <- tiles[tiles!=""]
			if (length(tiles) > 1) {
				tiles <- tiles[1]
				warn("plet", "only a single tileset can be used with raster data")
			}
			map <- baselayers(tiles, wrap)
		} else {
			tiles <- NULL
		}
		if (missing(col)) {
			col <- rev(grDevices::terrain.colors(255))
		}

		main <- gsub("\n", "</br>", main)
		if (length(main) != length(y)) {
			main <- rep_len(main, length(x))[y]
		}


		notmerc <- isTRUE(crs(x, describe=TRUE)$code != "3857")

		if (nlyr(x) == 1) {
			map <- leaflet::addRasterImage(map, x, colors=col, opacity=alpha, project=notmerc)
			if (!is.null(legend)) {
				if (!all(hasMinMax(x))) setMinMax(x)
				r <- minmax(x)
				v <- seq(r[1], r[2], length.out=5)
				pal <- leaflet::colorNumeric(col, v, reverse = TRUE)
				map <- leaflet::addLegend(map, legend, pal=pal, values=v, opacity=1, title=main[1],
					  labFormat = leaflet::labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
			}
			if (panel) {
				#map <- leaflet::addCircleMarkers(map, data=p, label=p$label, radius=1, opacity=1, col="red")
				map <- leaflet::addLabelOnlyMarkers(map, label=p$label, data=p, 
                      labelOptions = leaflet::labelOptions(noHide = T, textOnly = T))
			}
		} else {
			nms <- make.unique(names(x))
			many_legends <- one_legend <- FALSE
			if (!is.null(legend)) {
				if (!all(hasMinMax(x))) setMinMax(x)
				r <- minmax(x)
				if (shared) {
					rr <- range(r)
					pal <- leaflet::colorNumeric(col, rr, na.color="#00000000")
					one_legend <- TRUE
				} else {
					many_legends <- TRUE
				}
			} else {
				one_legend <- FALSE
			}
			for (i in 1:nlyr(x)) {
				if (one_legend) {
					map <- leaflet::addRasterImage(map, x[[i]], colors=pal, opacity=alpha, group=nms[i], project=notmerc)
				} else {
					map <- leaflet::addRasterImage(map, x[[i]], colors=col, opacity=alpha, group=nms[i], project=notmerc)
					if (many_legends) {
						v <- seq(r[1,i], r[2,i], length.out=5)
						pal <- leaflet::colorNumeric(col, v, reverse=TRUE)
						map <- leaflet::addLegend(map, position=legend, pal=pal, values=v, 
							  title=main[i], opacity=1, group=nms[i],
							  labFormat = leaflet::labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
					}
				}
			}
			map <- leaflet::addLayersControl(map, baseGroups=nms,
				options = leaflet::layersControlOptions(collapsed=collapse))
			if (many_legends) {	# show one legend at a time
				map <- htmlwidgets::onRender(map, 
					"function(el, x) {
						var updateLegend = function () {
						var selectedGroup = document.querySelectorAll('input:checked')[0].nextSibling.innerText.substr(1);
						document.querySelectorAll('.legend').forEach(a => a.hidden=true);
						document.querySelectorAll('.legend').forEach(l => {
							if (l.children[0].children[0].innerText == selectedGroup) l.hidden=false;
						});
					};
					updateLegend();
					this.on('baselayerchange', e => updateLegend());}")
			} else if (one_legend) {
				v <- seq(rr[1], rr[2], length.out=5)
				pal <- leaflet::colorNumeric(col, v, reverse = TRUE)
				map <- leaflet::addLegend(map, position=legend, pal=pal, values=v, opacity=1, group=nms[i],
						  labFormat = leaflet::labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
			}
		}
		map
	}
)


