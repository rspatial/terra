
# these methods require the dev version of leaflet


popUp <- function(x) {
	nms <- names(x)
	s <- sapply(1:length(nms), function(i) paste0(nms[i], ": ", x[[i, drop=TRUE]]))
	apply(s, 1, function(i) paste(i, collapse="<br>"))
}


setMethod("plet", signature(x="SpatVector"),
	function(x, y="", col, split=FALSE, tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), alpha=1, legend="bottomright", collapse=FALSE, cex=1, map=NULL)  {
	
		if (missing(col)) col = grDevices::rainbow
		
		alpha <- 1 - max(0, min(1, alpha))
		#stopifnot(packageVersion("leaflet") > "2.1.1")
		if (is.null(map)) {
			map <- leaflet::leaflet()
		} else {
			tiles <- ""
		}
		g <- geomtype(x)
		if (!all(tiles == "")) {
			if ("Streets" %in% tiles) {
				map <- leaflet::addTiles(map, group="Streets")			
			} 
			tiles2 <- tiles[tiles != "Streets"]
			if (length(tiles) > 0) {
				tiles2 <- tiles[tiles != "Streets"]
				for (i in 1:length(tiles2)) {
					map <- leaflet::addProviderTiles(map, tiles2[i], group=tiles2[i])
				}
			}
		}
		y <- y[1]
		if (y == "") {
			cols <- .getCols(nrow(x), col)
			if (g == "polygons") {
				map <- leaflet::addPolygons(map, data=x, label=1:nrow(x),  
							col=cols, fillOpacity=alpha, popup=popUp(x))
			} else if (g == "lines") {
				map <- leaflet::addPolylines(map, data=x, label=1:nrow(x),  
							col=cols, fillOpacity=alpha, popup=popUp(x))
			} else {
				map <- leaflet::addMarkers(map, data=x, label=1:nrow(x),  
							col=cols, fillOpacity=alpha, popup=popUp(x))			
			}
			if (!all(tiles == "")) {
				map <- leaflet::addLayersControl(map, baseGroups = tiles, 
						options = leaflet::layersControlOptions(collapsed=collapse))
			}
			map
		} else {
			if (is.numeric(y)) {
				y <- round(y)
				stopifnot((y > 0) && (y <= nlyr(x)))
				y <- names(x)[y]
			}
			stopifnot(y %in% names(x))
			u <- unique(x[[y, drop=TRUE]])
			cols <- .getCols(length(u), col)
			if (split) {
				for (i in seq_along(u)) {
					s <- x[x[[y]] == u[i], ]
					if (g == "polygons") {
						map <- leaflet::addPolygons(map, data=s, label=u[i], group=u[i], 
							col=cols[i], fillOpacity=alpha, popup=popUp(s))
					} else if (g == "lines") {
						map <- leaflet::addPolylines(map, data=s, label=u[i], group=u[i], 
							col=cols[i], fillOpacity=alpha, popup=popUp(s))
					} else {
						map <- leaflet::addCircleMarkers(map, data=s, label=u[i], group=u[i], 
							col=cols[i], fillOpacity=alpha, popup=popUp(s))					
					}
				}
				if (all(tiles == "")) {
					map <- leaflet::addLayersControl(map, overlayGroups = u, 
						options = leaflet::layersControlOptions(collapsed=collapse))
				} else {
					map <- leaflet::addLayersControl(map, baseGroups = tiles, overlayGroups = u, 
						options = leaflet::layersControlOptions(collapsed=collapse))
				}
			} else {
				values <- x[[y,drop=TRUE]]
				vcols <- cols[as.numeric(as.factor(values))]
				if (g == "polygons") {
					map <- leaflet::addPolygons(map, data=x, label=values,  
						col=vcols, fillOpacity=alpha, popup=popUp(x))
				} else if (g == "lines") {
					map <- leaflet::addPolylines(map, data=x, label=values,  
						col=vcols, popup=popUp(x), fillOpacity=alpha)
				} else {
					map <- leaflet::addCircleMarkers(map, data=x, label=values,  
						col=vcols, radius=cex, popup=popUp(x), fillOpacity=alpha)
				}
				if (!all(tiles == "")) {
					map <- leaflet::addLayersControl(map, baseGroups = tiles, 
						options = leaflet::layersControlOptions(collapsed=collapse))
				}

			}
			if (!is.null(legend)) {
				map <- leaflet::addLegend(map, position=legend, colors=cols, labels=u, opacity=1, title=y)
			}
			map
		}
	}
)



setMethod("lines", signature(x="leaflet"),
	function(x, y, col, lwd=3, alpha=1)  {
		stopifnot(inherits(y, "SpatVector"))
		if (nrow(y) == 0) return(x)
		if (missing(col)) col <- "black"
		if (!(geomtype(y) %in% c("lines", "polygons"))) {
			error("lines", "SpatVector y must have either lines or polygons geometry")
		}
		leaflet::addPolylines(x, data=y, weight=lwd, opacity=alpha, col=col)
	}
)

setMethod("points", signature(x="leaflet"),
	function(x, y, col, cex=1, alpha=1, popup=FALSE)  {
		stopifnot(inherits(y, "SpatVector"))
		if (nrow(y) == 0) return(x)
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
		leaflet::addCircleMarkers(x, data=y, radius=cex, popup=popup, label=1:nrow(y), opacity=alpha, col=col)
	}
)




setMethod("plet", signature(x="SpatRaster"),
	function(x, y=1, col, alpha=0.2, tiles="", maxcell=500000, legend="bottomright", map=NULL, ...)  {
		#stopifnot(packageVersion("leaflet") > "2.1.1")

		alpha <- 1 - max(0, min(1, alpha))

		if (is.null(map)) {
			map <- leaflet::leaflet()
		} else {
			tiles <- ""
		}
		if (!all(tiles == "")) {
			if ("Streets" %in% tiles) {
				map <- leaflet::addTiles(map, group="Streets")			
			} 
			tiles2 <- tiles[tiles != "Streets"]
			if (length(tiles) > 0) {
				for (i in 1:length(tiles2)) {
					map <- leaflet::addProviderTiles(map, tiles2[i], group=tiles2[i])
				}
			}
		}
		col <- rev(grDevices::terrain.colors(255))
		
		y <- y[1]
		x <- spatSample(x[[y]], maxcell, "regular", as.raster=TRUE)
	
		map <- leaflet::addRasterImage(map, x, colors = col, opacity=opacity)
		if (!is.null(legend)) {
			r <- minmax(x)
			v <- seq(r[1], r[2], 5)
			pal <- leaflet::colorNumeric(col, v, reverse = TRUE)
			map <- leaflet::addLegend(map, legend, pal=pal, values=v, 
                  labFormat = leaflet::labelFormat(transform = function(x) sort(x, decreasing = TRUE)))	
		}
		map
	}
)


