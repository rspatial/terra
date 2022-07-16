
# these methods require the dev version of leaflet


popUp <- function(x) {
	nms <- names(x)
	s <- sapply(1:length(nms), function(i) paste0(nms[i], ": ", x[[i, drop=TRUE]]))
	apply(s, 1, function(i) paste(i, collapse="<br>"))
}


setMethod("plotlet", signature(x="SpatVector"),
	function(x, y="", col=rainbow, split=FALSE, tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), map=NULL, fillOpacity=0, legend="bottomright", collapse=FALSE, ...)  {
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
		if (y == "") {
			cols <- .getCols(nrow(x), col)
			if (g == "polygons") {
				map <- leaflet::addPolygons(map, data=x, label=1:nrow(x),  
							col=cols, fillOpacity=fillOpacity, popup=popUp(x))
			} else if (g == "lines") {
				map <- leaflet::addPolylines(map, data=x, label=1:nrow(x),  
							col=cols, fillOpacity=fillOpacity, popup=popUp(x))
			} else {
				map <- leaflet::addMarkers(map, data=x, label=1:nrow(x),  
							col=cols, fillOpacity=fillOpacity, popup=popUp(x))			
			}
			if (!all(tiles == "")) {
				map <- leaflet::addLayersControl(map, baseGroups = tiles, 
						options = leaflet::layersControlOptions(collapsed=collapse))
			}
			map
		} else {
			stopifnot(y %in% names(x))
			u <- unique(x[[y, drop=TRUE]])
			cols <- .getCols(length(u), col)
			if (split) {
				for (i in seq_along(u)) {
					s <- x[x[[y]] == u[i], ]
					if (g == "polygons") {
						map <- leaflet::addPolygons(map, data=s, label=u[i], group=u[i], 
							col=cols[i], fillOpacity=fillOpacity, popup=popUp(s))
					} else if (g == "lines") {
						map <- leaflet::addPolylines(map, data=s, label=u[i], group=u[i], 
							col=cols[i], fillOpacity=fillOpacity, popup=popUp(s))
					} else {
						map <- leaflet::addCircleMarkers(map, data=s, label=u[i], group=u[i], 
							col=cols[i], fillOpacity=fillOpacity, popup=popUp(s))					
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
						col=vcols, fillOpacity=fillOpacity, popup=popUp(x))
				} else if (g == "lines") {
					map <- leaflet::addPolylines(map, data=x, label=values,  
						col=vcols, fillOpacity=fillOpacity, popup=popUp(x))
				} else {
					map <- leaflet::addMarkers(map, data=x, label=values,  
						col=vcols, fillOpacity=fillOpacity, popup=popUp(x))
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



setMethod("plotlet", signature(x="SpatRaster"),
	function(x, y="", col="Spectral", opacity=0.8, tiles="", map=NULL, ...)  {
		#stopifnot(packageVersion("leaflet") > "2.1.1")
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
		leaflet::addRasterImage(map, x, colors = col, opacity=opacity)
	}
)

