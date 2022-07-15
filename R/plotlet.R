
# these methods require the dev version of leaflet


popUp <- function(x) {
	nms <- names(x)
	s <- sapply(1:length(nms), \(i) paste0(nms[i], ": ", x[[i, drop=TRUE]]))
	apply(s, 1, function(i) paste(i, collapse="<br>"))
}


setMethod("plotlet", signature(x="SpatVector"),
	function(x, y="", col=rainbow, split=TRUE, tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), m=NULL ...)  {
		#stopifnot(packageVersion("leaflet") > "2.1.1")
		if (is.null(m)) {
			m <- leaflet::leaflet()
		} else {
			tiles <- ""
		}
		g <- geomtype(x)
		if (!all(tiles == "")) {
			if ("Streets" %in% tiles) {
				m <- leaflet::addTiles(m, group="Streets")			
			} 
			tiles2 <- tiles[tiles != "Streets"]
			if (length(tiles) > 0) {
				tiles2 <- tiles[tiles != "Streets"]
				for (i in 1:length(tiles2)) {
					m <- leaflet::addProviderTiles(m, tiles2[i], group=tiles2[i])
				}
			}
		}
		if (y == "") {
			if (g == "polygons") {
				m <- leaflet::addPolygons(m, x)
			} else if (g == "lines") {
				m <- leaflet::addPolylines(m, x)
			} else {
				m <- leaflet::addMarkers(m, x)			
			}
		} else {
			stopifnot(y %in% names(x))
			if (split) {
				u <- unique(x[[y, drop=TRUE]])
				cols <- .getCols(length(u), col)
				for (i in seq_along(u)) {
					s <- x[x[[y]] == u[i], ]
					m <- leaflet::addPolygons(m, data=s, label=u[i], group=u[i], 
						col=cols[i], fillOpacity=0, popup=popUp(s))
				}
				if (!all(tiles == "")) {
					leaflet::addLayersControl(m,
						baseGroups = tiles,
						overlayGroups = u, 
						options = layersControlOptions(collapsed = FALSE))
				} else {
					leaflet::addLayersControl(m,
						overlayGroups = u, 
						options = layersControlOptions(collapsed = FALSE))
				}
			} else {
				stop("tbd")
			}
		}
	}
)

#plotlet(v, "cropID")


setMethod("plotlet", signature(x="SpatRaster"),
	function(x, y="", col="Spectral", opacity=0.8, tiles="", m=NULL, ...)  {
		if (is.null(m)) {
			m <- leaflet::leaflet()
		} else {
			tiles <- ""
		}
		if (!all(tiles == "")) {
			if ("Streets" %in% tiles) {
				m <- leaflet::addTiles(m, group="Streets")			
			} 
			tiles2 <- tiles[tiles != "Streets"]
			if (length(tiles) > 0) {
				for (i in 1:length(tiles2)) {
					m <- leaflet::addProviderTiles(m, tiles2[i], group=tiles2[i])
				}
			}
		}
		leaflet::addRasterImage(m, x, colors = col, opacity=opacity)
	}
)

