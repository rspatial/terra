
# these methods require the dev version of leaflet


popUp <- function(x) {
	nms <- names(x)
	s <- sapply(1:length(nms), function(i) paste0(nms[i], ": ", x[[i, drop=TRUE]]))
	apply(s, 1, function(i) paste(i, collapse="<br>"))
}


setMethod("plet", signature(x="SpatVector"),
	function(x, y="", col, alpha=1, tiles=c("Streets", "Esri.WorldImagery", "OpenTopoMap"), legend="bottomright", popup=TRUE, split=FALSE, collapse=FALSE, cex=1, map=NULL)  {
	
		if (missing(col)) col = grDevices::rainbow		
		alpha <- max(0, min(1, alpha))
		pop <- NULL
		
		#stopifnot(packageVersion("leaflet") > "2.1.1")
		if (is.null(map)) {
			map <- leaflet::leaflet()
			tiles <- unique(tiles)
			tiles <- tiles[tiles!=""]
		} else {
			tiles <- NULL
		}
		g <- geomtype(x)
		if ((!is.null(tiles)) && (length(tiles) > 0)) {
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
			if (popup) {
				pop <- popUp(x)
			}
			if (g == "polygons") {
				map <- leaflet::addPolygons(map, data=x, label=1:nrow(x),  
							col=cols, fillOpacity=alpha, popup=pop)
			} else if (g == "lines") {
				map <- leaflet::addPolylines(map, data=x, label=1:nrow(x),  
							col=cols, fillOpacity=alpha, popup=pop)
			} else {
				map <- leaflet::addMarkers(map, data=x, label=1:nrow(x),  
							col=cols, fillOpacity=alpha, popup=pop)			
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
					if (popup) {
						pop <- popUp(s)
					}
					if (g == "polygons") {
						map <- leaflet::addPolygons(map, data=s, label=u[i], group=u[i], 
							col=cols[i], fillOpacity=alpha, popup=pop)
					} else if (g == "lines") {
						map <- leaflet::addPolylines(map, data=s, label=u[i], group=u[i], 
							col=cols[i], fillOpacity=alpha, popup=pop)
					} else {
						map <- leaflet::addCircleMarkers(map, data=s, label=u[i], group=u[i], 
							col=cols[i], fillOpacity=alpha, popup=pop))					
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
				if (popup) {
					pop <- popUp(x)
				}
				if (g == "polygons") {
					map <- leaflet::addPolygons(map, data=x, label=values,  
						col=vcols, fillOpacity=alpha, popup=pop)
				} else if (g == "lines") {
					map <- leaflet::addPolylines(map, data=x, label=values,  
						col=vcols, popup=pop, fillOpacity=alpha)
				} else {
					map <- leaflet::addCircleMarkers(map, data=x, label=values,  
						col=vcols, radius=cex, popup=pop, fillOpacity=alpha)
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
	function(x, y=1, col, alpha=0.8, tiles=NULL, maxcell=500000, legend="bottomright", shared=FALSE, collapse=TRUE, map=NULL)  {
		#stopifnot(packageVersion("leaflet") > "2.1.1")

		alpha <- max(0, min(1, alpha))
		if (is.null(map)) {
			map <- leaflet::leaflet()
			tiles <- unique(tiles)
			tiles <- tiles[tiles!=""]
		} else {
			tiles <- NULL
		}
		if (missing(col)) {
			col <- rev(grDevices::terrain.colors(255))
		}
		if (length(tiles) > 0) {
			tiles <- unique(tiles)
			if (length(tiles) > 1) {
				tiles <- tiles[1]
				warn("plet", "only a single tileset can be used with raster data")
			}
			if (tiles == "Streets") {
				map <- leaflet::addTiles(map)
			} else {
				map <- leaflet::addProviderTiles(map, tiles)
			}
		}
		
		x <- spatSample(x[[y]], maxcell, "regular", as.raster=TRUE)
		if (nlyr(x) == 1) {
			map <- leaflet::addRasterImage(map, x, colors = col, opacity=alpha)
			if (!is.null(legend)) {
				r <- minmax(x)
				v <- seq(r[1], r[2], 5)
				pal <- leaflet::colorNumeric(col, v, reverse = TRUE)
				map <- leaflet::addLegend(map, legend, pal=pal, values=v, opacity=1,
					  labFormat = leaflet::labelFormat(transform = function(x) sort(x, decreasing = TRUE)))	
			}
		} else {
			nms <- make.unique(names(x))
			many_legends <- one_legend <- FALSE
			if (!is.null(legend)) {
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
					map <- leaflet::addRasterImage(map, x[[i]], colors=pal, opacity=alpha, group=nms[i])
				} else {
					map <- leaflet::addRasterImage(map, x[[i]], colors=col, opacity=alpha, group=nms[i])
					if (many_legends) {
						v <- seq(r[1,i], r[2,i], length.out=5)
						pal <- leaflet::colorNumeric(col, v, reverse=TRUE)
						map <- leaflet::addLegend(map, position=legend, pal=pal, values=v, 
							  title=nms[i], opacity=1, group=nms[i],
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


