# Author: Jacob van Etten
# email jacobvanetten@yahoo.com
# Date :  May 2010
# Version 1.1
# Licence GPL v3

# RH: updated for igraph (from igraph0)
# sept 23, 2012

# RH: adapted from raster to terra 
# Jan 11, 2022


.calcDist <- function(r, origin=NULL, omit=NULL, filename, ...) {
	x <- rast(r)
	lonlat <- isTRUE(is.lonlat(x))
	chunkSize <- ncell(x)
	
	ftC <- values(r, mat=FALSE)
	if (is.null(origin)) {
		oC <- which(!is.na(ftC)) 
	} else {
		oC <- which(ftC %in% origin) 	
	}
	if (length(oC) == 0) {
		error("gridDistance", "no cells to compute distance to")
	} else if (length(oC) == ncell(x)) {
		return( init(x, 0, filename=filename, ...))
	}
	ftC <- which(!(ftC %in% omit)) 	
	
	adj <- adjacent(x, ftC, directions=8, pairs=TRUE)
	if (!is.null(ftC)) {
		adj <- adj[adj[,2] %in% ftC, ,drop=FALSE]
	}
	startNode <- max(adj)+1 #extra node to serve as origin
	adjP <- rbind(adj, cbind(rep(startNode, times=length(oC)), oC))
	distGraph <- igraph::graph.edgelist(adjP, directed=TRUE)
	perCell <- rep(0, times=length(oC))
	if (lonlat) {
		distance <- distance(xyFromCell(x,adj[,1]), xyFromCell(x,adj[,2]), lonlat=TRUE, pairwise=TRUE) 
		igraph::E(distGraph)$weight <- c(distance, perCell)

	} else {
		sameRow <- which(rowFromCell(x, adj[,1]) == rowFromCell(x, adj[,2]))
		sameCol <- which(colFromCell(x, adj[,1]) == colFromCell(x, adj[,2]))
		igraph::E(distGraph)$weight <- sqrt(xres(x)^2 + yres(x)^2)
		igraph::E(distGraph)$weight[sameRow] <- xres(x)
		igraph::E(distGraph)$weight[sameCol] <- yres(x)
		igraph::E(distGraph)$weight[(length(adj[,1])+1):(length(adj[,1])+length(oC))] <- perCell
	}
		
	shortestPaths <- igraph::shortest.paths(distGraph, startNode)
		
	if (length(shortestPaths) < chunkSize) { 
		#add Inf values where shortest.paths() leaves off before completing all nodes
		shortestPaths <- c(shortestPaths, rep(Inf, times=chunkSize-length(shortestPaths))) 
	}	
	shortestPaths <- shortestPaths[-(length(shortestPaths))] #chop startNode off
	x <- setValues(x, shortestPaths)
	classify(x, cbind(Inf, NA), filename, ...)
}




.calcDistChunk <- function(x, chunkSize, ftC, oC, perCell=0, startCell=0, lonlat) {

	if (length(oC) > 0) {
		#adj <- adjacency(x, fromCells=ftC, toCells=ftC, directions=8)
		adj <- adjacent(x, ftC, directions=8, pairs=TRUE)
		if (!is.null(ftC)) {
			adj <- adj[adj[,2] %in% ftC, ,drop=FALSE]
		}
		startNode <- max(adj)+1 #extra node to serve as origin
		adjP <- rbind(adj, cbind(rep(startNode, times=length(oC)), oC))
		distGraph <- igraph::graph.edgelist(adjP, directed=TRUE)
		if (length(perCell) == 1) {
			if (perCell == 0) {
				perCell <- rep(0, times=length(oC))
			}
		}

		if (lonlat) {
			distance <- distance(xyFromCell(x,adj[,1]+startCell), xyFromCell(x,adj[,2]+startCell), lonlat=TRUE, pairwise=TRUE) 
			igraph::E(distGraph)$weight <- c(distance, perCell)

		} else {
			sameRow <- which(rowFromCell(x, adj[,1]) == rowFromCell(x, adj[,2]))
			sameCol <- which(colFromCell(x, adj[,1]) == colFromCell(x, adj[,2]))
			igraph::E(distGraph)$weight <- sqrt(xres(x)^2 + yres(x)^2)
			igraph::E(distGraph)$weight[sameRow] <- xres(x)
			igraph::E(distGraph)$weight[sameCol] <- yres(x)
			igraph::E(distGraph)$weight[(length(adj[,1])+1):(length(adj[,1])+length(oC))] <- perCell
		}
		
		shortestPaths <- igraph::shortest.paths(distGraph, startNode)
		shortestPaths <- shortestPaths[-(length(shortestPaths))] #chop startNode off
		
		if (length(shortestPaths) < chunkSize) { 
			#add Inf values where shortest.paths() leaves off before completing all nodes
			shortestPaths <- c(shortestPaths, rep(Inf, times=chunkSize-length(shortestPaths))) 
		}
		
	} else {
		shortestPaths <- rep(Inf, times=chunkSize)
	}
	
	return(shortestPaths)
}



setMethod("gridDistance", signature(x="SpatRaster"), 
function(x, origin=NULL, omit=NULL, chunk=FALSE, filename="", overwrite=FALSE, ...) {

	if (! hasValues(x) ) {
		error("gridDistance", "SpatRaster has no cells values")
	}
	if (file.exists(filename) && (!isTRUE(overwrite))) {
		error("gridDistance", "file exists. Use 'overwrite=TRUE' if you want to overwrite it")
	}
	if( !requireNamespace("igraph")) {
		stop("the igraph package needs to be installed to use this function")
	}

	if (nlyr(x) > 1) {
		r <- list()
		for (i in 1:nlyr(x)) {
			r[[i]] <- gridDistance(x[[i]], origin=origin, omit=omit, chunk=chunk, filename="", ...) 
		}
		r <- rast(r)
		if (filename != "") {
			r <- writeRaster(r, filename, overwrite, ...)
		}
		return(r)
	}

	if (!chunk) {
		x <- try(.calcDist(x, origin, omit, filename, ...))
		if (inherits(x, "try-error")) {
			error("gridDistance", "gridDistance failed. Perhaps the raster is too large (set 'chunk=TRUE')?")
		}
		messages(x)

	} else {
	#going up		
		lonlat <- isTRUE(is.lonlat(x))
		r1 <- rast(x)
		b <- writeStart(r1, "")
		readStart(x)
		for (i in b$n:1) {
			chunk <- readValues(x, row=b$row[i], nrows=b$nrows[i]) 
			startCell <- (b$row[i]-1) * ncol(x)
			chunkSize <- length(chunk)
			oC <- which(chunk %in% origin) 
			ftC <- which(!(chunk %in% omit))
			if (length(ftC) != 0) {

				if (i < b$n) {
					firstRowftC <- firstRowftC + chunkSize 
					chunkDist <- .calcDistChunk(x, 
								chunkSize=chunkSize + ncol(x), 
								ftC=c(ftC, firstRowftC), 
								oC=c(oC, firstRowftC), 
								perCell=c(rep(0,times=length(oC)),firstRowDist), 
								startCell=startCell,
								lonlat=lonlat)[1:chunkSize]
				} else {
					chunkDist <- .calcDistChunk(x, chunkSize=chunkSize, 
								ftC=ftC, oC=oC, perCell=0,
								startCell=startCell, lonlat=lonlat)
				}
			} else {
				if (i < b$n) {
					firstRowftC <- firstRowftC + chunkSize 
				}
				chunkDist <- rep(NA, b$nrows[i] * ncol(r1))
			}
			firstRow <- chunk[1:ncol(x)]
			firstRowDist <- chunkDist[1:ncol(x)]
			firstRowftC <- which(!(firstRow %in% omit))
			firstRowDist <- firstRowDist[firstRowftC]
			chunkDist[is.infinite(chunkDist)] <- NA

			writeValues(r1, chunkDist, b$row[i], b$nrows[i])
		}
		writeStop(r1)
		readStart(r1)
		#going down
		
		out <- rast(x)
		b <- writeStart(out, filename=filename, overwrite=TRUE, ...)			
		for (i in 1:b$n) {
			chunk <- readValues(x, row=b$row[i], nrows=b$nrows[i]) 
			chunkSize <- length(chunk)
			startCell <- (b$row[i]-1) * ncol(x)
			oC <- which(chunk %in% origin) 
			ftC <- which(!(chunk %in% omit))
			
			if (length(ftC) != 0) {
			
				if (i > 1) {
					chunkDist <- readValues(r1, row=b$row[i], nrows=b$nrows[i]) 
					chunkDist[is.na(chunkDist)] <- Inf 
				
					chunkDist <- pmin(chunkDist,
						.calcDistChunk(x, chunkSize=chunkSize+ncol(x), 
							ftC = c(lastRowftC, ftC+ncol(x)), 
							oC = c(lastRowftC, oC+ncol(x)), 
							perCell = c(lastRowDist, rep(0,times=length(oC))), 
							startCell = startCell - ncol(x),
							lonlat=lonlat)[-(1:ncol(r1))])
							
				} else {
					chunkDist <- readValues(r1, row=b$row[i], nrows=b$nrows[i])
					chunkDist[is.na(chunkDist)] <- Inf
			
					chunkDist <- pmin(chunkDist,
						.calcDistChunk(x, chunkSize=chunkSize, 
							ftC=ftC, oC=oC, perCell=0, 
							startCell=startCell, lonlat=lonlat))
				}
			} else {			
				chunkDist <- rep(NA, b$nrows[i] * ncol(out))						
			}

			lastRow <- chunk[(length(chunk)-ncol(x)+1):length(chunk)]
			lastRowDist <- chunkDist[(length(chunkDist)-ncol(x)+1):length(chunkDist)]
			lastRowftC <- which(!(lastRow %in% omit))
			lastRowDist <- lastRowDist[lastRowftC]
			chunkDist[is.infinite(chunkDist)] <- NA
			writeValues(out, chunkDist, b$row[i], b$nrows[i])
		}
		readStop(x)
		writeStop(out)
	}
}
)



