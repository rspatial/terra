

.destPoint <- function (p, d, b=90, r=6378137) {
    toRad <- pi/180
    lon1 <- p[, 1] * toRad
    lat1 <- p[, 2] * toRad
    b <- b * toRad
    lat2 <- asin(sin(lat1) * cos(d/r) + cos(lat1) * sin(d/r) * cos(b))
    lon2 <- lon1 + atan2(sin(b) * sin(d/r) * cos(lat1), cos(d/r) - sin(lat1) * sin(lat2))
    lon2 <- (lon2 + pi)%%(2 * pi) - pi
    cbind(lon2, lat2)/toRad
}


.arrow <- function(d, xy=click(), head=0.1, ...) {
	graphics::arrows(xy[1], xy[2], xy[1], xy[2]+d, length=head, ...)
	lines(rbind(xy, rbind(cbind(xy[1], xy[2]-d))), ...)
	text(xy[1,1], xy[1,2]-(0.25*d), 'N')
}


.scalebar <- function(d, xy=NULL, type='line', divs=2, below='', lonlat=NULL, label, adj=c(0.5, -0.5), lwd=2, ...){

	stopifnot(type %in% c('line', 'bar'))
	pr <- graphics::par()
	if (is.null(lonlat)) {
		if ( pr$usr[1] > -181 & pr$usr[2] < 181 &  pr$yaxp[1] > -200 &  pr$yaxp[2] < 200  ) {
			lonlat <- TRUE
		} else {
			lonlat <- FALSE
		}
	}

	if (lonlat) {
		lat <- mean(pr$yaxp[1:2])
		if (missing(d)) {
			dx <- (pr$usr[2] - pr$usr[1]) / 10
			d <- pointDistance(cbind(0, lat), cbind(dx, lat), TRUE)
			d <- signif(d / 1000, 2) 
			label <- NULL
		}
		p <- cbind(0, lat)
		dd <- .destPoint(p, d * 1000)
		dd <- dd[1,1]
	} else {
		if (missing(d)) {
			d <- round(10*(pr$usr[2] - pr$usr[1])/10) / 10
			label <- NULL
		}
		dd <- d
	}
	
    if(is.null(xy)) {
		padding=c(5,5) / 100
		#defaults to a lower left hand position
		parrange <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
		xy <- c(pr$usr[1]+(padding[1]*parrange[1]), pr$usr[3]+(padding[2]*parrange[2]))
	}

	if (type == 'line') {
		lines(matrix(c(xy[1], xy[2], xy[1]+dd, xy[2]), byrow=T, nrow=2), lwd=lwd, ...)
		if (missing(label)) {
			label <- paste(d)
		}
		if (is.null(label)) {
			label <- paste(d)
		}
		if (missing(adj)) {
			adj <- c(0.5, -0.2-lwd/20 )
		}
		text(xy[1]+(0.5*dd), xy[2],labels=label, adj=adj,...)
		
		
	} else if (type == 'bar') {
		stopifnot(divs > 0)
		
		if (missing(adj)) {
			adj <- c(0.5, -1 )
		}
		lwd <- dd / 25
		
		if (divs==2) {
			half <- xy[1] + dd / 2
			graphics::polygon(c(xy[1], xy[1], half, half), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col='white')
			graphics::polygon(c(half, half, xy[1]+dd, xy[1]+dd ), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col='black')
			if (missing(label)) {
				label <- c('0', '', d)
			}
			if (is.null(label)) {
				label <- c('0', '', d)
			}
			
			text(xy[1], xy[2],labels=label[1], adj=adj,...)
			text(xy[1]+0.5*dd, xy[2],labels=label[2], adj=adj,...)
			text(xy[1]+dd, xy[2],labels=label[3], adj=adj,...)
		} else {
			q1 <- xy[1] + dd / 4
			half <- xy[1] + dd / 2
			q3 <- xy[1] + 3 * dd / 4
			end <- xy[1] + dd 
			graphics::polygon(c(xy[1], xy[1], q1, q1), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col='white')
			graphics::polygon(c(q1, q1, half, half), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col='black')
			graphics::polygon(c(half, half, q3, q3 ), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col='white')
			graphics::polygon(c(q3, q3, end, end), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col='black')
			if (missing(label)) {
				label <- c('0', round(0.5*d), d)
			}
			if (is.null(label)) {
				label <- c('0', round(0.5*d), d)
			}
			text(xy[1], xy[2], labels=label[1], adj=adj,...)
			text(half, xy[2], labels=label[2], adj=adj,...)
			text(end, xy[2],labels=label[3], adj=adj,...)
		}
			
		if (below != "") {
			adj[2] <- -adj[2]
			text(xy[1]+(0.5*dd), xy[2], labels=below, adj=adj,...)
		}
	}
}
