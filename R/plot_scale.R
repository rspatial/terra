
.assume_lonlat <- function(pr) {
	(pr$usr[1] > -181) && (pr$usr[2] < 181) && (pr$yaxp[1] > -200) && (pr$yaxp[2] < 200)
}

.get_dd <- function(pr, lonlat, d=NULL) {
	if (lonlat) {
		lat <- mean(pr$yaxp[1:2])
		if (is.null(d)) {
			dx <- (pr$usr[2] - pr$usr[1]) / 6
			d <- as.vector(distance(cbind(0, lat), cbind(dx, lat), TRUE))
			d <- signif(d / 1000, 2) 
		}
		p <- cbind(0, lat)
		dd <- .destPoint(p, d * 1000)
		dd <- dd[1,1]
	} else {
		if (is.null(d)) {
			d <- (pr$usr[2] - pr$usr[1]) / 6
			digits <- floor(log10(d)) + 1
			d <- round(d, -(digits-1))	
		}
		dd <- d
	}
	dd
}

.get_xy <- function(xy, dx=0, dy=0, pr, defpos="bottomleft", caller="") {
    if(is.null(xy)) {
		xy <- defpos
	}
	if (!is.character(xy)) {
		return( cbind(xy[1], xy[2]) )
	}	
	xy <- tolower(xy)
	parrange <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
	padding=c(5,5) / 100
	if (xy == "bottomleft") {
		xy <- c(pr$usr[1]+(padding[1]*parrange[1]), pr$usr[3]+(padding[2]*parrange[2])) + c(0,dy)
	} else if (xy == "bottomright") {
		xy <- c(pr$usr[2]-(padding[1]*parrange[1]), pr$usr[3]+(padding[2]*parrange[2])) - c(dx,dy)
	} else if (xy == "topright") {
		xy <- c(pr$usr[2]-(padding[1]*parrange[1]), pr$usr[4]-(padding[2]*parrange[2])) - c(dx,dy)
	} else if (xy == "topleft") {
		xy <- c(pr$usr[1]+(padding[1]*parrange[1]), pr$usr[4]-(padding[2]*parrange[2])) - c(0,dy)
	} else {
		error(caller, 'xy must be a coordinate pair (two numbers) or one of "bottomleft", "bottomright", topleft", "topright"')
	}
	xy
}

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


arrow <- function(d, xy=NULL, head=0.1, angle=0, label="N", lonlat=NULL, xpd=TRUE, ...) {
	pr <- graphics::par()
	if (is.null(lonlat)) {
		lonlat <- .assume_lonlat(pr)
	}

	if (missing(d))	{
		d <- .get_dd(pr, lonlat, NULL) / 3
	}

	xy <- .get_xy(xy, 0, d, pr, "topright", caller="arrow")	
	if (angle != 0) {
		b <- angle * pi / 180;
		p2 <- xy + c(d * sin(b), d * cos(b))
		b <- b + pi
		p1 <- xy + c(d * sin(b), d * cos(b))
	} else {
		p1 <- xy - c(0,d)
		p2 <- xy + c(0,d)
	}
	lwd <- list(...)$lwd + 2
	if (is.null(lwd)) lwd <- 3
	graphics::arrows(p1[1], p1[2], p2[1], p2[2], length=head, lwd=lwd, col="white", xpd=xpd)
	graphics::arrows(p1[1], p1[2], p2[1], p2[2], length=head, xpd=xpd, ...)
	if (label != "") {
		if (is.null(list(...)$hw)) {
			.halo(xy[1], xy[2], label, hw=.2, xpd=xpd, ... )
		} else {
			.halo(xy[1], xy[2], label, xpd=xpd, ... )		
		}
	}
}



sbar <- function(d, xy=NULL, type="line", divs=2, below="", lonlat=NULL, label, adj=c(0.5, -1), lwd=2, xpd=TRUE, ...){

	stopifnot(type %in% c("line", "bar"))
	pr <- graphics::par()
	if (is.null(lonlat)) {
		lonlat <- .assume_lonlat(pr)
	}

	if (missing(d)) {
		label <- NULL
		d <- NULL
	}
	dd <- .get_dd(pr, lonlat, d)
	if (is.null(d)) d <- dd
	xy <- .get_xy(xy, dd, 0, pr, "bottomleft", caller="sbar")

	if (type == "line") {
		lines(matrix(c(xy[1], xy[2], xy[1]+dd, xy[2]), byrow=T, nrow=2), lwd=lwd, xpd=xpd, ...)
		if (missing(label)) {
			label <- paste(d)
		}
		if (is.null(label)) {
			label <- paste(d)
		}
		if (missing(adj)) {
			adj <- c(0.5, -0.2-lwd/20 )
		}

		if (length(label) == 1) label =c("", label, "")
		text(xy[1], xy[2],labels=label[1], xpd=xpd, adj=adj,...)
		text(xy[1]+0.5*dd, xy[2],labels=label[2], xpd=xpd, adj=adj,...)
		text(xy[1]+dd, xy[2],labels=label[3], xpd=xpd, adj=adj,...)

		xy[2] <- xy[2] - dd/10

	} else if (type == "bar") {
		stopifnot(divs > 0)

		if (missing(adj)) {
			adj <- c(0.5, -1 )
		}
		lwd <- dd / 25

		if (divs==2) {
			half <- xy[1] + dd / 2
			graphics::polygon(c(xy[1], xy[1], half, half), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="white", xpd=xpd)
			graphics::polygon(c(half, half, xy[1]+dd, xy[1]+dd ), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="black", xpd=xpd)
			if (missing(label)) {
				label <- c("0", "", d)
			}
			if (is.null(label)) {
				label <- c("0", "", d)
			}

			text(xy[1], xy[2],labels=label[1], xpd=xpd, adj=adj,...)
			text(xy[1]+0.5*dd, xy[2],labels=label[2], xpd=xpd, adj=adj,...)
			text(xy[1]+dd, xy[2],labels=label[3], xpd=xpd, adj=adj,...)
		} else {
			q1 <- xy[1] + dd / 4
			half <- xy[1] + dd / 2
			q3 <- xy[1] + 3 * dd / 4
			end <- xy[1] + dd 
			graphics::polygon(c(xy[1], xy[1], q1, q1), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="white", xpd=xpd)
			graphics::polygon(c(q1, q1, half, half), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="black", xpd=xpd)
			graphics::polygon(c(half, half, q3, q3 ), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="white", xpd=xpd)
			graphics::polygon(c(q3, q3, end, end), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="black", xpd=xpd)
			if (missing(label)) {
				label <- c("0", round(0.5*d), d)
			}
			if (is.null(label)) {
				label <- c("0", round(0.5*d), d)
			}
			text(xy[1], xy[2], labels=label[1], xpd=xpd, adj=adj, ...)
			text(half, xy[2], labels=label[2], xpd=xpd, adj=adj,...)
			text(end, xy[2],labels=label[3], xpd=xpd, adj=adj,...)
		}
	}
	if (below != "") {
		adj[2] <- -adj[2]
		text(xy[1]+(0.5*dd), xy[2], xpd=xpd, labels=below, adj=adj,...)
	}
}
