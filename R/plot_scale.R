
..assume_lonlat <- function(pr) {
	(pr$usr[1] > -300) && (pr$usr[2] < 300) && (pr$yaxp[1] > -200) && (pr$yaxp[2] < 200)
}

.get_dd <- function(pr, lonlat, d=NULL) {
	if (lonlat) {
		lat <- mean(pr$usr[3:4])
		if (is.null(d)) {
			dx <- (pr$usr[2] - pr$usr[1]) / 6
			d <- as.vector(distance(cbind(0, lat), cbind(dx, lat), TRUE))
			d <- max(1, 5 * round(d/5000))
		}
		p <- cbind(0, lat)
		dd <- .destPoint(p, d * 1000)
		dd <- c(dd[1,1], d)
	} else {
		if (is.null(d)) {
			d <- (pr$usr[2] - pr$usr[1]) / 6
			digits <- floor(log10(d)) + 1
			d <- round(d, -(digits-1))
		}
		dd <- c(d, d)
	}
	dd
}

.get_xy <- function(xy, dx=0, dy=0, pr, defpos="bottomleft", caller="") {
    if (is.null(xy)) {
		xy <- defpos
	}
	if (!is.character(xy)) {
		return( cbind(xy[1], xy[2]) )
	}
	xy <- tolower(xy)
	parrange <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
	pad=c(5,5) / 100
	if (xy == "bottom") {
		xy <- c(pr$usr[1]+0.5*parrange[1]-0.5*dx, pr$usr[3]+(pad[2]*parrange[2])) + c(0,dy)
	} else if (xy == "bottomleft") {
		xy <- c(pr$usr[1]+(pad[1]*parrange[1]), pr$usr[3]+(pad[2]*parrange[2])) + c(0,dy)
	} else if (xy == "bottomright") {
		xy <- c(pr$usr[2]-(pad[1]*parrange[1]), pr$usr[3]+(pad[2]*parrange[2])) - c(dx,-dy)
	} else if (xy == "topright") {
		xy <- c(pr$usr[2]-(pad[1]*parrange[1]), pr$usr[4]-(pad[2]*parrange[2])) - c(dx,dy)
	} else if (xy == "top") {
		xy <- c(pr$usr[1]+0.5*parrange[1]-0.5*dx, pr$usr[4]-(pad[2]*parrange[2])) - c(0,dy)
	} else if (xy == "topleft") {
		xy <- c(pr$usr[1]+(pad[1]*parrange[1]), pr$usr[4]-(pad[2]*parrange[2])) - c(0,dy)
	} else if (xy == "left") {
		xy <- c(pr$usr[1]+(pad[1]*parrange[1]), pr$usr[3]+0.5*parrange[2]-0.5*dy)
	} else if (xy == "right") {
		xy <- c(pr$usr[2]-(pad[1]*parrange[1])-dx, pr$usr[3]+0.5*parrange[2]-0.5*dy)
	} else {
		error(caller, 'xy must be a coordinate pair (two numbers) or one of "bottomleft", "bottom", "bottomright", topleft", "top", "topright"')
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




add_N <- function(x, y, asp, label, type=0, user="", angle=0, cex=1, srt=0, xpd=TRUE, ...) {

	type <- type[1]
	if (type == 0) { symbol = user[1]
	} else if (type == 2) { symbol = "\u27A2"
	} else if (type == 3) { symbol = "\u2799"
	} else if (type == 4) { symbol = "\u27B2"
	} else if (type == 5) { symbol = "\u27BE"
	} else if (type == 6) { symbol = "\u27B8"
	} else if (type == 7) { symbol = "\u27BB"
	} else if (type == 8) { symbol = "\u27B5"
	} else if (type == 9) { symbol = "\u279F"
	} else if (type == 10) { symbol = "\u261B"
	} else if (type == 11) { symbol = "\u2708"
	} else { symbol = "\u2629"}
	if (type == 11) {
		rangle <- 45 - angle
		mcex <- 1.5
	} else {
		rangle <- 90 - angle
		mcex <- 3
	}
	text(x, y, symbol, cex=cex*mcex, srt=rangle, xpd=xpd, ...)
	xs <- graphics::strwidth(symbol,cex=cex*3)
	ys <- graphics::strheight(symbol,cex=cex*3)
	b <- pi * angle / 180
	rxs <- (abs(xs * cos(b)) + abs(ys * sin(b)))# / asp
	rys <- (abs(xs * sin(b)) + abs(ys * cos(b)))# * asp
#	xoff <- (rxs - xs) / 2
#	yoff <- rys + 0.05 * graphics::strheight(label,cex=cex)
	xoff = 0.1 * rxs
	yoff = 0.8 * rys *  max(0.5, abs(cos(angle)))

    if (type == 4) {
        .halo(x+xoff, y-0.2*yoff, label, cex = cex, srt = srt, xpd = xpd, ...)
	} else if (type == 10) {
        .halo(x+xoff, y-yoff, label, cex = cex, srt = srt, xpd = xpd, ...)
    } else {
        text(x+xoff, y+yoff, label, cex = cex, srt = srt, xpd = xpd, ...)
    }
}


north <- function(xy=NULL, type=1, label="N", angle=0, d, head=0.1, xpd=TRUE, ...) {
	pr <- graphics::par()
	pr$usr <- unlist(get.clip()[1:4])
	pa <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
	asp <- pa[2]/pa[1]
	if (missing(d))	{
		d <- 0.07 * pa[2]
	}
	xy <- .get_xy(xy, 0, d, pr, "topright", caller="arrow")

	if (inherits(type, "character")) {
		usertype <- type
		type = 0
	} else {
		type <- round(type)
		usertype <- ""
	}
	if (type == 1) {
		if (angle != 0) {
			b <- angle * pi / 180;
			p2 <- xy + c(d * sin(b), d * cos(b))
			b <- b + pi
			p1 <- xy + c(d * sin(b), d * cos(b))
			if ((p2[1] - p1[1]) > (d/asp)) {
				m <- xy[1] #p1[1] + (p2[1] - p1[1]) / 2
				slope = (p2[2] - p1[2])/(p2[1] - p1[1])
				newx <- m - 0.5 * d / asp
				p1[2] <- p1[2] + (newx-p1[1]) * slope
				p1[1] <- newx
				newx <- m + 0.5 * d / asp
				p2[2] <- p2[2] - (p2[1]-newx) * slope
				p2[1] <- newx
			}
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
	} else {
		add_N(xy[1], xy[2], asp=asp, label=label, angle=angle, type=type, user=usertype, xpd=xpd, ...)
	}
}


sbar <- function(d, xy=NULL, type="line", divs=2, below="", lonlat=NULL, labels, adj=c(0.5, -1), lwd=2, xpd=TRUE, ticks=FALSE, scaleby=1, halo=TRUE, ...){

	stopifnot(type %in% c("line", "bar"))
	pr <- graphics::par()
	clp <- get.clip()
	pr$usr <- unlist(clp[,1:4])
	if (is.null(lonlat)) {
		lonlat <- isTRUE(clp[[5]])
	}

	if (missing(d)) {
		labels <- NULL
		d <- NULL
	}
	dd <- .get_dd(pr, lonlat, d)
	d <- dd[2]
	dd <- dd[1]

	xy <- .get_xy(xy, dd, 0, pr, "bottomleft", caller="sbar")

	if (type == "line") {
		if (halo) {
			lines(matrix(c(xy[1], xy[2], xy[1]+dd, xy[2]), byrow=T, nrow=2), lwd=lwd+1, xpd=xpd, col="white")
		}
		lines(matrix(c(xy[1], xy[2], xy[1]+dd, xy[2]), byrow=T, nrow=2), lwd=lwd, xpd=xpd, ...)
		
		if (missing(labels) || is.null(labels)) {
			ds <- d / scaleby
			if (divs > 2) {
				labels <- c(0, round(ds/2, 1), ds)
			} else {
				labels <- paste(ds)
			}
		}
		if (missing(adj)) {
			adj <- c(0.5, -0.2-lwd/20 )
		}
		tadd <- 0
		if (!isFALSE(ticks)) {
			if (isTRUE(ticks)) {
				tadd <- dd / (15 * diff(pr$usr[1:2]) / diff(pr$usr[3:4]))
			} else {
				tadd <- ticks
			}
			if (length(labels) == 1) {
				xtick <- c(xy[1], xy[1]+dd)
			} else {
				xtick <- c(xy[1], xy[1]+dd/2, xy[1]+dd)			
			}
			for (i in 1:length(xtick)) {
				lines(rbind(c(xtick[i], xy[2]), c(xtick[i], xy[2]+tadd)), lwd=ceiling(lwd/2), ...)
			}
		}
		tadd <- max(0, tadd)
		if (length(labels) == 1) labels =c("", labels, "")
		if (halo) {
			.halo(xy[1], xy[2]+tadd,labels=labels[1], xpd=xpd, adj=adj, ...)
			.halo(xy[1]+0.5*dd, xy[2]+tadd,labels=labels[2], xpd=xpd, adj=adj,...)
			.halo(xy[1]+dd, xy[2]+tadd,labels=labels[3], xpd=xpd, adj=adj,...)
		} else {
			text(xy[1], xy[2]+tadd,labels=labels[1], xpd=xpd, adj=adj, ...)
			text(xy[1]+0.5*dd, xy[2]+tadd,labels=labels[2], xpd=xpd, adj=adj,...)
			text(xy[1]+dd, xy[2]+tadd,labels=labels[3], xpd=xpd, adj=adj,...)
		}
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
			if (missing(labels) || is.null(labels)) {
				labels <- c("0", "", d/scaleby)
			}

			text(xy[1], xy[2],labels=labels[1], xpd=xpd, adj=adj,...)
			text(xy[1]+0.5*dd, xy[2],labels=labels[2], xpd=xpd, adj=adj,...)
			text(xy[1]+dd, xy[2],labels=labels[3], xpd=xpd, adj=adj,...)
		} else {
			q1 <- xy[1] + dd / 4
			half <- xy[1] + dd / 2
			q3 <- xy[1] + 3 * dd / 4
			end <- xy[1] + dd
			graphics::polygon(c(xy[1], xy[1], q1, q1), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="white", xpd=xpd)
			graphics::polygon(c(q1, q1, half, half), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="black", xpd=xpd)
			graphics::polygon(c(half, half, q3, q3 ), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="white", xpd=xpd)
			graphics::polygon(c(q3, q3, end, end), c(xy[2], xy[2]+lwd, xy[2]+lwd, xy[2]), col="black", xpd=xpd)
			if (missing(labels) || is.null(labels)) {
				ds <- d / scaleby
				labels <- c("0", round(0.5*ds), ds)
			}
			text(xy[1], xy[2], labels=labels[1], xpd=xpd, adj=adj, ...)
			text(half, xy[2], labels=labels[2], xpd=xpd, adj=adj,...)
			text(end, xy[2],labels=labels[3], xpd=xpd, adj=adj,...)
		}
	}
	if (below != "") {
		adj[2] <- -adj[2]
		text(xy[1]+(0.5*dd), xy[2], xpd=xpd, labels=below, adj=adj,...)
	}
}


