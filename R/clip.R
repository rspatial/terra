
set.clip <- function(clp, geo) {
	# remove non-existing ones
	x <- grDevices::dev.list()
	x <- paste(names(x), x, sep="_")
	e <- .terra_environment$devs
	e <- e[e[,1] %in% x, ]

	graphics::clip(clp[1], clp[2], clp[3], clp[4])

	d <- grDevices::dev.cur()
	d <- data.frame(dev=paste(names(d), d[[1]], sep="_"), rbind(clp), geo=geo, row.names="")
	# remove one with the same name/number 
	e <- e[!(e[,1] %in% d[1]), ]
	e <- rbind(e, d)
	.terra_environment$devs <- e
}

get.clip <- function() {
	d <- grDevices::dev.cur()
	dev <- paste(names(d), d[[1]], sep="_")
	e <- .terra_environment$devs
	i <- match(dev, e[,1])[1]
	if (is.na(i)) {
		NULL
	} else {
		e[i[1],-1]
	}
}

reset.clip <- function() {
	g <- get.clip()
	if (!is.null(g)) {
		graphics::clip(g[[1]], g[[2]], g[[3]], g[[4]])
	}
}

