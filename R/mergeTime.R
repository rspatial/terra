
if (!isGeneric("fillTime")) {setGeneric("fillTime", function(x, ...) standardGeneric("fillTime"))}

setMethod("fillTime", signature(x="SpatRaster"),
	function(x, filename="", ...)  {
		tm <- time(x)
		if (any(is.na(tm))) {
			error("fillTime", "NA in time values")
		}
		if (any(table(tm)>1)) {
			error("fillTime", "duplicate time values")
		}
		if (is.unsorted(tm)) {
			warn("mergeTimelines", "sorting layers")
			ord <- order(tm)
			x <- x[[ ord ]]
			tm <- tm[ord]
		}
		d <- data.frame(time=seq(min(tm), max(tm), min(diff(tm))))
		d <- merge(d, data.frame(time=tm, tm=tm), by=1, all.x=TRUE)
		b <- (!is.na(d[,2])) + 0
		b <- cumsum(b) * b
		if (any(b==0)) {
			mx <- max(b)
			b[b==0] <- mx+1
			r <- init(rast(x, nlyr=1), NA)
			x <- c(x, r)
			x <- x[[b]]
		}
		time(x) <- d[,1]
		if (filename != "") {
			writeRaster(x, ...)
		} else {
			x
		}
	}
)


if (!isGeneric("mergeTime")) {setGeneric("mergeTime", function(x, ...) standardGeneric("mergeTime"))}

setMethod("mergeTime", signature(x="SpatRasterDataset"),
	function(x, fun="mean", filename="", ...)  {
		tim <- lapply(1:length(x), function(i) time(x[i]))
		if (any(sapply(tim, function(i) any(is.na(i))))) {
			error("mergeTime", "NA in time values")
		}
		if (any(sapply(tim, function(i) any(table(i)>1)))) {
			error("mergeTime", "duplicate time values")
		}
		us <- sapply(tim, is.unsorted)
		if (any(us)) {
			warn("mergeTime", paste("sorting layers of SpatRaster:", paste(us, collapse=", ")))
			us <- which(us)
			for (i in us) {
				ord <- order(tim[[i]])
				x[i] <- x[i][[ ord ]]
				tim[[i]] <- tim[[i]][ord]
			}
		}
		z <- data.frame(time=sort(unique(do.call(c, tim))))
		for (i in 1:length(tim)) {
			d <- data.frame(tim[[i]],tim[[i]])
			names(d) <- c("time", paste0("x", i))
			z <- merge(z, d, by=1, all.x=TRUE)
		}
		b <- (!is.na(z)) + 0
		y <- apply(b, 1, function(i) paste(i, collapse=""))
		r <- rep(1, length(y))
		for (i in 2:length(y)) {
			if (y[i] == y[i-1]) {
				r[i] <- r[i-1]
			} else {
				r[i] <- r[i-1] + 1
			}
		}
		u <- unique(r)
		out <- list()
		d <- apply(b, 2, cumsum) * b
		for (i in u) {
			zz <- z[r==i, ,drop=FALSE]
			dd <- d[r==i, -1, drop=FALSE]
			tim <- zz[,1]
			s <- which(colSums(is.na(zz[,-1])) == 0)
			if (length(s) == 1) {
				out[[i]] <- x[s][[ dd[,s] ]]
			} else {
				ss <- x[s]
				for (j in 1:length(s)) {
					ss[j] = ss[j][[ dd[, s[j]] ]]
				}
				out[[i]] <- app(ss, fun)
				time(out[[i]]) <- tim
			}
		}
		out <- rast(out)
		if (filename != "") {
			out <- writeRaster(out, ...)
		}
		out
	}
)


