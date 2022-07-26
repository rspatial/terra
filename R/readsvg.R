
transform <- function(xy, m) {
	newX = m[1] * xy[,1] + m[3] * xy[,2] + m[5]
	newY = m[2] * xy[,1] + m[4] * xy[,2] + m[6]
	cbind(newX, newY)
}

oneline <- function(x, id) {
	x <- trimws(gsub('^m', "", x))
	ss <- trimws(unlist(strsplit(x, "m ")))
	out <- list()
	for (j in 1:length(ss)) {
		v <- unlist(utils::read.table(text=ss[j], sep=" "))
		vv <- as.numeric(unlist(strsplit(v, ",")))
		vv <- matrix(vv, ncol=2, byrow=TRUE)
		if (j > 1) {
			vv[1,] <- vv[1,] + a[1,]
		}
		a <- apply(vv, 2, cumsum)
		out[[j]] <- a
	}
	out <- lapply(1:length(out), function(p) cbind(id=id, part=p, out[[p]], hole=0))
	out <- do.call(rbind, out)
	out[,4] <- -out[,4]
	#out[,3:4] <- transform(out[,3:4], m)
	out
}


oneline2 <- function(x, id) {
	x <- trimws(gsub('^M', "", x))
	ss <- trimws(unlist(strsplit(x, "ZM")))
	out <- list()
	for (j in 1:length(ss)) {
		v <- unlist(strsplit(ss[j], "L"))
		v <- unlist(strsplit(v, "C"))
		v <- unlist(strsplit(v, " "))
		v <- gsub("Z", "", v)
		vv <- as.numeric(unlist(strsplit(v, ",")))
		vv <- matrix(vv, ncol=2, byrow=TRUE)
		out[[j]] <- vv
	}
	out <- lapply(1:length(out), function(p) cbind(id=id, part=p, out[[p]], hole=0))
	do.call(rbind, out)
	#out[,3:4] <- transform(out[,3:4], m)
}

readSVG <- function(f) {
	doc <- XML::htmlParse(f)
	p <- XML::xpathSApply(doc, "//path", XML::xmlGetAttr, "d")
	s <- list()
	for (i in 1:length(p)) {
		s[[i]] <- oneline2(p[i], i)
	}
	ss <- do.call(rbind, s)
	v <- vect(ss, type="polygons")

	a <- XML::xpathSApply(doc, "//path", XML::xmlAttrs)
	a <- unique(unlist(sapply(a, names)))
	a <- a[-grep(":", a)]
	a <- a[a != "d"]
	if (length(a) > 0) {
		att <- list()
		for (i in 1:length(a)) {
			z <- XML::xpathSApply(doc, "//path", XML::xmlGetAttr, a[i])
			att[[i]] <- sapply(z, function(i) if (is.null(i)) NA else i, USE.NAMES = FALSE)
		}
		names(att) <- a
		values(v) <- data.frame(att)
	}
	v
}




