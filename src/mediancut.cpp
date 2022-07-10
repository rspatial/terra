// Copyright (c) 2018-2022  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

/*
#include "spatRaster.h"


make_cut <- function(x) {
	j <- length(x)
	out <- vector("list", 2*j)
	for (i in 1:j) {
		rgb <- x[[i]]
		if (NROW(rgb) <= 1) {
			out[[i]] <- rgb
			j <- j - 1
			next
		}
		rng <- apply(rgb[,-1], 2, function(i) diff(range(i)))
		if (max(rng) == 0) {
			out[[i]] <- rgb
			j <- j - 1
			next
		}
		p <- which.max(rng) + 1
		m <- median(rgb[,p])
		out[[i]] <- rgb[rgb[,p] >= m, ,drop=FALSE]
		out[[i+j]] <- rgb[rgb[,p] < m, ,drop=FALSE]
	}
	i <- sapply(out, is.null)
	out <- out[!i]
	i <- sapply(out, nrow) > 0
	out[i]
}

median_cut <- function(v) {
	v <- list(v)
	n <- 0
	while ((length(v) < 129) & (length(v) > n)) {
		n <- length(v)
		v <- make_cut(v)
	}
	s <- sapply(v, function(i) max(apply(i[,-1,drop=FALSE], 2, function(j) diff(range(j)))))
	n <- 256 - length(v)
	ss <- rev(sort(s))
	ss <- max(2, min(ss[1:n]))
	i <- which(s > ss)
	if (length(i) > 0) {
		vv <- make_cut(v[i])
		v <- c(v[-i], vv)
	}
	v <- lapply(1:length(v), function(i) cbind(group=i, v[[i]]))
	do.call(rbind, v)
}


SpatRatser SpatRaster::RGB2col(std::string stretch, SpatOptions &opt) {
		std::vector<int> idx = x.RGB();
		SpatRaster out = geometry(1);
		if (idx.size() != 3) {
			out.setError("x does not have a valid RGB attribute");
			return out;
		}
		if (vmax(idx) >= x.nlyr()) {
			out.setError("invalid RGB indices")
		}
		*this = subset(idx);

		if (stretch != "") {
			if (stretch == "lin") {

			} else if (stretch == "hist") {

			} else {
				out.addWarning("invalid stretch option");
			}
		}

		v <- cbind(id=1:ncell(x), values(x))
		v <- median_cut(stats::na.omit(v))

		a <- aggregate(v[,3:5], list(v[,1]), median)
		a$cols <- grDevices::rgb(a[,2], a[,3], a[,4], maxColorValue=255)
		m <- merge(v[,1:2], a[, c(1,5)], by=1)
		r <- rast(x, 1)
		r[m$id] <- m$group - 1
		coltab(r) <- a$cols
		if (filename != "") {
			r <- writeRaster(r, filename, overwrite, ...)
		}
		r
	}
)


*/
