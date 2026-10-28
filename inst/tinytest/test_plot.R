
# multiple sources
s <- rast(system.file("ex/logo.tif", package="terra"))
plotRGB(c(s[[1]], s[[2]], s[[3]]))

# blend: grayscale hill-like layer, RGB, user function, and resample
r <- rast(nrows=6, ncols=6, xmin=0, xmax=6, ymin=0, ymax=6, crs="local", vals=1:36)
h <- rast(r, vals=rep(seq(0.3, 1, length.out=6), each=6))
png(tempfile(fileext=".png"))
plot(r, blend=h, legend=FALSE, axes=FALSE, mar=NA)
plot(r, blend=h, blendfun="screen", legend=FALSE, axes=FALSE, mar=NA)
plot(r, blend=h, blendfun="overlay", legend=FALSE, axes=FALSE, mar=NA)
plot(r, blend=h, blendfun=function(a, b) a * b, legend=FALSE, axes=FALSE, mar=NA)
hrgb <- c(h, h, h)
plot(r, blend=hrgb, legend=FALSE, axes=FALSE, mar=NA)
h2 <- rast(nrows=3, ncols=3, ext=ext(r), crs="local", vals=0.5)
plot(r, blend=h2, legend=FALSE, axes=FALSE, mar=NA)
ha <- rast(r, vals=0.5)
plot(r, blend=h, alpha=ha, legend=FALSE, axes=FALSE, mar=NA)
h4 <- c(h, h, h, ha)
plot(r, blend=h4, legend=FALSE, axes=FALSE, mar=NA)
dev.off()
expect_error(plot(r, blend=1, legend=FALSE, axes=FALSE, mar=NA))
expect_error(plot(r, blend=h, blendfun="nope", legend=FALSE, axes=FALSE, mar=NA))

# colNA for continuous vector plot (#2198)
v <- vect(system.file("ex/lux.shp", package="terra"))
v[5, "AREA"] <- NA
png(tempfile(fileext=".png"))
plot(v, "AREA", colNA="grey", type="continuous")
dev.off()
out <- list(v=c(1, 2, NA, 4), range=c(1, 4), cols=grDevices::terrain.colors(5),
	leg=list(border="black"), colNA="grey")
out <- terra:::.vect.legend.continuous(out)
expect_equal(out$main_cols[3], "grey")
expect_true(all(!is.na(out$main_cols)))

# continuous legend: same size, different vertical alignment
mkleg <- function(yalign) {
	list(clip=TRUE, lim=c(0, 1, 0, 100),
		leg=list(reverse=FALSE, horiz=FALSE, size=0.75,
			ext=c(1, 1.1, 0, 100), yalign=yalign))
}
et <- terra:::.get.leg.coords(mkleg("top"))$leg$ext
ec <- terra:::.get.leg.coords(mkleg("center"))$leg$ext
eb <- terra:::.get.leg.coords(mkleg("bottom"))$leg$ext
expect_equal(et$dy, 75)
expect_equal(ec$dy, 75)
expect_equal(eb$dy, 75)
expect_equal(et$ymax, 100)
expect_equal(eb$ymin, 0)
expect_equal((ec$ymin + ec$ymax) / 2, 50)

