
# multiple sources
s <- rast(system.file("ex/logo.tif", package="terra"))
plotRGB(c(s[[1]], s[[2]], s[[3]]))

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

