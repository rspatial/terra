
## update — non-RGB raster: names, crs, extent persist
fname1 <- paste0(tempfile(), ".tif")
s <- rast(system.file("ex/logo.tif", package = "terra"))
RGB(s) <- NULL
x <- writeRaster(s, fname1)

ext(x) <- ext(x) + 1
crs(x) <- "+proj=utm +zone=1"
names(x) <- LETTERS[1:3]
update(x, names = TRUE, crs = TRUE, extent = TRUE)

y <- rast(fname1)
expect_equal(names(y), c("A", "B", "C"))
expect_equivalent(ext(y)[1], -1)
expect_equivalent(crs(y, describe = TRUE)$code, "32601")
file.remove(fname1)


## update — RGB raster: names persist
fname2 <- paste0(tempfile(), ".tif")
s2 <- rast(system.file("ex/logo.tif", package = "terra"))
x2 <- writeRaster(s2, fname2)

names(x2) <- c("X", "Y", "Z")
update(x2, names = TRUE)

y2 <- rast(fname2)
expect_equal(names(y2), c("X", "Y", "Z"))
aux <- paste0(fname2, ".aux.xml")
file.remove(fname2)
if (file.exists(aux)) file.remove(aux)
