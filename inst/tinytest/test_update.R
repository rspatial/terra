
## update meta — non-RGB raster: names, crs, extent persist
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


## update meta — RGB raster: names persist
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


## update values — single layer
f3 <- paste0(tempfile(), ".tif")
x3 <- writeRaster(rast(nrows = 10, ncols = 10, vals = 1:100), f3)
update(x3, cells = c(1, 50, 100), values = c(999, 888, 777))
y3 <- rast(f3)
expect_equal(as.numeric(values(y3)[1, ]),  999)
expect_equal(as.numeric(values(y3)[50, ]), 888)
expect_equal(as.numeric(values(y3)[100,]), 777)
expect_equal(as.numeric(values(y3)[2, ]),  2)
file.remove(f3)


## update values — all layers broadcast
f4 <- paste0(tempfile(), ".tif")
x4 <- writeRaster(rast(nrows = 5, ncols = 5, nlyr = 3, vals = as.double(1:75)),
                   f4, datatype = "FLT4S")
update(x4, cells = c(1, 13), values = c(-5, -10))
y4 <- rast(f4)
v4 <- values(y4)
expect_true(all(v4[1, ]  == -5))
expect_true(all(v4[13, ] == -10))
expect_equal(as.numeric(v4[2, 1]), 2)
file.remove(f4)


## update values — specific layer only
f5 <- paste0(tempfile(), ".tif")
x5 <- writeRaster(rast(nrows = 5, ncols = 5, nlyr = 3, vals = as.double(1:75)),
                   f5, datatype = "FLT4S")
update(x5, cells = c(1, 25), values = c(0, 0), layer = 2)
y5 <- rast(f5)
v5 <- values(y5)
expect_equal(as.numeric(v5[1, ]),  c(1, 0, 51))
expect_equal(as.numeric(v5[25, ]), c(25, 0, 75))
file.remove(f5)


## update values — per-layer values (cs * nlyrs)
f6 <- paste0(tempfile(), ".tif")
x6 <- writeRaster(rast(nrows = 5, ncols = 5, nlyr = 2, vals = as.double(1:50)),
                   f6, datatype = "FLT4S")
update(x6, cells = c(1, 2), values = c(100, 200, 300, 400))
y6 <- rast(f6)
v6 <- values(y6)
expect_equal(as.numeric(v6[1, ]), c(100, 300))
expect_equal(as.numeric(v6[2, ]), c(200, 400))
file.remove(f6)


## update values — single value broadcast to multiple cells
f7 <- paste0(tempfile(), ".tif")
x7 <- writeRaster(rast(nrows = 5, ncols = 5, vals = as.double(1:25)),
                   f7, datatype = "FLT4S")
update(x7, cells = c(1, 2, 3, 4, 5), values = 0)
y7 <- rast(f7)
expect_true(all(values(y7)[1:5, ] == 0))
expect_equal(as.numeric(values(y7)[6, ]), 6)
file.remove(f7)
