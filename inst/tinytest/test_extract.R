
set.seed(500)
r <- rast(ext=c(0, 3, 0, 3), ncol = 3, nrow = 3, nlyr=2, vals=runif(18))
pts <- cbind(c(0.5, 1.5), c(0.5, 1.5))
vct <- terra::vect(pts)
e <- terra::extract(r, vct, ID = FALSE)
expect_equal(e[,1], c(0.5121819, 0.81227813))
e <- terra::extract(r, vct, ID = FALSE, fun=sum)
expect_equal(e[,1], c(0.5121819, 0.81227813))


f <- system.file("ex/lux.shp", package="terra")
y <- vect(f)[1:2,]

elev <- rast(system.file("ex/elev.tif", package = "terra"))
e <- extract(elev, y, fun=mean, na.rm=TRUE, ID=TRUE)
expect_equal(e[,2], c(467.10517, 333.86294))

e <- extract(elev, y, fun=mean, exact=TRUE, na.rm=TRUE, ID=FALSE)
expect_equal(e[,1], c(467.379239, 334.685564))

e <- extract(elev, y, fun=mean, weights=TRUE, na.rm=TRUE, ID=TRUE)
expect_equal(e[,2], c(467.3933629, 334.65513085))

e <- extract(elev, 1, raw=TRUE)
expect_true(inherits(e, "array"))

x <- rast(y, res=.2)
values(x) <- 1:ncell(x)
expect_equal(cells(x, y), cbind(ID=c(1,2), cell=c(1,4)))
expect_equal(as.vector(cells(x, y, weights=TRUE)), c(1, 1, 1, 2, 2, 2, 1, 2, 4, 2, 3, 4, 0.521988, 0.414648, 0.043812, 0.000240, 0.072960, 0.491532))

expect_equivalent(unlist(extract(x, y, ID=TRUE)), c(1,2,1,4))
expect_equivalent(unlist(extract(x, y, cells=TRUE, weights=TRUE, ID=TRUE)), c(1, 1, 1, 2, 2, 2, 1, 2, 4, 2, 3, 4, 1, 2, 4, 2, 3, 4, 0.521988, 0.414648, 0.043812, 0.000240, 0.072960, 0.491532))


r <- rast(nrows=5, ncols=5, xmin=0, xmax=1, ymin=0, ymax=1, names="test")
r[c(2,7)] <- c(15, 20)
rr <- c(r, r/2)
names(rr)[2] <- "half"
xy <- cbind(x=0.3, y=c(0.9, 0.7))

v <- vect(xy)
e <- extract(r, v, ID=TRUE)
expect_equal(e, data.frame(ID=1:2, test=c(15,20)))
ee <- extract(rr, v, ID=TRUE)
expect_equal(ee, data.frame(ID=1:2, test=c(15,20), half=c(7.5, 10)))

e <- extract(r, v, cell=TRUE, ID=TRUE)
expect_equal(e, data.frame(ID=1:2, test=c(15,20), cell=c(2,7)))
ee <- extract(rr, v, cell=TRUE, ID=FALSE)
expect_equal(ee, data.frame(test=c(15,20), half=c(7.5, 10), cell=c(2,7)))

ee <- extract(rr, v, cell=TRUE, xy=TRUE, ID=TRUE)
expect_equal(ee, data.frame(ID=1:2, test=c(15,20), half=c(7.5, 10), cell=c(2,7), xy))

ee <- extract(rr, v, xy=TRUE, ID=TRUE)
expect_equal(ee, data.frame(ID=1:2, test=c(15,20), half=c(7.5, 10), xy))

f <- system.file("ex/meuse.tif", package="terra")
r <- rast(f)
xy <- cbind(179000, 330000)
xy <- rbind(xy-100, xy, xy+1000)
e <- extract(r, xy)
expect_equal(e[,1] , c(378, 251, 208))

vxy <- vect(xy)
e <- extract(r, vxy, ID=FALSE, method="bilinear")
expect_equal(e[,1] , c(378.00, 270.75, 197.25))

e <- extract(r, vxy, method="bilinear", cells=TRUE, ID=FALSE)
expect_equal(unlist(e, use.names=FALSE), c(378.00, 270.75, 197.25,8173.00,8016.00,6041.00))
 
r <- rast(nrows = 10, ncols = 10, nlyrs = 1, vals = 1:100, names = "temp")
x1 <- rbind(c(-145,-10), c(-145,-5), c(-140, -5), c(-140,-10))
x2 <- rbind(c(-10,0), c(140,60), c(160,0), c(140,-55))
z <- rbind(cbind(object=1, part=1, x1, hole=0), cbind(object=3, part=1, x2, hole=0))
colnames(z)[3:4] <- c('x', 'y')
p <- vect(z, "polygons", crs=crs(r))
rr <- c(r, r*2)
test <- terra::extract(r, p, ID=TRUE, fun = mean)
expect_equal(as.vector(as.matrix(test)), c(1,2,51.5,53))

test <- terra::extract(rr, p, fun = mean, ID=TRUE)
expect_equal(as.vector(as.matrix(test)), c(1,2,51.5,53,103,106))

test <- terra::extract(r, p, fun = mean, exact=TRUE, ID=TRUE)
expect_equal(round(as.vector(as.matrix(test)),5), c(1,2, 51.80006, 52.21312))

test <- terra::extract(rr, p, fun = mean, exact=TRUE, ID=TRUE)
expect_equal(round(as.vector(as.matrix(test)),5), c(1,2, 51.80006, 52.21312, 103.60012, 104.42623))

r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(100, 100, NaN, 100)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))


r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(100, NaN, 100, 100)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))

r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(100, 100, 100, NaN)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))


r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(NaN, 100, 100, 100)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))

r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(NaN, NaN, 100, 100)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))

r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(NaN, 100, NaN, 100)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))

r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(NaN, 100, 100, NaN)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))

r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(100, NaN, NaN, 100)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))

r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(100, NaN, 100, NaN)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))


r <- terra::rast(nrow = 2, ncol = 2, nlyrs = 1, xmin = -100, xmax = 100, ymin = -100, ymax = 100)
terra::values(r) <- c(100, 100, NaN, NaN)

expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 0.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = 45.0, Y = -45.0), method = "bilinear")$lyr.1), c(100.0))
expect_equal(round(terra::extract(r, data.frame(X = -45.0, Y = 45.0), method = "bilinear")$lyr.1), c(100.0))


# v <- vect(system.file("ex/lux.shp", package="terra"))[1:2,]
# z <- rast(v, resolution=.1, names="test")
# values(z) <- 1:ncell(z)
# ee <- extract(z, v, as.list=TRUE)
# expect_equal(rapply(ee, mean), c(5.20,16.75))


# x <- c(z, z*2, z/3)
# names(x) <- letters[1:3]
# ee <- extract(x, v, as.list=TRUE)
# m <- matrix(rapply(ee, mean), ncol=nlyr(x), byrow=TRUE)
# expect_equal(round(as.vector(t(m)),6), c(5.200000, 10.400000,  1.733333, 16.750000, 33.500000,  5.583333))

m = rast(nrow = 2, ncol = 2, vals=1:4)
expect_equal(cells(m), 1:4)


v <- vect(system.file("ex/lux.shp", package="terra"))
r <- rast(system.file("ex/elev.tif", package="terra"))

fun <- function(x, ...) c(min(x, ...), mean(x, ...), max(x, ...))
e <- extract(r, v, fun, na.rm=TRUE)[6,]
expect_equivalent(round(e, 4), c(6, 164, 314.9969, 403))

e <- extract(r, v, quantile, na.rm=TRUE)[12,]
expect_equivalent(round(e), c(12, 213, 274, 317, 352, 413))

fun <- c("min", "max")
e <- extract(r, v[1:2, ], fun, na.rm=TRUE)
expect_equivalent(unlist(e), c(1, 2, 339, 195, 547, 514))

x <- rast(ncol=4, nrow=4, nlyr=3, extent=ext(-2, 2, -2, 2))
values(x) <- rbind(11:13)
xy <- data.frame(x=-1:1,y=-1:1)
e <- extract(x, xy, ID=TRUE, layer=1:3)
expect_equal(e$value, 11:13)
