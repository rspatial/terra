
r1 <- terra::rast(xmin = 0, xmax = 1); r1$x = 1; r1$y = 2
r2 <- terra::rast(xmin = 1, xmax = 2); r2$x = 3; r2$y = 4
expect_equal(names(terra::merge(r1, r2)), c("x", "y"))

x <- rast(xmin=-110, xmax=-50, ymin=40, ymax=70, ncols=30, nrows=15)
y <- rast(xmin=-80, xmax=-20, ymax=60, ymin=30)
res(y) <- res(x)

set.seed(1)
vx <- 1:ncell(x)
vx[floor(runif(50, min = 300, max = 450))] <- NA
vy <- 1:ncell(y)
vy[floor(runif(50, min = 0, max = 150))] <- NA
values(x) <- vx
values(y) <- vy

m <- terra::merge(x, y)
a <- m[11,1:20][[1]]
w <- which(is.na(a))
expect_equal(w, c(2,3,9))
expect_equal(a[-w], c(301, 304, 305, 306, 307, 308, 310, 311, 312, 313, 314, 315, 151, 317, 153, 319, 320))


x <- rast(xmin=-110, xmax=-60, ymin=40, ymax=70, res=1, vals=1)
y <- rast(xmin=-95, xmax=-45, ymax=60, ymin=30, res=1, vals=2)
z <- rast(xmin=-80, xmax=-30, ymax=50, ymin=20, res=1, vals=3)

m0 <- mosaic(x, y, z, fun="mean")
d <- diag(as.matrix(m0, wide=T))
expect_equal(d, c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3))


m1 <- mosaic(x, y, z, fun="blend")
d <- diag(as.matrix(m1, wide=T))

expect_equal(d, c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.03333333333333, 1.1, 1.16666666666667, 1.23333333333333, 1.3, 1.36666666666667, 1.43333333333333, 1.5, 1.56666666666667, 1.63333333333333, 1.7, 1.76666666666667, 1.83333333333333, 1.88461538461538, 1.95454545454545, 2.05, 2.15, 2.25, 2.35, 2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3))
