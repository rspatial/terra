
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
