
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE) - 1
lv <- c("forest", "water", "urban")
levels(r) <- lv
names(r) <- "land cover"
v <- cats(r)[[1]]
#coltab(r) <- rainbow(4)

expect_equal(v$value, 0:2)
expect_equal(v$category, lv)

ftmp <- tempfile(fileext = ".tif")
z <- writeRaster(r, ftmp, overwrite=TRUE)
v <- cats(z)[[1]]
expect_equal(v$value, 0:2)
expect_equal(v$`land cover`, lv)

levels(r) = cats(r)[[1]][2:3,]
zz = writeRaster(r, ftmp, overwrite=TRUE)
v <- cats(zz)[[1]]

expect_equal(v$value, 1:2)
expect_equal(v$category, lv[-1])
