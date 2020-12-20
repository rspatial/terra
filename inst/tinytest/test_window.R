
#library(terra)
x <- rast(system.file("ex/logo.tif", package="terra"))
y <- rast(system.file("ex/logo.tif", package="terra")) * 1
e <- ext(c(35,55,35,55))
z <- crop(x, e)

window(x) <- e
window(y) <- e
a <- c(z, y, x) 

s <- spatSample(x, 4, cell=TRUE)
expect_equal(s, spatSample(y, 3, cell=TRUE))
expect_equal(s, spatSample(z, 3, cell=TRUE))
expect_equal(s, spatSample(a, 3, cell=TRUE))

expect_equal(values(x), values(y))
expect_equal(values(x), values(z))

xy <- 10 * cbind(-1:6, -1:6)

e1 <- extract(x, xy)
e2 <- extract(y, xy)
e3 <- extract(z, xy)
e4 <- extract(a, xy)
e <- cbind(e1, e2, e3) 
expect_equal(e1, e2)
expect_equal(e1, e3)
expect_equal(e, e4)





expect_equal(dim(lux), c(12,5)) 
expect_equivalent(unlist(lux[,"ID_2",drop=TRUE]), c(1:7, 12, 8:11))	

expect_equal(size(x), 4) 
expect_equal(nrow(x), 4) 
expect_equal(ncol(x), 2)
expect_equal(unique(x$NAME_1), "Luxembourg")
expect_equivalent(unlist(x[,2,drop=TRUE]), 8:11)	

