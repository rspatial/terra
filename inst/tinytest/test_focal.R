
m <- matrix(0,3,3)
m[c(4,6)] <- c(1,-1)
r <- rast(nrows=3, ncols=3, vals=1:9, crs="+proj=merc")
f <- focal(r, m)
e <- c(-4, -5, -6, -6, -6, -6, 4, 5, 6)
expect_equal(e, as.vector(values(f)))

f <- focal(r, t(m))
e <- c(-2, -2, 2, -5, -2, 5, -8, -2, 8)
expect_equal(e, as.vector(values(f)))

f <- focal(r, m, fun=mean)
e <- c(-4, -5, -6, -3, -3, -3, 4, 5, 6)
expect_equal(e, as.vector(values(f)))

m <- matrix(1,3,3)
f <- focal(r, m)
e <- c(12,21,16,27,45,33,24,39,28)
expect_equal(e, as.vector(values(f)))

f <- focal(r, 3)
expect_equal(e, as.vector(values(f)))

f <- focal(r, m, na.rm=FALSE, fillvalue=0)
expect_equal(e, as.vector(values(f)))

f <- focal(r, m, na.rm=FALSE)
e <- c(NA, NA, NA, NA, 45, NA, NA, NA, NA)
expect_equal(e, as.vector(values(f)))

f <- focal(r, 3, na.rm=FALSE)
expect_equal(e, as.vector(values(f)))

