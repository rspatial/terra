
m <- matrix(0,3,3)
m[c(4,6)] <- c(1,-1)
r <- rast(nrows=3, ncols=3, vals=1:9, crs="+proj=merc")
f <- focal(r, m, na.rm=TRUE)
e <- c(-4, -5, -6, -6, -6, -6, 4, 5, 6)
expect_equal(e, as.vector(values(f)))

f <- focal(r, t(m), na.rm=TRUE)
e <- c(-2, -2, 2, -5, -2, 5, -8, -2, 8)
expect_equal(e, as.vector(values(f)))

#expect_error(focal(r, m, fun=mean, na.rm=TRUE), pattern="[focal]")

#e <- c(-4, -5, -6, -3, -3, -3, 4, 5, 6)
#expect_equal(e, as.vector(values(f)))

m <- matrix(1,3,3)
f <- focal(r, m, na.rm=TRUE)
e <- c(12,21,16,27,45,33,24,39,28)
expect_equal(e, as.vector(values(f)))

f <- focal(r, 3, na.rm=TRUE)
expect_equal(e, as.vector(values(f)))

f <- focal(r, m, na.rm=FALSE, fillvalue=0)
expect_equal(e, as.vector(values(f)))

f <- focal(r, m, na.rm=FALSE)
e <- c(NA, NA, NA, NA, 45, NA, NA, NA, NA)
expect_equal(e, as.vector(values(f)))

f <- focal(r, 3, na.rm=FALSE)
expect_equal(e, as.vector(values(f)))



r <- rast(nrow=3, ncol=3)
values(r) <- 1:ncell(r)
rr <- rast(nrow=3, ncol=3, xmin=0)
values(rr) <- 1:ncell(rr)

f <- focalValues(r)[1,] 
e <- c(NA, NA, NA, 3, 1, 2, 6, 4, 5)
expect_equal(e, f)

f <- focalValues(rr)[1,]
e <- c(NA, NA, NA, NA, 1, 2, NA, 4, 5)
expect_equal(e, f)

f <- as.vector(values(focal(rr, 3, max, na.rm=TRUE)))
e <- c(5, 6, 6, 8, 9, 9, 8, 9, 9)
expect_equal(e, f)

f <- as.vector(values(focal(r, 3, max, na.rm=TRUE)))
e <- c(6, 6, 6, 9, 9, 9, 9, 9, 9)
expect_equal(e, f)

f <- as.vector(values(focal(rr, 3, sum, na.rm=TRUE)))
e <- c(12, 21, 16, 27, 45, 33, 24, 39, 28)
expect_equal(e, f)

f <- as.vector(values(focal(r, 3, sum, na.rm=TRUE)))
e <- c(21, 21, 21, 45, 45, 45, 39, 39, 39)
expect_equal(e, f)

f <- as.vector(values(focal(rr, 3, mean, na.rm=FALSE)))
e <- c(NA, NA, NA, NA, 5, NA, NA, NA, NA)
expect_equal(e, f)

f <- as.vector(values(focal(r, 3, mean, na.rm=FALSE)))
e <- c(NA, NA, NA, 5, 5, 5, NA, NA, NA)
expect_equal(e, f)


r <- rast(ncols=100, nrows=100, ext(0, 10, 0, 10))
values(r) = 1:ncell(r)
r[5,]=NA
f= focal(r, w=5, fun=mean, na.policy="only", na.rm=TRUE, wopt=list(steps=4))
x  = (f - r)
expect_equal(sum(values(x), na.rm=TRUE), 0)

