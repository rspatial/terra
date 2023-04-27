
v <- rast(matrix(c(1, 2, 3, NA), ncol = 2))
z <- rast(matrix(c(1, 1, 2, 2), ncol = 2))

a <- zonal(v, z, \(x) sum(is.na(x)))
expect_equal(unlist(a, use.names=FALSE), c(1,2,0,1))

a <- zonal(v, z, "isNA")
expect_equal(unlist(a, use.names=FALSE), c(1,2,0,1))

a <- zonal(v, z, \(x) mean(x))
expect_equal(unlist(a, use.names=FALSE), c(1,2,1.5,NA))

a <- zonal(v, z, "mean")
expect_equal(unlist(a, use.names=FALSE), c(1,2,1.5,NA))

a <- zonal(v, z, \(x) mean(x, na.rm=T))
expect_equal(unlist(a, use.names=FALSE), c(1,2,1.5,3))

a <- zonal(v, z, \(x, ...) mean(x, ...), na.rm=T)
expect_equal(unlist(a, use.names=FALSE), c(1,2,1.5,3))

a <- zonal(v, z, "mean", na.rm=T)
expect_equal(unlist(a, use.names=FALSE), c(1,2,1.5,3))
