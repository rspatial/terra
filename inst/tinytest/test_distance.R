
x <- y <- seq(0, 3, 1)
d <- terra::distance(cbind(x, y), lonlat=FALSE, sequential=TRUE)
sq2 <- sqrt(2)

expect_equal(d, c(0, rep(sq2, 3)))

d <- terra::distance(cbind(x, y), lonlat=FALSE, sequential=FALSE)
expect_equal(as.vector(d), c(sq2, 2*sq2, 3*sq2, sq2, 2*sq2, sq2))

d <- terra::distance(cbind(x, y), lonlat=TRUE, sequential=TRUE)
expect_equivalent(d, c(0.0, 156899.6, 156876.1, 156829.3), .000001)

d <- terra::distance(cbind(x, y), lonlat=TRUE, sequential=FALSE)
expect_equivalent(as.vector(d), c(156899.6, 313775.7, 470605.0, 156876.1, 313705.4, 156829.3), .000001)
