
m <- matrix(c(0,1,0,0,0,1), ncol=2)
v <- vect(m, type="polygons")
x <- as.vector(crds(v))

expect_equal(x, c(0,1,0,0,0,0,1,0))

