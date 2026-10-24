# intersection of disjoint envelopes is NULL
a <- ext(c(xmin = 0, xmax = 10, ymin = 0, ymax = 10))
b <- ext(c(xmin = 100, xmax = 101, ymin = 100, ymax = 101))
expect_null(intersect(a, b))

# shift(..., to=) recenters on an extent, object, or x, y
e1 <- ext(0, 2, 0, 4)
e2 <- ext(10, 14, 20, 26)
expect_equivalent(as.vector(shift(e1, 5, 5)), c(5, 7, 5, 9))
expect_equivalent(as.vector(shift(e1, to=e2)), c(11, 13, 21, 25))
expect_equivalent(as.vector(shift(e1, to=c(12, 23))), c(11, 13, 21, 25))
expect_equivalent(as.vector(shift(e1, to=c(12, NA))), c(11, 13, 0, 4))
expect_equivalent(as.vector(shift(e1, to=c(NA, 23))), c(0, 2, 21, 25))
expect_equivalent(as.vector(shift(e1, to=c(10, 14, 20, 26))), c(11, 13, 21, 25))

r1 <- rast(e1, ncols=2, nrows=2, vals=1:4)
r2 <- rast(e2, ncols=2, nrows=2)
expect_equivalent(as.vector(ext(shift(r1, to=r2))), c(11, 13, 21, 25))
expect_equivalent(as.vector(values(shift(r1, to=r2))), 1:4)

p1 <- as.polygons(e1)
p2 <- as.polygons(e2)
expect_equivalent(as.vector(ext(shift(p1, to=p2))), c(11, 13, 21, 25))
expect_equivalent(as.vector(ext(shift(p1, to=r2))), c(11, 13, 21, 25))

expect_error(shift(e1, dx=1, to=e2), pattern="both")
expect_error(shift(e1, to=12), pattern="to")
expect_error(shift(e1, to="nope"), pattern="to")