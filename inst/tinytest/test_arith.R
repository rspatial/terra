
r <- rast(matrix(0.5, 2, 2))
expect_equal(as.vector(values(2 - r)), rep(1.5, 4))
expect_equal(as.vector(values(r - 2)), rep(-1.5, 4))
expect_equal(as.vector(values(2 / r)), rep(4, 4))
expect_equal(as.vector(values(r / 2)), rep(0.25, 4))

