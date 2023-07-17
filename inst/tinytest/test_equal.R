
x <- sqrt(1:100)
mat <- matrix(x, 10, 10)
r1 <- rast(nrows=10, ncols=10, xmin=0, vals = x)
r2 <- rast(nrows=10, ncols=10, xmin=0, vals = mat)

expect_inherits(all.equal(r1, r2), "character")

expect_true(all.equal(r1, r1))

# works for multiple layers
r3 <- c(r1, r2)
r4 <- c(r1, r2)

expect_true(all.equal(r3, r4))

