
r <- rast(matrix(c(0,1,1,0), ncol=2))
p <- patches(r, directions = 8, zeroAsNA=TRUE)
expect_equal(as.vector(values(p)), c(NaN, 1, 1, NaN))

p <- patches(r, directions = 4, zeroAsNA=TRUE)
expect_equal(as.vector(values(p)), c(NaN, 1, 2, NaN))

r <- rast(matrix(c(1,0,0,1), ncol=2))
p <- patches(r, directions = 8, zeroAsNA=TRUE)
expect_equal(as.vector(values(p)), c(1, NaN, NaN, 1))

p <- patches(r, directions = 4, zeroAsNA=TRUE)
expect_equal(as.vector(values(p)), c(1, NaN, NaN, 2))

p <- patches(r, directions = 4, zeroAsNA=FALSE)
expect_equal(as.vector(values(p)), c(1, 1, 1, 1))


r <- rast(nrows=18, ncols=36)
r[1:2, 5:8] <- 11
r[7:8, 1:6] <- 12
r[5:6, 22:36] <- 13
r[15:16, 18:29] <- 14
p <- patches(r)
expect_equal(as.vector(unique(values(p))), c(NaN, 1:4))
p <- patches(r, directions=8)
expect_equal(as.vector(unique(values(p))), c(NaN, 1:3))

xmin(r) <- 0
p <- patches(r, directions=8)
expect_equal(as.vector(unique(values(p))), c(NaN, 1:4))

