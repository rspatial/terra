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

# tests for values=TRUE
m <- matrix(c(1, 1, 2, 2,
              1, 2, 2, 1,
              3, 3, 1, 1), nrow=3, byrow=TRUE)
r <- rast(m)

p4 <- patches(r, directions=4, values=TRUE)
expect_equal(as.vector(values(p4)), c(1, 1, 2, 2, 1, 2, 2, 4, 5, 5, 4, 4))

p8 <- patches(r, directions=8, values=TRUE)
expect_equal(as.vector(values(p8)), c(1, 1, 2, 2, 1, 2, 2, 3, 4, 4, 3, 3))

m_u <- matrix(c(1, 1, 1, 1,
                1, 2, 2, 1,
                1, 1, 1, 1), nrow=3, byrow=TRUE)
r_u <- rast(m_u)
p_u <- patches(r_u, directions=4, values=TRUE)
expect_equal(as.vector(values(p_u)), c(1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1))

r_large <- rast(nrows=100, ncols=100, vals=1)
r_large[50, 1:100] <- 2
p_large <- patches(r_large, directions=4, values=TRUE)
expect_equal(as.vector(unique(values(p_large))), c(1, 2, 3))

m_na <- matrix(c(1, 1, NA, 1,
                 1, NA, 1, 1,
                 1, 1, 1, 1), nrow=3, byrow=TRUE)
r_na <- rast(m_na)
p_na <- patches(r_na, directions=4, values=TRUE)
expect_equal(as.vector(values(p_na)), c(1, 1, NaN, 1, 1, NaN, 1, 1, 1, 1, 1, 1))

# zeroAsNA should be ignored if values=TRUE
r_z <- rast(matrix(c(1, 0, 0, 1), ncol=2))
p_z1 <- patches(r_z, directions=4, values=TRUE, zeroAsNA=TRUE)
expect_equal(as.vector(values(p_z1)), c(1, 2, 3, 4))
p_z2 <- patches(r_z, directions=4, values=TRUE, zeroAsNA=FALSE)
expect_equal(values(p_z1), values(p_z2))
