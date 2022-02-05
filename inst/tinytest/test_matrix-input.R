
expect_equal(class(rast(volcano))[1], "SpatRaster")
## degeneracies
dm1 <- rast(matrix(1:10))
dm2 <- rast(matrix(1:10, nrow = 1))
expect_equal(class(dm1)[1], "SpatRaster")
expect_equal(class(dm2)[1], "SpatRaster")

expect_equal(dim(dm1), c(10L, 1L, 1L))
expect_equal(dim(dm2), c(1L, 10L, 1L))

r <- rast(matrix(1))
expect_equal(class(r)[1], "SpatRaster")

expect_error(rast(matrix(1:4)[0, , drop = FALSE]))

# test_that matrix input type xyz works
r <- rast(cbind(as.matrix(expand.grid(1:3, 1:4)), 1:12), type = "xyz")
expect_equal(class(r)[1], "SpatRaster")

m <- cbind(x=c(1,1,2,2), y=c(1,2,1,2), z1=1:4, z2=5:8)
r <- rast(m, type="xyz")
expect_equal(as.vector(values(r)), c(2,4,1,3,6,8,5,7))

#----------------------------------------------------------------
## ensure matrix and array input result in extent 0,nrow,0,ncol
## no input
r <- rast(ncol = 4, nrow = 4, xmin=0, xmax=1, ymin=0, ymax=1)
expect_equal(c(xmin(r), xmax(r), ymin(r), ymax(r)), c(0, 1, 0, 1))

## matrix
r <- rast(matrix(1:12, 3))
expect_equal(c(xmin(r), xmax(r), ymin(r), ymax(r)), c(0, 4, 0, 3))
## array (degenerate)
r <- rast(array(1:12, c(3, 4)))
expect_equal(c(xmin(r), xmax(r), ymin(r), ymax(r)), c(0, 4, 0, 3))
## array
r <- rast(array(1:24, c(3, 4, 2)))
expect_equal(c(xmin(r), xmax(r), ymin(r), ymax(r)), c(0, 4, 0, 3))

## unless we provide on creation
r <- rast(matrix(1:12, 3), extent = c(-2, 10, -1, 20))
expect_equal(c(xmin(r), xmax(r), ymin(r), ymax(r)), c(-2, 10, -1, 20))
r <- rast(array(1:12, c(3, 4)), extent = c(-2, 10, -1, 20))
expect_equal(c(xmin(r), xmax(r), ymin(r), ymax(r)), c(-2, 10, -1, 20))

r <- rast(array(1:24, c(3, 4, 2)), extent = c(-2, 10, -1, 20))
expect_equal(c(xmin(r), xmax(r), ymin(r), ymax(r)), c(-2, 10, -1, 20))


#----------------------------------------------------------------
