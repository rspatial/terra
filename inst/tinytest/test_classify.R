
r <- rast(ncols=3, nrows=3)
set.seed(68)
values(r) <- runif(ncell(r)) # from example
# all values > 0 and <= 0.25 become 1, etc.
m <- c(0, 0.25, 1,  0.25, 0.5, 2,  0.5, 1, 3) # from example
rclmat <- matrix(m, ncol=3, byrow=TRUE)

rc <- classify(r, rclmat, include.lowest=TRUE)
expect_equal(as.vector(values(rc)), c(3, 3, 2, 3, 2, 3, 3, 1, 2))

values(r) <- (1:ncell(r))-1
m <- c(0, 2, 1,  2, 3, 2,  3, 8, 3) 
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- classify(r, rclmat, right=FALSE)
expect_equal(as.vector(values(rc)), c(1, 1, 2, 3, 3, 3, 3, 3, 8))

rc <- classify(r, rclmat, right=TRUE)
expect_equal(as.vector(values(rc)), c(0, 1, 1, 2, 3, 3, 3, 3, 3))

rc <- classify(r, rclmat, right=TRUE, include.lowest=TRUE)
expect_equal(as.vector(values(rc)), c(1, 1, 1, 2, 3, 3, 3, 3, 3))

 
