
r <- rast(nrow=2, ncol=2, vals=1:4)
x <- c(r, r)
x[1] <- cbind(5, 6)
x[2] <- c(7, 8)
expect_equal(as.vector(values(x)), c(5,7,3,4,6,8,3,4))



