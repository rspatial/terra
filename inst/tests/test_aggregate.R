
r <- rast(ncol = 4, nrow = 4, xmin=0, xmax=1, ymin=0, ymax=1)
values(r) <- c(NA, 2:16)

expect_equal(as.vector(values(aggregate(r, 2, mean))), c(NaN, 5.5, 11.5, 13.5))
expect_equal(as.vector(values(aggregate(r, 2, mean, na.rm=TRUE))), c(4 + 1/3, 5.5, 11.5, 13.5))

rr = c(r,r*2)
expect_equal(as.vector(values(aggregate(rr, 2, mean))), c(NaN, 5.5, 11.5, 13.5, NaN, 11.0, 23.0, 27.0))
expect_equal(as.vector(values(aggregate(rr, 2, mean, na.rm=TRUE))), c(4 + 1/3, 5.5, 11.5, 13.5, 8+2/3, 11.0, 23.0, 27.0))
expect_equal(as.vector(values(aggregate(rr, 2, min, na.rm=TRUE))), c(2, 3, 9, 11, 4, 6, 18, 22))


