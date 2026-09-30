
r <- rast(matrix(0.5, 2, 2))
expect_equal(as.vector(values(2 - r)), rep(1.5, 4))
expect_equal(as.vector(values(r - 2)), rep(-1.5, 4))
expect_equal(as.vector(values(2 / r)), rep(4, 4))
expect_equal(as.vector(values(r / 2)), rep(0.25, 4))

r <- rast(nrows=3, ncols=3, nlyr=2, vals=1:18)
s <- scale(r)
expect_equal(attr(s, "scaled:center"), as.vector(global(r, "mean")[,1]))
rr <- r - attr(s, "scaled:center")
expect_equal(attr(s, "scaled:scale"), as.vector(global(rr, "rms")[,1]))
expect_equal(as.vector(values(s)), as.vector(values(rr / attr(s, "scaled:scale"))))

s <- scale(r, center=FALSE)
expect_null(attr(s, "scaled:center"))
expect_equal(attr(s, "scaled:scale"), as.vector(global(r, "rms")[,1]))

s <- scale(r, scale=FALSE)
expect_null(attr(s, "scaled:scale"))
expect_equal(attr(s, "scaled:center"), as.vector(global(r, "mean")[,1]))

s <- scale(r, center=c(1, 2), scale=10)
expect_equal(attr(s, "scaled:center"), c(1, 2))
expect_equal(attr(s, "scaled:scale"), 10)

