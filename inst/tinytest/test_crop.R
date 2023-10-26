
r <- rast(nrow=10, ncol=10, extent=ext(0, 10, 0, 10))

e <- ext(-10, 15, -10, 15)
expect_equivalent(as.vector(ext(crop(r, e))), c(0,10,0,10))
expect_equivalent(as.vector(ext(crop(r, e, extend=TRUE))), c(-10,15,-10,15))

e <- ext(-10, 5, -10, 5)
expect_equivalent(as.vector(ext(crop(r, e))), c(0,5,0,5))
expect_equivalent(as.vector(ext(crop(r, e, extend=TRUE))), c(-10,5,-10,5))

e <- ext(0, 5, 0, 5)
expect_equivalent(as.vector(ext(crop(r, e))), c(0,5,0,5))
expect_equivalent(as.vector(ext(crop(r, e, extend=TRUE))), c(0,5,0,5))


values(r) <- 1
e <- ext(-10, 15, -10, 15)
expect_equivalent(as.vector(ext(crop(r, e))), c(0,10,0,10))
expect_equivalent(as.vector(ext(crop(r, e, extend=TRUE))), c(-10,15,-10,15))

e <- ext(-10, 5, -10, 5)
expect_equivalent(as.vector(ext(crop(r, e))), c(0,5,0,5))
expect_equivalent(as.vector(ext(crop(r, e, extend=TRUE))), c(-10,5,-10,5))

e <- ext(0, 5, 0, 5)
expect_equivalent(as.vector(ext(crop(r, e))), c(0,5,0,5))
expect_equivalent(as.vector(ext(crop(r, e, extend=TRUE))), c(0,5,0,5))

