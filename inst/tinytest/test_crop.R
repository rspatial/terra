
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

# INT2U without NoData: 65535 stays data when crop writes to disk (#2207)
set.seed(1)
v <- matrix(sample(1000:60000, 500, replace=TRUE), ncol=5)
v[c(12, 47, 83), ] <- 65535
ru <- rast(nrows=10, ncols=10, nlyrs=5, extent=ext(0, 10, 0, 10), crs="", vals=v)
f2207 <- tempfile(fileext=".tif")
writeRaster(ru, f2207, datatype="INT2U", NAflag=NA, overwrite=TRUE)
xu <- rast(f2207)
eu <- ext(0, 8, 0, 10)
terraOptions(todisk=FALSE)
a2207 <- crop(xu, eu)
terraOptions(todisk=TRUE)
b2207 <- crop(xu, eu)
terraOptions(todisk=FALSE)
expect_equal(global(a2207, "mean")[,1], global(b2207, "mean")[,1])
expect_false(anyNA(values(a2207)))
expect_false(anyNA(values(b2207)))
expect_true(any(values(b2207) == 65535))
unlink(f2207)

