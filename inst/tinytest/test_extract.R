
f <- system.file("ex/lux.shp", package="terra")
y <- vect(f)
y <- y[1:2,]
x <- rast(y, res=.2)
values(x) <- 1:ncell(x)
expect_equal(cells(x, y), cbind(id=1, cell=1))
expect_equal(as.vector(cells(x, y, weights=TRUE)), c(1,1,1,1,2,2,1,2,3,4,3,4,0.6,0.19,0.02,0.01,0.16,0.39))
expect_equivalent(unlist(extract(x, y)), c(1,2,1,NA))
expect_equivalent(unlist(extract(x, y, cells=TRUE, weights=TRUE)), c(1,1,1,1,2,2,1,2,3,4,3,4,1,2,3,4,3,4,0.6,0.19,0.02,0.01,0.16,0.39))
