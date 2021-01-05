
#f <- system.file("ex/lux.shp", package="terra")
#y <- vect(f)
#y <- y[1:2,]
#x <- rast(y, res=.2)
#values(x) <- 1:ncell(x)
#expect_equal(cells(x, y), cbind(id=c(1,2), cell=c(1,5)))
#expect_equal(as.vector(cells(x, y, weights=TRUE)), c(1, 1, 1, 2, 2, 2, 1, 2, 5, 4, 5, 6, 0.55, 0.41, 0.04, 0.08, 0.51, 0.11))

#expect_equivalent(unlist(extract(x, y)), c(1,2,1,5))
#expect_equivalent(unlist(extract(x, y, cells=TRUE, weights=TRUE)), c(1, 1, 1, 2, 2, 2, 1, 2, 5, 4, 5, 6, 1, 2, 5, 4, 5, 6, 0.55, 0.41, 0.04, 0.08, 0.51, 0.11))

