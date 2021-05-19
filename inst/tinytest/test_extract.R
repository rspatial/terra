
#f <- system.file("ex/lux.shp", package="terra")
#y <- vect(f)
#y <- y[1:2,]
#x <- rast(y, res=.2)
#values(x) <- 1:ncell(x)
#expect_equal(cells(x, y), cbind(id=c(1,2), cell=c(1,5)))
#expect_equal(as.vector(cells(x, y, weights=TRUE)), c(1, 1, 1, 2, 2, 2, 1, 2, 5, 4, 5, 6, 0.55, 0.41, 0.04, 0.08, 0.51, 0.11))

#expect_equivalent(unlist(extract(x, y)), c(1,2,1,5))
#expect_equivalent(unlist(extract(x, y, cells=TRUE, weights=TRUE)), c(1, 1, 1, 2, 2, 2, 1, 2, 5, 4, 5, 6, 1, 2, 5, 4, 5, 6, 0.55, 0.41, 0.04, 0.08, 0.51, 0.11))


r <- rast(nrows=5, ncols=5, xmin=0, xmax=1, ymin=0, ymax=1, names="test")
r[c(2,7)] <- c(15, 20)
rr <- c(r, r/2)
names(rr)[2] <- "half"
xy <- cbind(x=0.3, y=c(0.9, 0.7))
v <- vect(xy)
e <- extract(r, v)
expect_equal(e, data.frame(ID=1:2, test=c(15,20)))
ee <- extract(rr, v)
expect_equal(ee, data.frame(ID=1:2, test=c(15,20), half=c(7.5, 10)))

e <- extract(r, v, cell=TRUE)
expect_equal(e, data.frame(ID=1:2, test=c(15,20), cell=c(2,7)))
ee <- extract(rr, v, cell=TRUE)
expect_equal(ee, data.frame(ID=1:2, test=c(15,20), half=c(7.5, 10), cell=c(2,7)))

ee <- extract(rr, v, cell=TRUE, xy=TRUE)
expect_equal(ee, data.frame(ID=1:2, test=c(15,20), half=c(7.5, 10), cell=c(2,7), xy))

ee <- extract(rr, v, xy=TRUE)
expect_equal(ee, data.frame(ID=1:2, test=c(15,20), half=c(7.5, 10), xy))

