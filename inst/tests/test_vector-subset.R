
f <- system.file("exdata/lux.shp", package="terra")
lux <- vect(f)

expect_equal(dim(lux), c(12,5)) 
expect_equivalent(unlist(lux[,"ID_2",drop=TRUE]), c(1:7, 12, 8:11))	

x <- lux[lux$NAME_1 == "Luxembourg", 2:3]

expect_equal(size(x), 4) 
expect_equal(nrow(x), 4) 
expect_equal(ncol(x), 2)
expect_equal(unique(x$NAME_1), "Luxembourg")
expect_equivalent(unlist(x[,2,drop=TRUE]), 8:11)	

