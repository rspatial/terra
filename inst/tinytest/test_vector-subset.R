
f <- system.file("ex/lux.shp", package="terra")
lux <- vect(f)

expect_equal(dim(lux), c(12,6)) 
expect_equivalent(unlist(lux[,"ID_2",drop=TRUE]), c(1:7, 12, 8:11))	

x <- lux[lux$NAME_1 == "Luxembourg", 2:3]

expect_equal(length(x), 4) 
expect_equal(nrow(x), 4) 
expect_equal(ncol(x), 2)
expect_equal(unique(x$NAME_1), "Luxembourg")
expect_equivalent(unlist(x[,2,drop=TRUE]), 8:11)	

lux$ID_1 <- factor(LETTERS[1:12])
expect_equal(class(lux$ID_1), "factor")
