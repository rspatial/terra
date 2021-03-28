
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)
r[1:50,] <- NA
v <- c(717627, 304.33715, 141, 432, 308.22562, 48.40079, 141, 432)
f <- c("sum", "mean", "min", "max", "rms", "sd", "range")
x <- unlist(sapply(f, function(s) global(r, s, na.rm=TRUE)))
#terraOptions(steps=4, todisk=T)
y <- unlist(sapply(f, function(s) global(r, s, na.rm=TRUE)))
expect_equivalent(x,v)
expect_equal(x, y)
