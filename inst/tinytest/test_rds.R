
f <- system.file("ex/rds_tst.rds", package="terra")
x <- rast(f)

y <- rast(nrows=2, ncols=2, vals=1:4, names="random")
time(y) <- as.POSIXct("2025-01-25 21:28:51", tz="UTC") 

expect_true(all.equal(x, y))

