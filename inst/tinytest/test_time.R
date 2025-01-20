
x <- rast(nrows=2, ncols=2, vals=1:4, names="random")
time(x) <- as.POSIXct("2000-01-02 12:12:12", tz="America/Los_Angeles")
expect_equal(timeInfo(x)$zone, "America/Los_Angeles")
