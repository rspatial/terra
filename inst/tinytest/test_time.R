
x <- rast(nrows=2, ncols=2, vals=1:4, names="random")
time(x) <- as.POSIXct("2000-01-02 12:12:12", tz="America/Los_Angeles")
expect_equal(timeInfo(x)$zone, "America/Los_Angeles")

# time from any source is kept if step/zone match or can be unified (#2160)
a <- rast(nrows=2, ncols=2, nlyr=1, vals=1)
b <- rast(nrows=2, ncols=2, nlyr=1, vals=2)
time(a) <- as.Date("2000-01-01")
ab <- c(a, b)
expect_true(has.time(ab))
expect_equal(as.Date(time(ab)[1]), as.Date("2000-01-01"))
expect_true(is.na(time(ab)[2]))
ba <- c(b, a)
expect_true(has.time(ba))
expect_true(is.na(time(ba)[1]))
expect_equal(as.Date(time(ba)[2]), as.Date("2000-01-01"))

# Calendar steps combine at the coarsest shared type (#2160)
time(b) <- as.POSIXct("2000-01-02 12:00:00", tz="UTC")
sd <- c(a, b)
expect_true(has.time(sd))
expect_equal(timeInfo(sd)$step, "days")

time(b, tstep="years") <- 2001
yd <- c(a, b)
expect_true(has.time(yd))
expect_equal(timeInfo(yd)$step, "years")
expect_equal(time(yd), c(2000, 2001))

time(b, tstep="yearmonths") <- as.Date("2001-06-15")
ym <- c(a, b)
expect_true(has.time(ym))
expect_equal(timeInfo(ym)$step, "yearmonths")

# months (1-12) do not share a calendar axis with days
time(b, tstep="months") <- 6
md <- suppressWarnings(c(a, b))
expect_false(has.time(md))
