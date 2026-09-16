
## points: fun="modal"
rm <- rast(xmin=0, ncols=3, nrows=3)
xy <- rbind(
	xyFromCell(rm, c(1, 1, 1)),
	xyFromCell(rm, c(2, 2, 2, 2)),
	xyFromCell(rm, 5)
)
vals <- c(1, 1, 2,  5, 3, 5, 5,  7)
mm <- rasterize(xy, rm, values=vals, fun="modal")
vm <- as.vector(values(mm))
expect_equal(vm[1], 1)
expect_equal(vm[2], 5)
expect_equal(vm[5], 7)
expect_true(is.na(vm[3]))
# ties: equal counts -> lowest value
xy2 <- xyFromCell(rm, c(1, 1, 1, 1))
mt <- rasterize(xy2, rm, values=c(2, 2, 1, 1), fun="modal")
expect_equal(as.vector(values(mt))[1], 1)
# value= is an alias for values=
set.seed(1)
r <- rast(xmin=0, ncols=18, nrows=18)
p <- spatSample(r, 1000, xy=TRUE, replace=TRUE)
y <- rasterize(p, r, value=sample(10, 1000, replace=TRUE), fun="modal")
expect_false(all(is.na(values(y))))
expect_true(all(values(y) %in% c(1:10, NA)))

v <- vect(system.file("ex/lux.shp", package = "terra"))
r <- rast(v, ncols = 75, nrows = 100)
z <- rasterize(v, r, cover = TRUE, by = "ID_2")
v <- unlist(z[30*75+28])
e <- c(0.01538462, NA, NA, NA, 0.9846154, NA, NA, NA, NA, NA, NA, NA)

expect_equivalent(v, e, tolerance=2e-07)
