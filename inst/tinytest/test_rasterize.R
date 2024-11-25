
v <- vect(system.file("ex/lux.shp", package = "terra"))
r <- rast(v, ncols = 75, nrows = 100)
z <- rasterize(v, r, cover = TRUE, by = "ID_2")
v <- unlist(z[30*75+28])
e <- c(0.01538462, NA, NA, NA, 0.9846154, NA, NA, NA, NA, NA, NA, NA)

expect_equivalent(v, e, tolerance=2e-07)
