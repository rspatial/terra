
r1 <- terra::rast(xmin = 0, xmax = 1); r1$x = 1; r1$y = 2
r2 <- terra::rast(xmin = 1, xmax = 2); r2$x = 3; r2$y = 4
expect_equal(names(terra::merge(r1, r2)), c("x", "y"))
