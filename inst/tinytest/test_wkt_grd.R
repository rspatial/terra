fl <- system.file("ex/test.grd", package="terra")
tst <- terra::rast(fl)
terra::crs(tst) <- "EPSG:28992"
tf <- tempfile(fileext=".grd")
terra::writeRaster(tst, tf)
tst1 <- terra::rast(tf)
expect_identical(terra::crs(tst, describe=TRUE)[1:3], terra::crs(tst1, describe=TRUE)[1:3])

