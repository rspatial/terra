
#https://github.com/rspatial/terra/issues/1410

r <- terra::rast(ncols=100, nrows=100)
values(r) <- 1:ncell(r)
r2 <- r*2
file_name_1 <- paste0(tempfile(), ".tif")
file_name_2 <- paste0(tempfile(), ".tif")
vrt_filename <- file.path(tempdir(),"test.vrt")
terra::writeRaster(r, file_name_1)
terra::writeRaster(r2, file_name_2)
x <- terra::vrt(c(file_name_1,file_name_2), filename = vrt_filename, overwrite=TRUE, options="-separate")

expect_equal(nlyr(x), 2)

