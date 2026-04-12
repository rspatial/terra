
## basic: returns a data.frame with expected columns
pp <- proj_pipelines("EPSG:4326", "EPSG:32632")

## first pipeline should be instantiable with 0 grids needed
expect_true(pp$instantiable[1])
expect_equal(pp$grid_count[1], 0)
expect_true(nchar(pp$definition[1]) > 0)

## accuracy == 0 for an exact same-datum transformation
#CI OSX fail
#expect_equal(pp$accuracy[1], 0)

## accepts SpatRaster input
r <- rast(ncols=10, nrows=10, crs="EPSG:4326")
pp2 <- proj_pipelines(r, "EPSG:3857")
expect_true(nrow(pp2) >= 1)

## grid_availability="DISCARD" reduces results
pp_all  <- proj_pipelines("EPSG:5070", "EPSG:4326")
pp_disc <- proj_pipelines("EPSG:5070", "EPSG:4326", grid_availability="DISCARD")
expect_true(nrow(pp_disc) <= nrow(pp_all))
expect_true(nrow(pp_disc) >= 1)

## AOI constrains results
pp_aoi <- proj_pipelines("EPSG:5070", "EPSG:4326", AOI=ext(-100, -90, 30, 40))
expect_true(nrow(pp_aoi) >= 1)
expect_true(nrow(pp_aoi) <= nrow(pp_all))

## invalid CRS
expect_error(proj_pipelines("INVALID", "EPSG:4326"))
expect_error(proj_pipelines("EPSG:4326", "INVALID"))
