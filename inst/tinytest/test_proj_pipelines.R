
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


### 
a <- rast(ncols=40, nrows=40, xmin=-120, xmax=-60, ymin=20, ymax=60, crs="EPSG:4326")
values(a) <- 1:ncell(a)
b <- rast(ncols=94, nrows=124, ext=c(-2562000, 3768000, -364000, 4470000), crs="EPSG:5070")

pp <- proj_pipelines("EPSG:4326", "EPSG:5070")
w1 <- project(a, b, pipe=pp[1,])
w2 <- project(a, b, pipe=pp[2,])
w3 <- project(a, crs(b), pipe=pp[2,])
w4 <- project(a, pipe=pp[2,])

expect_true(all.equal(w1, w2))
expect_true(all.equal(w3, w4))


library(terra)
r <- rast(ncols=40, nrows=40, xmin=-80, xmax=-56.5, ymin=45, ymax=62, crs="EPSG:4326")
values(r) <- 1:ncell(r)
vr <- as.polygons(r)
p1 <- project(vr, pipe=pp[1,])
p2 <- project(vr, pipe=pp[2,])
p3 <- project(vr, pipe=pp[3,])

p1b <- project(vr, "EPSG:5070", pipe=pp[1,])
p3b <- project(vr, "EPSG:5070", pipe=pp[3,])

