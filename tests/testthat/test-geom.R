context("test-geom")

for (i in 1:4) {
	if (i == 1) {
		r <- rast(ncol = 20, nrow = 10, xmin=-10, xmax=10, ymin=-5, ymax=6)
		rx <- raster(ncol = 20, nrow = 10, xmn=-10, xmx=10, ymn=-5, ymx=6)
		values(r) <- ncell(r):1
		values(rx) <- ncell(rx):1
	} else if (i == 2) {
		r <- c(r, r/2, r*2)
		rx <- stack(rx, rx/2, rx*2)
	} else if (i == 3) {
		f <- system.file("exdata/logo.grd", package="terra")
		r <- rast(f)
		rx <- brick(f)
	} else {
		f <- system.file("exdata/test.tif", package="terra")
		r <- rast(f)
		rx <- raster(f)
	}
	
	test_that("raster and terra geometries work the same", { 
		expect_equal(res(r), res(rx)) 
		expect_equal(dim(r), dim(rx)) 
		expect_equivalent(as.vector(ext(r)), as.vector(extent(rx))) 
		expect_equal(xyFromCell(r, 10), xyFromCell(rx, 10))
		expect_equal(xyFromCell(r, 1), xyFromCell(rx, 1))
		expect_equal(xyFromCell(r, ncell(r)), xyFromCell(rx, ncell(rx)))
		expect_equal(as.vector(is.na(xyFromCell(r, 0))), c(TRUE, TRUE))
		expect_equal(as.vector(is.na(xyFromCell(r, ncell(r)+1))), c(TRUE, TRUE))
	} )


	test_that("raster and terra subsets work the same", { 
		expect_equivalent(r[10], rx[10]) 
		expect_equivalent(r[2,], rx[2,]) 
		expect_equivalent(r[,2], rx[,2]) 
		expect_equivalent(r[3:4, 2:3], rx[3:4, 2:3]) 
		expect_equivalent(r[[1]][5,6], rx[[1]][5,6]) 
	} )

}


