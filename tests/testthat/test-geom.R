context("test-geom")

for (i in 1:5) {
	if (i == 1) {
		r <- rast(ncol = 20, nrow = 10, xmin=-10, xmax=10, ymin=-5, ymax=6)
		rx <- raster(ncol = 20, nrow = 10, xmn=-10, xmx=10, ymn=-5, ymx=6)
		values(r) <- ncell(r):1
		values(rx) <- ncell(rx):1
	} else if (i == 2) {
		r <- c(r, r/2, r*2)
		rx <- stack(rx, rx/2, rx*2)
	} else if (i == 3) {
		f <- system.file("exdata/test.tif", package="terra")
		r <- rast(f)
		rx <- raster(f)
	} else if (i == 4) {
		f <- system.file("exdata/logo.tif", package="terra")
		r <- rast(f)
		rx <- brick(f)
	} else {
		f <- system.file("exdata/logo.grd", package="terra")
		r <- rast(f)
		rx <- brick(f)
	}
	
	test_that("raster and terra geometries", { 
		expect_equal(res(r), res(rx)) 
		expect_equal(dim(r), dim(rx)) 
		expect_equivalent(as.vector(ext(r)), as.vector(extent(rx))) 
		expect_equal(xyFromCell(r, 10), xyFromCell(rx, 10))
		expect_equal(xyFromCell(r, 1), xyFromCell(rx, 1))
		expect_equal(xyFromCell(r, ncell(r)), xyFromCell(rx, ncell(rx)))
		expect_equal(as.vector(is.na(xyFromCell(r, 0))), c(TRUE, TRUE))
		expect_equal(as.vector(is.na(xyFromCell(r, ncell(r)+1))), c(TRUE, TRUE))
	} )


	test_that("raster and terra subsets", {
		a <- ifelse(i>2,20,0)
		if (i>3) {
			j <- 4999
			expect_equivalent(r[j],	cbind(r[[1]][j], r[[2]][j], r[[3]][j]))
		}
		expect_equivalent(r[10+a], rx[10+a]) 
		expect_equivalent(r[2+a,], rx[2+a,]) 
		expect_equivalent(r[,2+a], rx[,2+a]) 
		expect_equivalent(r[(3:4)+a, (2:3)+a], rx[(3:4)+a, (2:3)+a]) 
		expect_equivalent(r[[1]][5+a,6+a], rx[[1]][5+a,6+a]) 
	} )

}


