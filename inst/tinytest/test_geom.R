
r <- rast(ncol = 20, nrow = 10, xmin=-10, xmax=10, ymin=-5, ymax=6)
values(r) <- ncell(r):1
r <- c(r, r/2, r*2)
	
# compare geometries
expect_equal(res(r), c(1.0, 1.1)) 
expect_equal(dim(r), c(10, 20, 3)) 
expect_equivalent(as.vector(ext(r)), c(-10, 10, -5, 6)) 
expect_equal(xyFromCell(r, 10),  cbind(x=-0.5, y=5.45))
expect_equal(xyFromCell(r, 1), cbind(x=-9.5, y=5.45))

expect_equal(xyFromCell(r, ncell(r)), cbind(x=9.5, y=-4.45))
expect_equal(as.vector(is.na(xyFromCell(r, 0))), c(TRUE, TRUE))
expect_equal(as.vector(is.na(xyFromCell(r, ncell(r)+1))), c(TRUE, TRUE))

expect_equivalent(r[10], data.frame(cbind(191, 95.5, 382)))
expect_equivalent(r[2,], data.frame(matrix(c(180,179,178,177,176,175,174,173,172,171,170,169,168,167,166,165,164,163,162,161,90,89.5,89,88.5,88,87.5,87,86.5,86,85.5,85,84.5,84,83.5,83,82.5,82,81.5,81,80.5,360,358,356,354,352,350,348,346,344,342,340,338,336,334,332,330,328,326,324,322), ncol=3)))
expect_equivalent(r[,2], data.frame(matrix(c(199,179,159,139,119,99,79,59,39,19,99.5,89.5,79.5,69.5,59.5,49.5,39.5,29.5,19.5,9.5,398,358,318,278,238,198,158,118,78,38), ncol=3)))
expect_equivalent(r[3:4, 2:3], data.frame(matrix(c(159,158,139,138,79.5,79,69.5,69,318,316,278,276),ncol=3)))
expect_equivalent(r[[1]][5,6][1], data.frame(lyr.1=115))

