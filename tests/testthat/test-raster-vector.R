#context("test-raster-vector")
#
#p1 <- rbind(c(-180, -20), c(-140, 55), c(10, 0), c(-140, -60), c(-180, -20))
#hole <- rbind(c(-150, -20), c(-100, -10), c(-110, 20), c(-150, -20))
#p2 <- rbind(c(-10, 0), c(140, 60), c(160, 0), c(140, -55), c(-10, 0))
#p3 <- rbind(c(-125, 0), c(0, 60), c(40, 5), c(15, -45), c(-125, 0))
#z <- rbind(cbind(object=1, part=1, p1, hole=0), 
#		   cbind(object=1, part=1, hole, hole=1),
#           cbind(object=2, part=1, p2, hole=0),
#		   cbind(object=3, part=1, p3, hole=0))   
#colnames(z)[3:4] <- c('x', 'y')
#z <- data.frame(z)
#v <- vect(z, "polygons")
#vx <- as(z, "SpatialPolygons")
#crs(vx) <- "+proj=longlat +datum=WGS84"
#
#r <- rast(ncol = 20, nrow = 10)
#rx <- raster(ncol = 20, nrow = 10)
#rv <- rasterize(v, r)
#rvx <- rasterize(vx, rx)
#
#test_that("rasterize", {
#	expect_equivalent(values(rv), values(rvx)) 
#})
#
#test_that("crop", {
#	#ex <- extent(-90,90,-60,60)
#	#e <- ext(-90,90,-60,60)
#	rc <- crop(rv, v)
#	rcx <- crop(rvx, vx)
#	expect_equivalent(values(rc), values(rcx)) 
#})
#
#
#test_that("mask", {
#	values(r) <- values(rx) <- 1:ncell(r)
#	m <- mask(r, v)
#	mx <- mask(rx, vx)
#	expect_equivalent(values(m), values(mx)) 
#})
#
#
#set.seed(0)
#values(r) <- values(rx) <- runif(ncell(r))
#
#
#test_that("classify", {
#	rcl <- cbind(from=seq(0,0.9,.1), to=seq(0.1,1,.1), becomes=1:10) 
#	rc  <- classify(r, rcl)
#	rcx <- reclassify(rx, rcl)
#	expect_equivalent(values(rc), values(rcx)) 
#})
#
#
#test_that("app", {
#	s <- c(r, log(r), sqrt(r))
#	sx <- stack(rx, log(rx), sqrt(rx))
#	rd  <- app(s, mean)
#	rdx  <- calc(sx, mean)
#	expect_equivalent(values(rd), values(rdx), tolerance=0.00001) 
#})
#
#