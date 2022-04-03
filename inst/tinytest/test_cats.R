
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE) - 1
r0 <- r * 1

# categories in a single layer
lv <- c("forest", "water", "urban")
levels(r) <- lv
names(r) <- "land cover"
v <- cats(r)[[1]]
#coltab(r) <- rainbow(4)

expect_equal(v$value, 0:2)
expect_equal(v$category, lv)

# reading/writing categories 
ftmp <- tempfile(fileext = ".tif")
z <- writeRaster(r, ftmp, overwrite=TRUE)
v <- cats(z)[[1]]
expect_equal(v$value, 0:2)
expect_equal(v$`land cover`, lv)

# reading/writing subset of categories 
levels(r) = cats(r)[[1]][2:3,]
zz = writeRaster(r, ftmp, overwrite=TRUE)
v <- cats(zz)[[1]]

expect_equal(v$value, 1:2)
expect_equal(v$category, lv[-1])

# categories in multiple layers
r2 <- rast(list(a=r0, b=r0))

# avoid warnings related to assumed lon/lat 
crs(r2) <- "EPSG:4326"

# values are numeric initially
expect_equal(values(r2)[1,], c(a = 1, b = 1))

# set layer a categories
levels(r2) <- lv
expect_equal(levels(r2), list(lv, ""))

# verify samples for layer a have categories
set.seed(0)
expect_equal(spatSample(r2, 2), data.frame(a = factor(c("forest", "forest"), levels = lv), 
                                           b = c(0, 0)))

# set all layer categories
levels(r2) <- rep(list(lv), 2)
expect_equal(levels(r2), list(lv, lv))

# verify samples for layer a and layer b have categories
set.seed(1)
expect_equal(spatSample(r2, 2), data.frame(a = factor(c("forest", "water"), levels = lv),
                                           b = factor(c("forest", "water"), levels = lv)))

# make sure no errors when show()ing factors
expect_silent(show(r2))
