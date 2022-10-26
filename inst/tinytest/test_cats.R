
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE)
r0 <- r * 1

# categories in a single layer
lv <- data.frame(id=1:4, cover=c("forest", "water", "urban", "other"), stringsAsFactors=FALSE)
levels(r) <- lv
names(r) <- "land cover"
v <- cats(r)[[1]]
#coltab(r) <- rainbow(4)

expect_equal(v$id, lv[,1])
expect_equal(v$cover, lv[,2])

r <- droplevels(r)
expect_equal(levels(r)[[1]][,2], c("forest", "water", "urban"))
lv <- lv[-4,]

# reading/writing categories 
ftmp <- tempfile(fileext = ".tif")
z <- writeRaster(r, ftmp, overwrite=TRUE)
v <- cats(z)[[1]]
expect_equal(v$value, 1:3)
expect_equal(v$`land cover`, lv[,2])

# reading/writing subset of categories 
levels(r) = cats(r)[[1]][2:3,]
zz = writeRaster(r, ftmp, overwrite=TRUE)
v <- cats(zz)[[1]]

expect_equal(v$value, 2:3)
expect_equal(v$cover, lv[-1,2])

# categories in multiple layers
r2 <- rast(list(a=r0, b=r0))

# values are numeric initially
expect_equal(values(r2)[1,], c(a = 2, b = 2))

# set layer a categories
levels(r2) <- lv
expect_equal(levels(r2), list(lv, ""))

# verify samples for layer a have categories
expect_equal( r2[c(15,3)], data.frame(cover = factor(c("forest", "urban"), levels = lv[,2]), 
                                           b = c(1, 3), stringsAsFactors=FALSE ))

# set all layer categories
levels(r2) <- rep(list(lv), 2)
expect_equal(levels(r2), list(lv, lv))

# verify samples for layer a and layer b have categories
expect_equal(r2[c(15,3)], data.frame(cover = factor(c("forest", "urban"), levels = lv[,2]),
                                           cover = factor(c("forest", "urban"), levels = lv[,2]), stringsAsFactors=FALSE, check.names = FALSE))

# make sure no errors when show()ing factors
expect_silent(show(r2))
