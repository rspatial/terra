
f <- system.file("ex/lux.shp", package="terra")
lux <- vect(f)

listofraw <- geom(lux[1:2, ], wkb = TRUE)
wkb <- listofraw[[1]]

expect_equal(vect(listofraw, type = NULL, crs = crs(lux)), lux[1:2, 0])

hex <- geom(lux[1, ], hex = TRUE)

expect_equal(typeof(listofraw), "list") 
expect_equal(typeof(wkb), "raw")

expect_equal(tolower(paste0(as.character(wkb), collapse = "")), tolower(hex))

