
f <- system.file("ex/lux.shp", package="terra")
lux <- vect(f)

listofraw <- geom(lux[1:2, ], wkb = TRUE)
wkb <- listofraw[[1]]
hex <- geom(lux[1, ], hex = TRUE)


expect_equal(typeof(listofraw), "list") 
expect_equal(typeof(wkb), "raw")

expect_equal(tolower(paste0(as.character(wkb), collapse = "")), tolower(hex))

