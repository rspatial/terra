
# multiple sources
s <- rast(system.file("ex/logo.tif", package="terra"))
plotRGB(c(s[[1]], s[[2]], s[[3]]))

