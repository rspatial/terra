
f <- system.file("ex/logo.tif", package = "terra")
r <- rast(f)
km <- k_means(r, centers=5)
r[1] <- NA
km <- k_means(r, centers=5)
km <- k_means(r, centers=5, maxcell=10)


