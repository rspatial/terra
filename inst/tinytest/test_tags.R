
x <- r <- rast(system.file("ex/elev.tif", package = "terra"))
metags(r) <- "FOO=BAR"
metags(r, layer = 1) <- "this=that"

metags(x) <- data.frame(name="FOO", value="BAR", domain="")

expect_equal(metags(r), metags(x))
expect_equal(metags(r, layer=1), data.frame(layer=1, name="this", value="that"))



