
r <- rast(system.file("ex/elev.tif", package = "terra"))
metags(r) <- "FOO=BAR"
metags(r, layer = 1) <- "this=that"

expect_equal(metags(r), data.frame(name=c("AREA_OR_POINT", "FOO"), value=c("Area", "BAR"), domain=""))
expect_equal(metags(r, layer=1), data.frame(layer=1, name="this", value="that"))



