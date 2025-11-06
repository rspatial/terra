# empty spatvectors

f <- system.file("ex", "lux.shp", package="terra")
v <- vect(f)
empty_v <- v[0,]

outfile <- file.path(tempdir(), "empty.gpkg")
writeVector(empty_v, outfile, overwrite=TRUE)

read_v <- vect(outfile)
expect_equal(nrow(read_v), 0)
expect_equal(ncol(read_v), ncol(empty_v))
expect_equal(names(read_v), names(empty_v))

expect_equal(nrow(vect(outfile, proxy=TRUE)), 0)
expect_equal(names(vect(outfile, proxy=TRUE)), names(empty_v))

unlink(outfile)


# insert=TRUE with empty layer
outfile <- file.path(tempdir(), "empty.gpkg")
writeVector(v, outfile, layer="data", overwrite=TRUE)
writeVector(empty_v, outfile, layer="empty", insert=TRUE)

expect_equal(nrow(vect(outfile, layer="data")), nrow(v))
expect_equal(nrow(vect(outfile, layer="empty")), 0)
expect_equal(names(vect(outfile, layer="empty")), names(empty_v))

unlink(outfile)
