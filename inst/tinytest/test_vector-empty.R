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

outfile <- file.path(tempdir(), "schema_mismatch.gpkg")

a <- vect("POLYGON ((0 0, 0 1, 1 1, 1 0, 0 0))")
b <- vect("POLYGON ((10 10, 10 11, 11 11, 11 10, 10 10))")
a$id <- 1L
b$id <- "text"

suppressWarnings({
  empty_intersect <- intersect(a, b)
})

expect_equal(nrow(empty_intersect), 0)
expect_equal(names(empty_intersect), c("id_1", "id_2"))

writeVector(empty_intersect, outfile, overwrite=TRUE)

read_result <- vect(outfile)
expect_equal(nrow(read_result), 0)
expect_equal(names(read_result), c("id_1", "id_2"))

proxy_result <- vect(outfile, proxy=TRUE)
expect_equal(nrow(proxy_result), 0)
expect_equal(names(proxy_result), c("id_1", "id_2"))

unlink(outfile)
