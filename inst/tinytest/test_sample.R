# Value checks for spatSample() with different methods and options (see ?spatSample).

set.seed(42)
r <- rast(
  nrows = 10, ncols = 10,
  xmin = 0, xmax = 1, ymin = 0, ymax = 1,
  crs = "local", vals = 1:100
)
nc <- ncell(r)
vv <- values(r)[, 1]

## method = "random" — row count, valid cells, no duplicates without replacement
sr <- spatSample(r, 20, method = "random", cells = TRUE)
expect_equal(nrow(sr), 20)
expect_true(all(sr[, "cell"] %in% seq_len(nc)))
expect_equal(length(unique(sr[, "cell"])), 20)

## method = "random", replace = TRUE — sample size can exceed unique cells
set.seed(1)
r5 <- rast(nrows = 5, ncols = 5, xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "local", vals = 1:25)
srep <- spatSample(r5, 30, method = "random", replace = TRUE, cells = TRUE)
expect_equal(nrow(srep), 30)
expect_true(length(unique(srep[, "cell"])) < nrow(srep))

## method = "random", na.rm = TRUE — only non-NA cells
rna <- r
rna[1:70] <- NA
sn <- spatSample(rna, 8, method = "random", na.rm = TRUE, cells = TRUE)
expect_equal(nrow(sn), 8)
expect_true(all(!is.na(vv[sn[, "cell"]])))

## method = "regular" — coordinates inside extent, values match source raster
sreg <- spatSample(r, 9, method = "regular", cells = TRUE, xy = TRUE, values = TRUE)
expect_true(all(sreg[, "x"] >= 0 & sreg[, "x"] <= 1))
expect_true(all(sreg[, "y"] >= 0 & sreg[, "y"] <= 1))
expect_equal(sreg[, "lyr.1"], vv[sreg[, "cell"]])

## method = "regular", as.raster = TRUE — coarser grid, same extent
srast <- spatSample(r, 50, method = "regular", as.raster = TRUE)
expect_true(ncell(srast) < nc)
expect_equal(as.vector(ext(srast)), as.vector(ext(r)))

## method = "regular", exact = TRUE — requested number of rows (projected crs)
sreg25 <- spatSample(r, 25, method = "regular", cells = TRUE, exact = TRUE)
expect_equal(nrow(sreg25), 25)

## method = "regular" — size as c(nrow, ncol) grid dimensions
sgrid <- spatSample(r, c(4, 5), method = "regular", cells = TRUE)
expect_equal(nrow(sgrid), 20)

## method = "stratified", each = TRUE — size per stratum
set.seed(1)
rr <- rast(
  nrows = 10, ncols = 10,
  xmin = 0, xmax = 1, ymin = 0, ymax = 1,
  crs = "local"
)
values(rr) <- sample(1:3, ncell(rr), replace = TRUE)
st1 <- spatSample(rr, 3, method = "stratified", cells = TRUE, values = TRUE, each = TRUE)
expect_equal(nrow(st1), 9)
tb1 <- table(st1[, "lyr.1"])
expect_equal(as.vector(tb1), c(3L, 3L, 3L))

## method = "stratified", each = FALSE — total sample size across strata
st2 <- spatSample(rr, 9, method = "stratified", cells = TRUE, values = TRUE, each = FALSE)
expect_equal(nrow(st2), 9)

## method = "weights" — draws follow cell index raster (positive weights)
rw <- init(r, "cell") / nc
sw <- spatSample(rw, 24, method = "weights", cells = TRUE)
expect_equal(nrow(sw), 24)
expect_true(all(sw[, "cell"] %in% seq_len(nc)))

## method = "spread" — approximately uniform coverage, requested count
set.seed(0)
ssp <- spatSample(r, 12, method = "spread", cells = TRUE)
expect_equal(nrow(ssp), 12)
expect_true(all(ssp[, "cell"] %in% seq_len(nc)))

## ext — samples restricted to sub-extent
sube <- ext(0, 0.5, 0, 0.5)
sext <- spatSample(r, 20, method = "random", cells = TRUE, ext = sube)
xye <- xyFromCell(r, sext[, "cell"])
expect_true(all(xye[, 1] <= 0.5 + 1e-9))
expect_true(all(xye[, 2] <= 0.5 + 1e-9))

## SpatExtent — random points in bounding box
e <- ext(r)
sm <- spatSample(e, 12, method = "random", lonlat = FALSE)
expect_equal(nrow(sm), 12)
expect_true(all(sm[, "x"] >= as.vector(e)[1] & sm[, "x"] <= as.vector(e)[2]))
expect_true(all(sm[, "y"] >= as.vector(e)[3] & sm[, "y"] <= as.vector(e)[4]))

## SpatVector — random points inside data extent
f <- system.file("ex/lux.shp", package = "terra")
v <- vect(f)
set.seed(2)
pv <- spatSample(v, 6)
expect_equal(nrow(pv), 6)
ev <- as.vector(ext(v))
cr <- crds(pv)
expect_true(all(cr[, 1] >= ev[1] & cr[, 1] <= ev[2]))
expect_true(all(cr[, 2] >= ev[3] & cr[, 2] <= ev[4]))
