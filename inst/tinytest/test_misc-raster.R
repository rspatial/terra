# Value checks for raster spatial-context and distance helpers (see man/*.Rd).

## adjacent — rook neighbours of centre cell in a 1x5 raster
ra <- rast(nrows = 1, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 1, vals = 1:5)
adj <- adjacent(ra, cells = 3, directions = "rook", pairs = TRUE)
expect_equivalent(sort(adj[, "to"]), c(2L, 4L))

## boundaries — edges between two classes (see ?boundaries)
rb <- rast(nrows = 1, ncols = 5, vals = c(1, 1, 2, 2, 2))
expect_equivalent(as.vector(boundaries(rb, classes = TRUE)[]), c(0, 1, 1, 0, 0))

## patches — 2x2 clump of 1s on NA background (four cells same patch id)
rp <- rast(nrows = 4, ncols = 4, vals = NA_real_)
rp[2:3, 2:3] <- 1
pp <- patches(rp)
expect_equal(sum(pp[] == 1, na.rm = TRUE), 4)

## sieve — isolated single cell removed (see ?sieve)
rsv <- rast(nrows = 18, ncols = 18, xmin = 0, vals = 0, crs = "local")
rsv[2, 5] <- 1
expect_equal(as.numeric(sieve(rsv, 8)[2, 5]), 0)

## cellSize + expanse (elev.tif layer)
f <- system.file("ex/elev.tif", package = "terra")
r <- rast(f)
r <- if (nlyr(r) > 1) r[[1]] else r
cs <- cellSize(r)
expect_true(all(cs[] > 0, na.rm = TRUE))
expect_true(sum(unlist(expanse(r))) > 0)

## surfArea — projected synthetic (lon/lat not supported)
rpj <- rast(
  nrows = 10, ncols = 10,
  crs = "+proj=utm +zone=32 +datum=WGS84",
  xmin = 0, xmax = 1000, ymin = 0, ymax = 1000,
  vals = runif(100)
)
sa <- surfArea(rpj)
expect_true(all(sa[] > 0, na.rm = TRUE))


## direction — ?direction example
rd <- rast(ncol = 36, nrow = 18, crs = "+proj=merc")
values(rd) <- NA
rd[306] <- 1
di <- direction(rd, degrees = TRUE)
expect_true(all(di[] >= 0 & di[] <= 360, na.rm = TRUE))

## shade — hillshade in [0, 1]
slope <- terrain(r, v = "slope", unit = "radians")
asp <- terrain(r, v = "aspect", unit = "radians")
sh <- shade(slope, asp)
expect_true(min(sh[], na.rm = TRUE) >= 0 && max(sh[], na.rm = TRUE) <= 1)

## global Moran's I in [-1, 1]
aci <- autocor(r, global = TRUE)
expect_true(all(aci >= -1 & aci <= 1, na.rm = TRUE))

## where.min / where.max — global cell for single layer (see ?where)
r9 <- rast(nrows = 3, ncols = 3, vals = 1:9)
wm <- where.min(r9)
wx <- where.max(r9)
expect_equal(as.numeric(wm[1, "cell"]), 1)
expect_equal(as.numeric(wm[1, "value"]), 1)
expect_equal(as.numeric(wx[1, "cell"]), 9)
expect_equal(as.numeric(wx[1, "value"]), 9)

## NIDP — neighbour counts from flow directions
fd <- terrain(r, v = "flowdir")
nid <- NIDP(fd)
expect_true(all(nid[] >= 0 & nid[] <= 9, na.rm = TRUE))
