# Value checks for local / cell-based raster functions (man/*.Rd).

r <- rast(nrows = 9, ncols = 9, xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "local")
values(r) <- 1:ncell(r)
v <- values(r)[, 1]

## clamp
cl <- clamp(r, lower = 10, upper = 50)
expect_equal(as.vector(values(cl)), pmin(pmax(v, 10), 50))

## subst
sb <- subst(r, from = c(1, 2), to = c(NA, NA))
vv <- v
vv[v %in% c(1, 2)] <- NA
expect_equal(as.vector(values(sb)), vv)

## init — cell index
ini <- init(r, fun = "cell")
expect_equal(as.vector(values(ini)), as.numeric(1:ncell(r)))

## segregate — three binary layers for r < 20
sg <- segregate(r < 20)
expect_equal(nlyr(sg), 2L)
expect_equal(sum(values(sg)[]), ncell(r))

## cover — first layer wins where not NA
r2 <- r * 2
cv <- cover(r, r2)
expect_equal(as.vector(values(cv)), v)

## roll — constant stack, rolling mean equals constant
r1 <- rast(nrows = 3, ncols = 3, crs = "local", vals = 1)
stk <- c(r1, r1, r1, r1, r1)
rl <- roll(stk, 3, "mean", circular = TRUE)
expect_true(max(abs(values(rl)[] - 1), na.rm = TRUE) < 1e-10)

## diff — second minus first layer
expect_equal(as.vector(values(diff(c(r, r * 2)))), v)

## modal — identical layers
md <- modal(c(r, r, r))
expect_equal(as.vector(values(md)), v)

## thresh — mean split vs manual classify
tr <- thresh(r, method = "mean", as.raster = TRUE)
mu <- unlist(global(r, "mean", na.rm = TRUE))[1]
br <- ifel(r > mu, 1, 0)
expect_equal(as.numeric(values(tr)), as.numeric(values(br)))

## scale_linear — [0, 1] range
scl <- scale_linear(r, min = 0, max = 1)
mn <- min(v)
mx <- max(v)
expect_equal(as.vector(values(scl)), (v - mn) / (mx - mn))

## focalMat — circle weights (d in map units relative to resolution)
fm <- focalMat(r, d = 3 * res(r)[1], type = "circle")
expect_equal(dim(fm), c(7L, 7L))
expect_true(sum(fm) > 0)

## which.lyr — first layer where condition is true
wly <- which.lyr(c(r > 0, r > 10000))
expect_true(all(as.vector(values(wly)) == 1, na.rm = TRUE))

## selectRange / extractRange
idx <- classify(r, rbind(c(-Inf, 40, 1), c(40, Inf, 2)))
sr <- selectRange(c(r, r * 2, r * 3), idx)
exp_sr <- ifelse(values(idx)[, 1] == 1, v, 2 * v)
expect_equal(as.vector(values(sr)), exp_sr)
xy <- data.frame(x = c(0.15, 0.5), y = c(0.15, 0.5))
stk <- c(r, r * 2, r * 3)
er <- extractRange(stk, xy, first = 1, last = 3)
ex <- extract(stk, xy)
expect_equal(as.numeric(unlist(er[[1]][1, ])), as.numeric(ex[1, -1]))

## selectHighest — five largest cells marked
sh <- selectHighest(r, n = 5)
expect_equal(sum(sh[] == 1, na.rm = TRUE), 5)

## classify — two-column lookup (from test_classify.R pattern)
set.seed(1)
rcin <- rast(nrows = 10, ncols = 10, vals = sample(1:5, 100, replace = TRUE))
rcl <- data.frame(from = 1:3, to = 11:13)
rco <- classify(rcin, rcl)
vi <- as.vector(values(rcin))
vo <- as.vector(values(rco))
expect_equal(vo[vi == 1], rep(11, sum(vi == 1)))
expect_equal(vo[vi == 2], rep(12, sum(vi == 2)))
expect_equal(vo[vi == 3], rep(13, sum(vi == 3)))
expect_equal(vo[vi > 3], vi[vi > 3])

## not.na / as.bool
expect_equal(sum(not.na(r)[]), ncell(r))
expect_equal(sum(as.bool(r > 40)[]), sum(v > 40))

## divide — total cell value mass conserved across zones
set.seed(0)
dv <- divide(r, n = 2, as.raster = TRUE)
zs <- zonal(r, dv, sum, na.rm = TRUE)
expect_equal(sum(zs[, 2]), sum(v))
expect_true(nrow(zs) >= 2)

## approximate — NA filled from neighbouring layers
ra <- rast(ncols = 3, nrows = 3, crs = "local", vals = 1:9)
rb <- ra * 1.1
ra[5] <- NA
ap <- approximate(c(ra, rb))
expect_false(is.na(values(ap)[5, 1]))

## sum across layers (per-cell)
rs <- app(c(r, r2), sum)
expect_equal(as.vector(values(rs)), v + 2 * v)
