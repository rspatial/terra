# Verify SpatRaster$focal2 produces identical output to $focal
# for every code path (plain box, separable, generic, NA policies).

.spatOpts <- function() terra:::spatOptions("", FALSE, list())

f1 <- function(x, w, m, fillvalue, narm, naonly, naomit, fun, expand) {
	r <- x
	r@pntr <- x@pntr$focal(w, m, fillvalue, narm, naonly, naomit, fun, expand, .spatOpts())
	as.vector(values(r))
}
f2 <- function(x, w, m, fillvalue, narm, naonly, naomit, fun, expand) {
	r <- x
	r@pntr <- x@pntr$focal2(w, m, fillvalue, narm, naonly, naomit, fun, expand, .spatOpts())
	as.vector(values(r))
}

cmp <- function(x, w, m, fun, ..., fillvalue=NA, narm=FALSE, naonly=FALSE, naomit=FALSE, expand=FALSE) {
	a <- f1(x, w, m, fillvalue, narm, naonly, naomit, fun, expand)
	b <- f2(x, w, m, fillvalue, narm, naonly, naomit, fun, expand)
	expect_equal(a, b, info=paste("fun=", fun, "narm=", narm, "expand=", expand, ...))
}


# ── small 3x3 raster with the original tinytest scenarios ─────────────────────

set.seed(1)

m <- matrix(0,3,3)
m[c(4,6)] <- c(1,-1)
r <- rast(nrows=3, ncols=3, vals=1:9, crs="+proj=merc")

# rank-1 weighted (m has positive and negative weights, but is rank-1)
cmp(r, c(3,3), as.vector(t(m)), "sum", narm=TRUE)
cmp(r, c(3,3), as.vector(t(m)), "sum", narm=FALSE, fillvalue=0)
cmp(r, c(3,3), as.vector(t(t(m))), "sum", narm=TRUE)

# plain uniform
m1 <- matrix(1,3,3)
cmp(r, c(3,3), as.vector(t(m1)), "sum", narm=TRUE)
cmp(r, c(3,3), as.vector(t(m1)), "sum", narm=FALSE)
cmp(r, c(3,3), as.vector(t(m1)), "sum", narm=FALSE, fillvalue=0)
cmp(r, c(3,3), as.vector(t(m1)), "mean", narm=TRUE)
cmp(r, c(3,3), as.vector(t(m1)), "mean", narm=FALSE)


# ── lon/lat 3x3 (boundary should fill with NA, not wrap because not global) ──

r <- rast(nrow=3, ncol=3, vals=1:9)
m1 <- as.vector(t(matrix(1,3,3)))
cmp(r, c(3,3), m1, "sum",  narm=TRUE)
cmp(r, c(3,3), m1, "sum",  narm=FALSE)
cmp(r, c(3,3), m1, "mean", narm=TRUE)
cmp(r, c(3,3), m1, "mean", narm=FALSE)
cmp(r, c(3,3), m1, "max",  narm=TRUE)


# ── global raster: column wrap ───────────────────────────────────────────────

rg <- rast(nrow=3, ncol=4, vals=1:12)        # default ext is global lonlat
cmp(rg, c(3,3), m1, "sum",  narm=TRUE)
cmp(rg, c(3,3), m1, "sum",  narm=FALSE)
cmp(rg, c(3,3), m1, "mean", narm=TRUE)


# ── expand boundary ──────────────────────────────────────────────────────────

r <- rast(nrows=4, ncols=5, vals=1:20, crs="+proj=merc")
cmp(r, c(3,3), m1, "sum",  narm=FALSE, expand=TRUE)
cmp(r, c(3,3), m1, "mean", narm=FALSE, expand=TRUE)
cmp(r, c(3,3), m1, "sum",  narm=TRUE,  expand=TRUE)


# ── larger raster, larger window, with NAs and different policies ────────────

r <- rast(ncols=20, nrows=15, vals=1:300, crs="+proj=merc")
v <- values(r); v[c(7, 31, 122, 200, 250)] <- NA; values(r) <- v

m_plain <- rep(1, 25)                                  # 5x5 uniform
m_gauss <- as.vector(outer(c(1,2,4,2,1)/10, c(1,2,4,2,1)/10))   # rank-1 weighted

# plain box
for (narm in c(TRUE, FALSE)) {
	for (fillvalue in list(NA, 0)) {
		cmp(r, c(5,5), m_plain, "sum",  narm=narm, fillvalue=fillvalue)
		cmp(r, c(5,5), m_plain, "mean", narm=narm, fillvalue=fillvalue)
	}
}
# separable rank-1
for (narm in c(TRUE, FALSE)) {
	cmp(r, c(5,5), m_gauss, "sum",  narm=narm, fillvalue=0)
	cmp(r, c(5,5), m_gauss, "mean", narm=narm, fillvalue=0)
}
# expand boundary, plain
cmp(r, c(5,5), m_plain, "sum",  narm=TRUE,  expand=TRUE)
cmp(r, c(5,5), m_plain, "mean", narm=FALSE, expand=TRUE)

# na.policy="only" / "omit" pass-through
cmp(r, c(5,5), m_plain, "sum",  narm=TRUE, naonly=TRUE)
cmp(r, c(5,5), m_plain, "mean", narm=TRUE, naomit=TRUE)

# generic fallback (max) — plain window but non-sum/mean
cmp(r, c(5,5), m_plain, "max",  narm=TRUE)
cmp(r, c(5,5), m_plain, "modal", narm=TRUE)


# ── multi-layer ──────────────────────────────────────────────────────────────

r2 <- c(r, r * 2, r + 5)
cmp(r2, c(5,5), m_plain, "sum",  narm=TRUE)
cmp(r2, c(5,5), m_plain, "mean", narm=FALSE, fillvalue=0)
cmp(r2, c(5,5), m_gauss, "sum",  narm=TRUE,  fillvalue=0)
