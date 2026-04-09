
r <- rast(nrows=10, ncols=10, xmin=0)
values(r) <- 1:ncell(r)

w <- matrix(c(1,1,1,1,0,1,1,1,1), 3)

m1 <- autocor(r, w)
m2 <- autocor(r, w, standardize=TRUE)

expect_equal(as.numeric(m1), 0.8362573, tolerance=1e-6)
expect_equal(as.numeric(m2), 0.9330909, tolerance=1e-6)
expect_true(m2 >= -1 && m2 <= 1)

g1 <- autocor(r, w, method="geary")
g2 <- autocor(r, w, method="geary", standardize=TRUE)

expect_equal(as.numeric(g1), 0.04421053, tolerance=1e-6)
expect_equal(as.numeric(g2), 0.04384, tolerance=1e-5)

# rook case -- Moran's I
wr <- matrix(c(0,1,0,1,0,1,0,1,0), 3)
mr1 <- autocor(r, wr)
mr2 <- autocor(r, wr, standardize=TRUE)

expect_equal(as.numeric(mr1), 0.8888889, tolerance=1e-6)
expect_equal(as.numeric(mr2), 0.96, tolerance=1e-6)
expect_true(mr2 >= -1 && mr2 <= 1)

if (require(spdep, quietly=TRUE)) {

	nb <- cell2nb(nrow=10, ncol=10, type="queen")
	v <- values(r)[,1]
	n <- length(v)

	lw_B <- nb2listw(nb, style="B")
	lw_W <- nb2listw(nb, style="W")

	m_B <- moran(v, lw_B, n, Szero(lw_B))$I
	m_W <- moran(v, lw_W, n, Szero(lw_W))$I

	expect_equal(as.numeric(m1), m_B, tolerance=1e-6)
	expect_equal(as.numeric(m2), m_W, tolerance=1e-6)

	g_B <- geary(v, lw_B, n, n - 1, Szero(lw_B))$C
	g_W <- geary(v, lw_W, n, n - 1, Szero(lw_W))$C

	expect_equal(as.numeric(g1), g_B, tolerance=1e-6)
	expect_equal(as.numeric(g2), g_W, tolerance=1e-6)

	nb_r <- cell2nb(nrow=10, ncol=10, type="rook")
	lw_rB <- nb2listw(nb_r, style="B")
	lw_rW <- nb2listw(nb_r, style="W")

	mr_B <- moran(v, lw_rB, n, Szero(lw_rB))$I
	mr_W <- moran(v, lw_rW, n, Szero(lw_rW))$I

	expect_equal(as.numeric(mr1), mr_B, tolerance=1e-6)
	expect_equal(as.numeric(mr2), mr_W, tolerance=1e-6)
}
