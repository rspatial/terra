
library(terra)

# subst integer lookup
r <- rast(nrows=10, ncols=10, vals=sample(1:5, 100, replace=TRUE))
r_sub <- subst(r, 1:2, 21:22)
v_sub <- as.vector(values(r_sub))
v_in <- as.vector(values(r))
expect_equal(v_sub[v_in == 1], rep(21, sum(v_in == 1)))
expect_equal(v_sub[v_in == 2], rep(22, sum(v_in == 2)))
expect_equal(v_sub[v_in > 2], v_in[v_in > 2])

# subst floating point lookup
r <- rast(nrows=2, ncols=2, vals=c(0.1, 0.2, 0.5, 0.9))
r_sub <- subst(r, c(0.1, 0.5), c(99.1, 99.5))
v_sub <- as.vector(values(r_sub))
v_in <- as.vector(values(r))
expect_equal(v_sub[v_in == 0.1], 99.1)
expect_equal(v_sub[v_in == 0.5], 99.5)
expect_equal(v_sub[v_in == 0.2], 0.2)
expect_equal(v_sub[v_in == 0.9], 0.9)

# subst 'others' with non-integers
r <- rast(nrows=2, ncols=2, vals=c(1, 2, 3, 4))
r_sub <- subst(r, 1:2, 11:12, others=1.5)
v_sub <- as.vector(values(r_sub))
expect_equal(v_sub, c(11, 12, 1.5, 1.5))
