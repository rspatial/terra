# Validate that the TBB-parallel paths in classify() (lookup_apply +
# reclass_vector_range) and subst() (lookup_apply + replaceValues) produce the
# same results as the sequential fallback. The C++ kernels only switch to the
# parallel path when the per-block work is >= 16384 cells, so we use a 200x200
# raster (40000 cells) to make sure the TBB branch is actually exercised.

set.seed(1)
nr <- 200L; nc <- 200L
v_in <- sample(1L:6L, nr*nc, replace=TRUE)
r <- rast(nrows=nr, ncols=nc, vals=v_in)


run_with_parallel <- function(par, threads, expr) {
	old <- terraOptions(print=FALSE)
	on.exit(terraOptions(parallel = old$parallel, threads = old$threads), add=TRUE)
	terraOptions(parallel = par, threads = threads)
	force(eval(expr, parent.frame()))
}


# ---- subst (lookup_apply path) ----

ref_sub <- run_with_parallel(FALSE, 0L,
	quote(subst(r, 1:3, 101:103)))
par_sub <- run_with_parallel(TRUE, 4L,
	quote(subst(r, 1:3, 101:103)))
expect_equal(values(par_sub), values(ref_sub))


# subst with 'others' (forces the alternate branch in lookup_apply)
ref_sub2 <- run_with_parallel(FALSE, 0L,
	quote(subst(r, 1:3, 101:103, others=-1)))
par_sub2 <- run_with_parallel(TRUE, 4L,
	quote(subst(r, 1:3, 101:103, others=-1)))
expect_equal(values(par_sub2), values(ref_sub2))


# ---- classify is-becomes (lookup_apply path) ----

rcl_ib <- matrix(c(1,11, 2,12, 3,13), ncol=2, byrow=TRUE)
ref_cl_ib <- run_with_parallel(FALSE, 0L,
	quote(classify(r, rcl_ib, others=NA)))
par_cl_ib <- run_with_parallel(TRUE, 4L,
	quote(classify(r, rcl_ib, others=NA)))
expect_equal(values(par_cl_ib), values(ref_cl_ib))


# ---- classify from-to-becomes (reclass_vector_range path) ----

rcl_ftb <- matrix(c(0,2,1, 2,5,2, 5,8,3), ncol=3, byrow=TRUE)
ref_cl_ftb <- run_with_parallel(FALSE, 0L,
	quote(classify(r, rcl_ftb, include.lowest=TRUE, right=TRUE)))
par_cl_ftb <- run_with_parallel(TRUE, 4L,
	quote(classify(r, rcl_ftb, include.lowest=TRUE, right=TRUE)))
expect_equal(values(par_cl_ftb), values(ref_cl_ftb))


# left-closed variant (different code path inside reclass_vector_range)
ref_cl_lc <- run_with_parallel(FALSE, 0L,
	quote(classify(r, rcl_ftb, include.lowest=TRUE, right=FALSE)))
par_cl_lc <- run_with_parallel(TRUE, 4L,
	quote(classify(r, rcl_ftb, include.lowest=TRUE, right=FALSE)))
expect_equal(values(par_cl_lc), values(ref_cl_lc))


# ---- multi-output subst (replaceValues 'mout' path) ----

ref_mout <- run_with_parallel(FALSE, 0L,
	quote(subst(r, from=1:3, to=cbind(c(101,102,103), c(201,202,203)))))
par_mout <- run_with_parallel(TRUE, 4L,
	quote(subst(r, from=1:3, to=cbind(c(101,102,103), c(201,202,203)))))
expect_equal(values(par_mout), values(ref_mout))


# ---- multi-input subst (replaceValues 'min' path) ----

r2 <- c(r, r * 2L)
from_mat <- matrix(c(1,2, 2,4, 3,6), ncol=2, byrow=TRUE)
ref_min <- run_with_parallel(FALSE, 0L,
	quote(subst(r2, from=from_mat, to=c(901, 902, 903))))
par_min <- run_with_parallel(TRUE, 4L,
	quote(subst(r2, from=from_mat, to=c(901, 902, 903))))
expect_equal(values(par_min), values(ref_min))
