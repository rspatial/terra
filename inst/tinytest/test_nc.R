
f <- system.file("ex/nouragues.nc", package = "terra")
if (f == "") exit_file("nouragues.nc not available")

xy <- cbind(-52.67, 4.08)

## --- md=FALSE, single subds (t2m) -------------------------------------------
r1 <- rast(f, subds = "t2m")
expect_equal(nlyr(r1), 48)
expect_equal(dim(r1), c(3L, 3L, 48L))

v1 <- values(r1)
expect_equivalent(v1[1, 1], 296.9571, tolerance = 1e-3)
expect_equivalent(v1[9, 1], 297.0769, tolerance = 1e-3)
expect_equivalent(v1[1, 48], 296.9695, tolerance = 1e-3)

e1 <- extract(r1, xy)
expect_equivalent(as.numeric(e1[1, 1]), 297.0172, tolerance = 1e-3)


## --- md=TRUE, single subds (t2m) --------------------------------------------
r2 <- rast(f, subds = "t2m", md = TRUE)
expect_equal(nlyr(r2), 48)
expect_equal(dim(r2), c(3L, 3L, 48L))

v2 <- values(r2)
expect_equivalent(v2[1, 1], 296.9571, tolerance = 1e-3)
expect_equivalent(v2[9, 1], 297.0769, tolerance = 1e-3)

e2 <- extract(r2, xy)
expect_equivalent(as.numeric(e2[1, 1]), 297.0172, tolerance = 1e-3)

## md=TRUE and md=FALSE must give identical values and extracts
expect_equal(v1, v2)
expect_equal(e1, e2)


## --- md=FALSE, all vars (no subds) ------------------------------------------
suppressWarnings(rall <- rast(f))
expect_equal(nlyr(rall), 336)
expect_equal(dim(rall), c(3L, 3L, 336L))
expect_true(all(c("d2m", "sp", "ssrd", "t2m", "tp", "u10", "v10") %in%
                unique(varnames(rall))))

va <- values(rall)
expect_equivalent(va[1, 1], 295.4049, tolerance = 1e-3)

ea <- extract(rall, xy)
expect_equal(ncol(ea), 336)


## --- md=TRUE, all vars -------------------------------------------------------
suppressWarnings(rall2 <- rast(f, md = TRUE))
expect_equal(nlyr(rall2), 336)
expect_equal(dim(rall2), dim(rall))

va2 <- values(rall2)
expect_equal(va, va2)

ea2 <- extract(rall2, xy)
expect_equal(ea, ea2)


## --- SpatRasterDataset -------------------------------------------------------
d <- sds(f)
expect_equal(length(d), 7)
expect_true(all(c("u10", "v10", "d2m", "t2m", "sp", "tp", "ssrd") %in% names(d)))

du <- d["u10"]
expect_equal(nlyr(du), 48)
expect_equal(dim(du), c(3L, 3L, 48L))

dv <- values(du)
expect_equivalent(dv[1, 1], -1.2086, tolerance = 1e-3)
expect_equivalent(dv[9, 48], -1.0656, tolerance = 1e-3)

de <- extract(du, xy)
expect_equivalent(as.numeric(de[1, 1]), -1.1612, tolerance = 1e-3)
