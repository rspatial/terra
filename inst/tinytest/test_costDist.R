
## costDist / gridDist — target cell distance zero (?costDist)
rc <- rast(
  ncols = 5, nrows = 5, crs = "+proj=utm +zone=1 +datum=WGS84",
  xmin = 0, xmax = 5, ymin = 0, ymax = 5, vals = 1
)
rc[13] <- 0
cd <- costDist(rc, target = 0)
gd <- gridDist(rc, 0)
expect_equal(as.numeric(cd[13]), 0)
expect_equal(as.numeric(gd[13]), 0)


## --- costDist: basic distance (no nearest) -----------------------------------
r <- rast(ncols=5, nrows=5, crs="+proj=utm +zone=1 +datum=WGS84",
		xmin=0, xmax=5, ymin=0, ymax=5, vals=1)
r[13] <- 0
d <- costDist(r)
expect_equal(nlyr(d), 1)
expect_equal(values(d)[13], 0)
expect_true(all(values(d) >= 0, na.rm=TRUE))


## --- costDist: nearest=TRUE returns 2 layers ---------------------------------
dn <- costDist(r, nearest=TRUE)
expect_equal(nlyr(dn), 2)
expect_equal(names(dn), c("distance", "nearest"))

## distance layer must be identical to costDist without nearest
expect_equivalent(values(d), values(dn[[1]]))

## target cell 13 should map to itself
expect_equal(values(dn[[2]])[13], 13)

## all non-NA cells must have a nearest target ID
v_id <- values(dn[[2]])
expect_true(all(!is.na(v_id)))
expect_true(all(v_id == 13))


## --- costDist: multiple targets ----------------------------------------------
r2 <- rast(ncols=5, nrows=5, crs="+proj=utm +zone=1 +datum=WGS84",
		xmin=0, xmax=5, ymin=0, ymax=5, vals=1)
r2[3] <- 0
r2[23] <- 0
dn2 <- costDist(r2, nearest=TRUE)
expect_equal(nlyr(dn2), 2)

## target cells must map to themselves
expect_equal(values(dn2[[2]])[3], 3)
expect_equal(values(dn2[[2]])[23], 23)

## cell 1 is closer to target 3 than to target 23
expect_equal(values(dn2[[2]])[1], 3)

## cell 25 is closer to target 23 than to target 3
expect_equal(values(dn2[[2]])[25], 23)


## --- gridDist: nearest=TRUE --------------------------------------------------
g <- rast(ncols=5, nrows=5, crs="+proj=utm +zone=1 +datum=WGS84",
		xmin=0, xmax=5, ymin=0, ymax=5, vals=1)
g[13] <- 0
gd <- gridDist(g)
gn <- gridDist(g, nearest=TRUE)
expect_equal(nlyr(gn), 2)
expect_equal(names(gn), c("distance", "nearest"))
expect_equivalent(values(gd), values(gn[[1]]))
expect_equal(values(gn[[2]])[13], 13)
expect_true(all(values(gn[[2]]) == 13))


## --- gridDist: multiple targets ----------------------------------------------
g2 <- rast(ncols=5, nrows=5, crs="+proj=utm +zone=1 +datum=WGS84",
		xmin=0, xmax=5, ymin=0, ymax=5, vals=1)
g2[3] <- 0
g2[23] <- 0
gn2 <- gridDist(g2, nearest=TRUE)
expect_equal(values(gn2[[2]])[3], 3)
expect_equal(values(gn2[[2]])[23], 23)
expect_equal(values(gn2[[2]])[1], 3)
expect_equal(values(gn2[[2]])[25], 23)


## --- costDist: with NA barriers and different friction -----------------------
r3 <- rast(ncols=5, nrows=5, crs="+proj=utm +zone=1 +datum=WGS84",
		xmin=0, xmax=5, ymin=0, ymax=5, vals=2)
r3[1] <- 0
r3[25] <- 0
r3[11:15] <- NA
dn3 <- costDist(r3, nearest=TRUE)
expect_equal(nlyr(dn3), 2)
expect_equal(values(dn3[[2]])[1], 1)
expect_equal(values(dn3[[2]])[25], 25)

## cells above the NA barrier reach target 1; cells below reach target 25
expect_true(all(values(dn3[[2]])[1:10] == 1, na.rm=TRUE))
expect_true(all(values(dn3[[2]])[16:25] == 25, na.rm=TRUE))

## NA cells stay NA in distance layer
expect_true(all(is.na(values(dn3[[1]])[11:15])))
