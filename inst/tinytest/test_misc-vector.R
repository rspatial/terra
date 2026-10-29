# Value checks for SpatVector helpers documented in man/*.Rd.

f <- system.file("ex/lux.shp", package = "terra")
v <- vect(f)

## identity / size
expect_equal(unname(geomtype(v)[1]), "polygons")
expect_equal(nseg(v[1, ]), 330)
expect_equal(perim(v[1, ]), 117100.1, tolerance = 0.1)
expect_equal(expanse(v[1, ]), 312283206, tolerance = 1)
expect_equal(sum(expanse(v)), 2564810595, tolerance = 1)

## nearby — each feature has at least one neighbour within 50 km (lon/lat metres)
nb <- nearby(v, distance = 50000)
expect_true(length(nb) >= nrow(v))

## nearest — known distances (m) from first three to next set
nr <- nearest(v[1:3, ], v[4:8, ])
expect_equal(nr$distance[1], 16142.4, tolerance = 0.1)
expect_equal(nr$distance[2], 7482.58, tolerance = 0.05)

## buffer increases area
expect_true(expanse(buffer(v[1, ], width = 500))[1] > expanse(v[1, ])[1])

## simplifyGeom reduces vertices
simp <- simplifyGeom(v[1, ], tolerance = 100)
expect_true(nseg(simp) < nseg(v[1, ]))

## validity
expect_true(all(is.valid(v)))

## set operations change feature count or geometry measure predictably
e <- erase(v[1, ], v[2, ])
expect_true(expanse(e)[1] <= expanse(v[1, ])[1])
intr <- intersect(v[1:2, ], v[2:3, ])
expect_true(nrow(intr) >= 1)
un <- union(v[1:2, ])
expect_true(expanse(un)[1] >= max(expanse(v[1, ])[1], expanse(v[2, ])[1]))

## densify increases vertices
dns <- densify(v[1, ], 500)
expect_true(nseg(dns) >= nseg(v[1, ]))

## fillHoles does not increase area
fh <- fillHoles(v[1, ])
expect_true(expanse(fh)[1] <= expanse(v[1, ])[1] + 1)

## empty extent
expect_true(isTRUE(is.empty(ext(0, 0, 1, 1))))

## relate — same geometry is equal to itself
expect_true(relate(v[1, ], v[1, ], "equals"))

## .central_loc — shrink until convex hull fits; point is inside
bean <- vect("POLYGON ((5.9 50, 6.3 50, 6.5 49.8, 6.3 49.6, 5.9 49.5, 6.3 49.7, 6.3 49.9, 5.9 50))")
expect_false(is.related(bean, centroids(bean), "contains"))
lp <- terra:::.central_loc(bean)
expect_true(is.related(bean, vect(lp), "contains"))

# hole / C-shape: hull can never fit; must still return a point inside
donut <- vect("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0), (3 3, 7 3, 7 7, 3 7, 3 3))")
ld <- terra:::.central_loc(donut)
expect_true(is.related(donut, vect(ld), "contains"))

## aggregate by two fields keeps groups with NA in a by-column
vna <- v
vna$NAME_2[1] <- NA
ana <- aggregate(vna, c("NAME_1", "NAME_2"))
expect_equal(nrow(ana), nrow(unique(data.frame(vna$NAME_1, vna$NAME_2))))
expect_true(any(is.na(ana$NAME_2)))

## hull percent keeps the features closest to the centroid
pts <- centroids(v)
h100 <- hull(pts)
h60 <- hull(pts, percent=60)
expect_true(relate(h100, h60, "contains"))
expect_true(expanse(h60) < expanse(h100))
d <- as.vector(distance(pts, centroids(aggregate(pts))))
h60b <- hull(pts[d <= stats::quantile(d, 0.6)])
expect_true(relate(h60, h60b, "equals"))
expect_error(hull(pts, percent=0))
expect_error(hull(pts, percent=101))
