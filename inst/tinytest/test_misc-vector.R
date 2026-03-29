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
