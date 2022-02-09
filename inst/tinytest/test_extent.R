# intersection of disjoint envelopes is NULL
a <- ext(c(xmin = 0, xmax = 10, ymin = 0, ymax = 10))
b <- ext(c(xmin = 100, xmax = 101, ymin = 100, ymax = 101))
expect_null(intersect(a, b))