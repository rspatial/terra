selectSize <- function(x,
                       n = Inf,
                       largest = TRUE,
                       threshold = NULL,
                       dissolve = FALSE,
                       singleparts = TRUE
) {

  stopifnot(inherits(x, "SpatVector"),
            geomtype(x) %in% c("polygons", "lines"))

  if (dissolve) {
    x <- aggregate(x)
  }

  if (singleparts) {
    x <- disagg(x)
  }

  if (is.polygons(x)) {
    size <- expanse(x)
  } else if (is.lines(x)) {
    size <- perim(x)
  }

  if (!is.null(threshold)) {
    if (largest) {
      keep <- unname(size >= threshold)
    } else {
      keep <- unname(size <= threshold)
    }

    x <- x[keep, ]
    size <- size[keep]
  }

  size_rank <- order(size)

  if (is.finite(n)) {
    if (largest) {
      x <- x[rev(size_rank)[1:n], ]
    } else {
      x <- x[size_rank[1:n], ]
    }
  }

  x
}
