test_that("matrix input works", {
  expect_s4_class(rast(volcano), "SpatRaster")

  ## degeneracies
  expect_s4_class(dm1 <- rast(matrix(1:10)), "SpatRaster")
  expect_s4_class(dm2 <- rast(matrix(1:10, nrow = 1)), "SpatRaster")
  expect_equal(dim(dm1), c(10L, 1L, 1L))

  expect_equal(dim(dm2), c(1L, 10L, 1L))
  expect_s4_class(rast(matrix(1)), "SpatRaster")


  expect_s4_class(rast(matrix()), "SpatRaster")

  expect_error(rast(matrix(1:4)[0, , drop = FALSE]))
#
#   Error in `values<-`(`*tmp*`, value = integer(0)) :
#     (lv%%nc) == 0 is not TRUE
#
  expect_error(rast(matrix(letters, 13)),
               "value must be numeric, integer, or logical")

})


test_that("matrix input type xyz works", {
expect_s4_class(rast(cbind(as.matrix(expand.grid(1:3, 1:4)), 1:12), type = "xyz"), "SpatRaster")
})
