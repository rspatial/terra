###
library(terra)
#> terra 1.7.29

elev <- c(78, 72, 69, 71, 58, 49,
          74, 67, 56, 49, 46, 50,
          69, 53, 44, 37, 38, 48,
          64, 58, 55, 22, 31, 24,
          68, 61, 47, 21, 16, 19,
          74, 53, 34, 12, 11, 12) |> matrix(ncol = 6, byrow = TRUE) |> rast()


flowdir0 <- c(002,002,002,004,004,008,
             002,002,002,004,004,008,
             001,001,002,004,008,004,
             128,128,001,002,004,008,
             002,002,001,004,004,004,
             001,001,001,001,000,016) |> matrix(ncol = 6, byrow = TRUE) |> rast()

flowacc0 <- c(001,001,001,001,001,001,
              001,002,002,003,003,001,
              001,004,008,006,005,001,
              001,001,001,021,001,002,
              001,001,001,002,025,001,
              001,003,005,008,036,002) |> matrix(ncol = 6, byrow = TRUE) |> rast()



flowdir1 <- terrain(elev,"flowdir")
flowacc1 <- flowAccumulation(flowdir1)

result <- (flowacc1==flowacc0) & (flowdir1==flowdir0)

expect_equal(all(result[]),TRUE)
