
x <- rast(nrows=10, ncols=10, xmin=0, xmax=10, vals = 1)
y <- rast(nrows=10, ncols=10, xmin=0, xmax=10, vals = 2)
z <- rast(nrows=10, ncols=10, xmin=0, xmax=10, vals = 3)
wt <- c(x,y,z)
z[1,1] <- NA
xt <- c(x,y,z)

wm <- terra::weighted.mean(xt, wt, na.rm = TRUE) 
terra_wm <- as.numeric(global(wm, fun = "min"))
# |> not available on older R
#  global(fun = "min") |>
#  as.numeric()
stats_wm <- stats::weighted.mean(x = c(1,2,NA), w = 1:3, na.rm = TRUE)
expect_equal(terra_wm, stats_wm)

