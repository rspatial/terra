
x <- rast(nrows=2, ncols=2, vals=1:4, names="random")
time(x) <- Sys.time()
saveRDS(x, "test.rds")
y <- readRDS("test.rds")
expect_true(all.equal(x, y))

