
r <- rast(nrows=10, ncols=10)
set.seed(2)
v <- sample(5, ncell(r), replace=TRUE)
values(r) <- v

f <- freq(r)
expect_true(is.data.frame(f))
expect_equal(colnames(f), c("layer", "value", "count"))
tbl <- table(v)
expect_equivalent(f$count, as.vector(tbl))
expect_equivalent(f$value, as.numeric(names(tbl)))

# freq bylayer=FALSE
x <- c(r, r)
f2 <- freq(x, bylayer=FALSE)
expect_equal(colnames(f2), c("value", "count"))
# Should be double the counts
expect_equivalent(f2$count, as.vector(tbl) * 2)

# freq digits
r_float <- r + 0.1
f3 <- freq(r_float, digits=0)
expect_equivalent(f3$count, as.vector(tbl))
expect_equivalent(f3$value, as.numeric(names(tbl))) # rounded values should match integers

# freq value argument
f4 <- freq(r, value=5)
expect_equal(f4$count, as.integer(tbl["5"]))

# freq zones
z <- r
values(z) <- 1
f5 <- freq(r, zones=z)
expect_true(is.data.frame(f5))
