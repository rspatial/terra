
r <- rast(nrows=10, ncols=10)
set.seed(2)
v <- sample(5, ncell(r), replace=TRUE)
values(r) <- v

s <- r
values(s) <- sample(3, ncell(s), replace=TRUE)
x_cross <- c(r, s)

ct <- crosstab(x_cross)
expect_true(inherits(ct, "table"))
tbl_cross <- table(values(r), values(s))
expect_equivalent(as.vector(ct), as.vector(tbl_cross))

# crosstab long=TRUE
ct_long <- crosstab(x_cross, long=TRUE)
expect_true(is.data.frame(ct_long))
expect_equal(ncol(ct_long), 3) # layer1, layer2, Freq/count
# Order might differ, but sum of counts should match
expect_equal(sum(ct_long[,3]), ncell(r))

# crosstab useNA
r_na <- r
r_na[1] <- NA
s_na <- s
s_na[2] <- NA
x_na <- c(r_na, s_na)

ct_na <- crosstab(x_na, useNA=TRUE)
expect_true(any(is.na(rownames(ct_na)) | rownames(ct_na) == "NA")) # depending on representation
expect_true(sum(ct_na) == ncell(r)) # Should sum to total cells if NAs included

ct_no_na <- crosstab(x_na, useNA=FALSE)
expect_true(sum(ct_no_na) < ncell(r))
