
r <- rast(nrow=2, ncol=2, vals=1:4)
x <- c(r, r)
x[1] <- cbind(5, 6)
x[2] <- c(7, 8)
expect_equal(as.vector(values(x)), c(5,7,3,4,6,8,3,4))


r <- rast(nrow=20, ncol=20)
x <- setValues(r, 1:400)
y <- setValues(r, rep(1:5, 80))
set.seed(1)
x[sample(ncell(r), 100)] <- NA
y <- mask(y, x)
z1 <- rast(x)
z1[y==1] <- x[y==1]

y <- y == 1
z1[y] <- x[y]

z2 <- mask(x, y, maskvalue=FALSE)
z3 <- ifel(y, x, NA)
z <- as.data.frame(c(z1, z2, z3))
expect_true(all(z[,1] == z[,2]) && all(z[,1] == z[,3]))
 
