
## Creation of a Digital Elevation Model 

elev <- array(NA,c(9,9))
dx <- 1
dy <- 1 
for (r in 1:nrow(elev)) {
  x <- (r-5)*dx
  for (c in 1:ncol(elev)) {
    
    y <- (c-5)*dy
    elev[r,c] <- 10+5*(x^2+y^2)
  }
} 

elev <- cbind(elev,elev) 
elev <- rbind(elev,elev) 
elev <- rast(elev)





## Flow Directions

flowdir<- terrain(elev,v="flowdir")


#t(array(flowdir[],rev(dim(flowdir)[1:2])))

## Pit Detect

pits1 <- pitfinder(flowdir)

xypit <- as.data.frame(pits1,xy=TRUE)
names(xypit) <- c("x","y","pit")
xypit$icell <- 1:nrow(xypit)
xypit[which(xypit$pit!=0),]


xypit2 <-  xypit[which(xypit$pit!=0),]
# > xypit2
# x    y pit icell
# 77   4.5 13.5   1    77
# 86  13.5 13.5   2    86
# 239  4.5  4.5   3   239
# 248 13.5  4.5   4   248
pits0 <- elev*0
pits0[c(77,86,239,248)] <- c(1,2,3,4)

result <- (pits1==pits0) 

expect_equal(all(result[]),TRUE)