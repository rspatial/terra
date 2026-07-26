library(terra) 


cc <- system.time({

elev <- rast(system.file("ex/elev.tif", package = "terra"))
elev <- project(elev, "epsg:32632")

lambda <- 1 ## try also 0 (default)
flowdir <- flowDir(elev, lambda = lambda)
pits <- pitfinder(flowdir, pits_on_boundary = FALSE)
###
flowdir |> project(y="epsg:3586") |> plet(tile="OpenTopoMap") 
###

## Convergence after 10000 itarations 
elev2 <- pitfiller(x = elev, pit = pits,lambda=lambda,niter=10000) ##,niter=500,D=3000)

##flowdir2 <- terrain(elev2, "flowdir")
flowdir2 <- flowDir(elev2, lambda = lambda)
pits2 <- pitfinder(flowdir2, pits_on_boundary = FALSE)

})
((elev2-elev)+(pits2>0)) |> project(y="epsg:3586",method="near") |> plet(tile="OpenTopoMap") 

mm <- pits ##((elev2-elev))
mm[pits2==0] <- NA 
mm |> project(y="epsg:3586",method="near") |> plet(tile="OpenTopoMap") 

