``` r
library(terra)
#> terra 1.7.33
x <- rast(matrix(rep(1:3, 3), ncol = 3))

# existing behavior: nominal factors -> character
levels(x) <- data.frame(ID=1:3, category=factor(letters[3:1]))
str(levels(x)[[1]])
#> 'data.frame':    3 obs. of  2 variables:
#>  $ ID      : int  1 2 3
#>  $ category: chr  "c" "b" "a"
plot(x)
```

![](https://i.imgur.com/TKlArnN.png)<!-- existing behavior is that nominal factors are converted to character -->
  
``` r
# new behavior: ordinal factors are retained as factor
levels(x) <- data.frame(ID=1:3, category=factor(c(letters[3:1]),
                                                levels = letters[1:3], 
                                                ordered = TRUE))
str(levels(x)[[1]])
#> 'data.frame':    3 obs. of  2 variables:
#>  $ ID      : int  1 2 3
#>  $ category: Factor w/ 3 levels "a","b","c": 3 2 1
plot(x)
```

![](https://i.imgur.com/NdPApST.png)<!-- the primary new feature is that ordinal factors are retained as factor in the category table, and used for plotting/legends -->
  
``` r
# new behavior: ordinal factors with unused levels are retained in legend
#               (possibly should only happen with all_levels=TRUE?)
levels(x) <- data.frame(ID=1:3, category=factor(c(letters[3:1]),
                                                levels = letters[1:5], 
                                                ordered = TRUE))
str(levels(x)[[1]])
#> 'data.frame':    3 obs. of  2 variables:
#>  $ ID      : int  1 2 3
#>  $ category: Factor w/ 5 levels "a","b","c","d",..: 3 2 1
plot(x)
```

![](https://i.imgur.com/rSWSRnp.png)<!-- also supported are ordinal factors with more levels than are used in the data/present in the category table itself -->
  
``` r
# problem: cannot roundtrip categories with order/extra levels
x <- wrap(x)
saveRDS(x, file="test.rds")
x <- unwrap(readRDS("test.rds"))
plot(x)
```
  
``` r
# problem: cannot roundtrip categories with order/extra levels
x <- writeRaster(x, "test.tif", overwrite = TRUE)
plot(x)
```

![](https://i.imgur.com/FEoiQZ8.png)<!-- unfortunately, as currently implemented (perhaps in general) category order and extra levels not used in the category table are not supported in file/GeoTIFF output; fundamental GDAL limitation (?)-->
