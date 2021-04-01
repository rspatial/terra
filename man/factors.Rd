\name{factors}

\docType{methods}

\alias{is.factor}
\alias{is.factor,SpatRaster-method}
\alias{levels}
\alias{levels,SpatRaster-method}
\alias{levels<-}
\alias{levels<-,SpatRaster-method}
\alias{cats}
\alias{cats,SpatRaster-method}
\alias{setCats}
\alias{setCats,SpatRaster-method}
  
\title{Factors}

\description{

A SpatRaster layer can be a categorical variable (factor). Like \code{\link{factor}}s, categories are stored as indices (integers) that have an associated label. For a SpatRaster, the index starts at 0, and cannot exceed 255. 

The categories can be inspected with \code{levels} and \code{cats}. With \code{levels<-} you can set the categories of the first layer by providing a vector of labels (the first value will be for cells with value 0, the second for 1, etc). You can also provide a \code{data.frame} that must have two or more columns, the first one identifying the cell values and the other column(s) providing the category labels. To set categories for multiple layers you can provide \code{levels<-} with a list with one element for each layer. 

With \code{setCats} you can set it for any layer and you can also set the variable of interest if there are multiple categories. 

}

\usage{
\S4method{is.factor}{SpatRaster}(x)

\S4method{levels}{SpatRaster}(x)

\S4method{levels}{SpatRaster}(x)<-value

\S4method{cats}{SpatRaster}(x, layer)

\S4method{setCats}{SpatRaster}(x, layer=1, value, index)
}


\arguments{
  \item{x}{SpatRaster}
  \item{layer}{positive integer, the layer number or name}
  \item{value}{a data.frame (ID, category) or vector with category names}
  \item{index}{positive integer, indicating the column in \code{value} to be used as the category}
}  


\value{
list (levels, cats) or data.frame (cats for a single layer); logical (is.factor, setCats)
}

\examples{
set.seed(0)
r <- rast(nrow=10, ncol=10)
values(r) <- sample(3, ncell(r), replace=TRUE)
is.factor(r)

cls <- c("forest", "water", "urban")
# make the raster start at zero
x <- r - 1
levels(x) <- cls
names(x) <- "land cover"
is.factor(x)
x

plot(x, col=c("green", "blue", "light gray"))
text(x, digits=3, cex=.75, halo=TRUE)

# raster starts at 3
x <- r + 2
is.factor(x)

# approach 1
levels(x) <- c("", "", "", "forest", "water", "urban")

# approach 2, also showing the use of two categories
d <- data.frame(id=3:5, cover=cls, letters=letters[1:3])
levels(x) <- d
x

## switch categories
head(cats(x, 1))
# get current index
setCats(x, 1)
# set index 
setCats(x, 1, index=3)
setCats(x, 1)
plot(x, col=c("green", "blue", "light gray"))
text(x, digits=3, cex=.75, halo=TRUE)

#setCats(x, 1, index=2)
#p <- as.polygons(x)
#plot(p, "land cover", col=c("green", "blue", "light gray"))
}

\keyword{methods}
\keyword{spatial}

