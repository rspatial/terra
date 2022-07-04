\name{factors}

\docType{methods}

\alias{is.factor}
\alias{is.factor,SpatRaster-method}
\alias{as.factor}
\alias{as.factor,SpatRaster-method}
\alias{levels}
\alias{levels,SpatRaster-method}
\alias{levels<-}
\alias{levels<-,SpatRaster-method}
\alias{cats}
\alias{cats,SpatRaster-method}
\alias{categories}
\alias{categories,SpatRaster-method}
  
  
\title{Categorical rasters}

\description{

A SpatRaster layer can be a categorical variable (factor). Like \code{\link{factor}}s, categories are stored as indices (integer values) that have an associated label. 

The categories can be inspected with \code{levels} and \code{cats}. They are represented by a \code{data.frame} that must have two or more columns, the first one identifying the cell values and the other column(s) providing the category labels. 

If there are multiple columns with categories, you can set the "active" category. 

\code{cats} returns the entire data.frame, whereas \code{levels} only return two columns: the index and the active category.

To set categories for the first layers you can provide \code{levels<-} with a data.frame or a list with a data.frame. To set categories for multiple layers you can provide \code{levels<-} with a list with one element for each layer. Use \code{categories} to set the categories for a specific layer or specific layers.
}

\usage{
\S4method{levels}{SpatRaster}(x)

\S4method{levels}{SpatRaster}(x)<-value

\S4method{cats}{SpatRaster}(x, layer, active=FALSE)

\S4method{categories}{SpatRaster}(x, layer=1, value, index)

\S4method{is.factor}{SpatRaster}(x)

\S4method{as.factor}{SpatRaster}(x)
}


\arguments{
  \item{x}{SpatRaster}
  \item{layer}{positive integer, the layer number or name}
  \item{active}{logical. If \code{TRUE}, only return the active category}
  \item{value}{a data.frame (ID, category) or vector with category names}
  \item{index}{positive integer, indicating the column in \code{data.frame} \code{value} to be used as the (active) category, (not counting the first column with the cell values)}
}  


\value{
list of data.frames (levels, cats) or logical (is.factor)
}

\seealso{\code{\link{activeCat}}, \code{\link{catalyze}}, \code{\link{set.cats}}}

\examples{
set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE)
is.factor(r)

cls <- data.frame(id=1:3, cover=c("forest", "water", "urban"))
levels(r) <- cls
is.factor(r)
r

plot(r, col=c("green", "blue", "light gray"))
text(r, digits=3, cex=.75, halo=TRUE)

# raster starts at 3
x <- r + 2
is.factor(x)

# Multiple categories
d <- data.frame(id=3:5, cover=cls[,2], letters=letters[1:3], value=10:12)
levels(x) <- d
x

# get current index
activeCat(x)
# set index 
activeCat(x) <- 3
activeCat(x)
activeCat(x) <- "letters"
plot(x, col=c("green", "blue", "light gray"))
text(x, digits=3, cex=.75, halo=TRUE)

r <- as.numeric(x)
r

p <- as.polygons(x)
plot(p, "letters", col=c("green", "blue", "light gray"))
}

\keyword{methods}
\keyword{spatial}

