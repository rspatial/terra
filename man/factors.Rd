\name{factors}

\docType{methods}

\alias{droplevels}
\alias{droplevels,SpatRaster-method}
\alias{levels}
\alias{levels,SpatRaster-method}
\alias{levels<-}
\alias{levels<-,SpatRaster-method}
\alias{cats}
\alias{cats,SpatRaster-method}
\alias{categories}
\alias{categories,SpatRaster-method}
\alias{addCats}
\alias{addCats,SpatRaster-method}
  
  
\title{Categorical rasters}

\description{

A SpatRaster layer can represent a categorical variable (factor). Like \code{\link{factor}}s, SpatRaster categories are stored as integers that have an associated label. 

The categories can be inspected with \code{levels} and \code{cats}. They are represented by a \code{data.frame} that must have two or more columns, the first one identifying the (integer) cell values and the other column(s) providing the category labels. 

If there are multiple columns with categories, you can set the "active" category to choose the one you want to use. 

\code{cats} returns the entire data.frame, whereas \code{levels} only return two columns: the index and the active category.

To set categories for the first layer of a SpatRaster, you can provide \code{levels<-} with a data.frame or a list with a data.frame. To set categories for multiple layers you can provide \code{levels<-} with a list with one element (that either has a \code{data.frame} or is \code{NULL}) for each layer. Use \code{categories} to set the categories for a specific layer or specific layers.

\code{droplevels} removes categories that are not used (declared but not present as values in the raster) if \code{levels=NULL}.

\code{addCats} adds additional categories to a layer that already is categorical. It adds new variables, not new levels of an existing categorical variable.
}

\usage{
\S4method{levels}{SpatRaster}(x)

\S4method{levels}{SpatRaster}(x)<-value

\S4method{cats}{SpatRaster}(x, layer)

\S4method{categories}{SpatRaster}(x, layer=1, value, active=1, ...)

\S4method{droplevels}{SpatRaster}(x, level=NULL, layer=1)

\S4method{addCats}{SpatRaster}(x, value, merge=FALSE, layer=1)
}


\arguments{
  \item{x}{SpatRaster}
  \item{layer}{the layer name or number (positive integer); or 0 for all layers}
  \item{value}{a data.frame (ID, category) that define the categories. Or \code{NULL} to remove them}
  \item{active}{positive integer, indicating the column in \code{value} to be used as the active category (zero based to skip the first column with the cell values; that is 1 is the second column in \code{value})}
  \item{level}{the categories to remove for the layer specified with \code{layer}}
  \item{merge}{logical. If \code{TRUE}, the categories are combined with \code{\link[base]{merge}} using the first column of \code{value} as ID. If \code{FALSE} the categories are combined with \code{cbind}}
  \item{...}{additional arguments (none)}
}  


\value{
list of data.frames (levels, cats) or logical (is.factor)
}

\seealso{\code{\link{activeCat}}, \code{\link{catalyze}}, \code{\link{set.cats}}, \code{\link{as.factor}}, \code{\link{is.factor}}}

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

