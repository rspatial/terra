\name{xyRowColCell}

\alias{xFromCol}
\alias{xFromCol,SpatRaster,numeric-method}
\alias{xFromCol,SpatRaster,missing-method}
\alias{yFromRow}
\alias{yFromRow,SpatRaster,numeric-method}
\alias{yFromRow,SpatRaster,missing-method}
\alias{xyFromCell}
\alias{xyFromCell,SpatRaster,numeric-method}
\alias{xFromCell}
\alias{xFromCell,SpatRaster,numeric-method}
\alias{yFromCell}
\alias{yFromCell,SpatRaster,numeric-method}
\alias{cellFromRowCol}
\alias{cellFromRowCol,SpatRaster,numeric,numeric-method}
\alias{cellFromRowColCombine}
\alias{cellFromRowColCombine,SpatRaster,numeric,numeric-method}
\alias{rowColCombine}
\alias{rowColCombine,SpatRaster,numeric,numeric-method}
\alias{colFromX}
\alias{colFromX,SpatRaster,numeric-method}
\alias{rowFromY}
\alias{rowFromY,SpatRaster,numeric-method}
\alias{cellFromXY}
\alias{cellFromXY,SpatRaster,matrix-method}
\alias{cellFromXY,SpatRaster,data.frame-method}
\alias{rowFromCell}
\alias{rowFromCell,SpatRaster,numeric-method}
\alias{colFromCell}
\alias{colFromCell,SpatRaster,numeric-method}
\alias{rowColFromCell}
\alias{rowColFromCell,SpatRaster,numeric-method}


\title{Coordinates from a row, column or cell number and vice versa}

\description{
Get coordinates of the center of raster cells for a row, column, or cell number of a SpatRaster.
Or get row, column, or cell numbers from coordinates or from each other.

Cell numbers start at 1 in the upper left corner, and increase from left to right, and then from top to bottom.
The last cell number equals the number of cells of the SpatRaster.
row numbers start at 1 at the top, column numbers start at 1 at the left.
}

\usage{
\S4method{xFromCol}{SpatRaster,numeric}(object, col)

\S4method{yFromRow}{SpatRaster,numeric}(object, row)

\S4method{xyFromCell}{SpatRaster,numeric}(object, cell)

\S4method{xFromCell}{SpatRaster,numeric}(object, cell)

\S4method{yFromCell}{SpatRaster,numeric}(object, cell)

\S4method{colFromX}{SpatRaster,numeric}(object, x)

\S4method{rowFromY}{SpatRaster,numeric}(object, y)

\S4method{cellFromRowCol}{SpatRaster,numeric,numeric}(object, row, col)

\S4method{cellFromRowColCombine}{SpatRaster,numeric,numeric}(object, row, col)

\S4method{rowColCombine}{SpatRaster,numeric,numeric}(object, row, col)

\S4method{rowFromCell}{SpatRaster,numeric}(object, cell)

\S4method{colFromCell}{SpatRaster,numeric}(object, cell)

\S4method{rowColFromCell}{SpatRaster,numeric}(object, cell)

\S4method{cellFromXY}{SpatRaster,matrix}(object, xy)
}


\arguments{
  \item{object}{SpatRaster}
  \item{cell}{integer. cell number(s)}
  \item{col}{integer. column number(s) or missing (equivalent to all columns)}
  \item{row}{integer. row number(s) or missing (equivalent to all rows)}
  \item{x}{x coordinate(s)}
  \item{y}{y coordinate(s)}
  \item{xy}{matrix of x and y coordinates}
}

\details{
Cell numbers start at 1 in the upper left corner, and increase from left to right, and then from top to bottom.
The last cell number equals the number of cells of the SpatRaster (see \code{\link{ncell}}).
}

\value{
xFromCol, yFromCol, xFromCell, yFromCell: vector of x or y coordinates

xyFromCell: matrix(x,y) with coordinate pairs

colFromX, rowFromY, cellFromXY, cellFromRowCol, rowFromCell, colFromCell: vector of row, column, or cell numbers

rowColFromCell, rowColCombine: matrix of row and column numbers
}


\seealso{
\code{\link{crds}}
}



\examples{
r <- rast()

xFromCol(r, c(1, 120, 180))
yFromRow(r, 90)
xyFromCell(r, 10000)
xyFromCell(r, c(0, 1, 32581, ncell(r), ncell(r)+1))

cellFromRowCol(r, 5, 5)
cellFromRowCol(r, 1:2, 1:2)
cellFromRowCol(r, 1, 1:3)

# all combinations
cellFromRowColCombine(r, 1:2, 1:2)

colFromX(r, 10)
rowFromY(r, 10)
xy <- cbind(lon=c(10,5), lat=c(15, 88))
cellFromXY(r, xy)

# if no row/col specified all are returned
range(xFromCol(r))
length(yFromRow(r))
}

\keyword{spatial}
