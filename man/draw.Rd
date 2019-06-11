\name{draw}

\alias{drawPol}
\alias{drawLin}
\alias{drawPts}
\alias{drawExt}


\title{ Draw a line, polygon, extent, or points }

\description{
Draw on a plot (map) to get a SpatVector or SpatExtent object for later use. After calling the function, start clicking on the map. To finish, right-click and select 'stop'.
}

\usage{
drawPol(col="red", lwd=2, ...)
drawLin(col="red", lwd=2, ...)
drawPts(col="red", lwd=2, ...)
drawExt(col="red", lwd=2, ...) 

}

\arguments{
  \item{col}{the color of the lines to be drawn}
  \item{lwd}{the width of the lines to be drawn}
  \item{...}{additional arguments padded to locator}
}
 
\value{
SpatVector or SpatExtent 
}

\seealso{ \code{\link[graphics]{locator}}  }
 

\keyword{ spatial }
