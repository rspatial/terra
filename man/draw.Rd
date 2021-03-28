\name{draw}

\alias{draw}
\alias{draw,character-method}
\alias{draw,missing-method}


\title{ Draw a polygon, line, extent, or points }

\description{
Draw on a plot (map) to get a SpatVector or SpatExtent object for later use. After calling the function, start clicking on the map. To finish, press \code{ESC} (in R-Studio). In RGui you can right-click and select "stop". You can also preset the maximum number of clicks.
}

\usage{
\S4method{draw}{character}(x="extent", col="red", lwd=2, n=10000, ...)
}

\arguments{
  \item{x}{character. The type of object to draw. One of "extent", "polygon", "line", or "points"}
  \item{col}{the color to be used}
  \item{lwd}{the width of the lines to be drawn}
  \item{n}{the maximum number of clicks (does not apply when \code{x=="extent"} in which case \code{n} is always 2}
  \item{...}{additional graphics arguments for drawing}
}
 
\value{
SpatVector or SpatExtent 
}

\seealso{ \code{\link[graphics]{locator}}  }
 

\keyword{ spatial }

