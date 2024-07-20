\name{add_mtext}

\alias{add_mtext}

\title{draw a box}

\description{ 
Similar to \code{\link[graphics]{mtext}} allowing adding a text to the margins of a map. This function useds the margins around the mapped area; not the margins that R would use.
}

\usage{
add_mtext(text, side=3, line=0, ...)
}

\arguments{
  \item{text}{character or expression vector specifying the text to be written}
  \side{side} integer indicating the margin to use (1=bottom, 2=left, 3=top, 4=right).
  \side{line} numeric to move the text in or outwards.
  \item{...}{arguments passed to \code{\link{text}}}
}


\seealso{\code{\link{add_legend}}, \code{\link{add_grid}}, \code{\link{add_box}}}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

plot(r, axes=FALSE, legend=FALSE)
add_box()
for (i in 1:4) add_mtext("margin text", i, cex=i, col=i, line=2-i)
}

\keyword{methods}
\keyword{spatial}

