\name{add_mtext}

\alias{add_mtext}

\title{add a margin text}

\description{ 
Similar to \code{\link[graphics]{mtext}} allowing adding a text to the margins of a map. This function uses the margins around the mapped area; not the margins that R would use.
}

\usage{
add_mtext(text, side=3, line=0, adj=0.5, ...)
}

\arguments{
  \item{text}{character or expression vector specifying the text to be written}
  \item{side}{ integer indicating the margin to use (1=bottom, 2=left, 3=top, 4=right)}
  \item{line}{ numeric to move the text in or outwards}
  \item{adj}{numeric between 0 and 1 to align the text along the margin. Use 0 for left (sides 1 and 3) or bottom (sides 2 and 4), 1 for right or top, and 0.5, the default, to center it}
  \item{...}{arguments passed to \code{\link{text}}}
}


\seealso{\code{\link{add_legend}}, \code{\link{add_grid}}, \code{\link{add_box}}}

\examples{
f <- system.file("ex/elev.tif", package="terra")
r <- rast(f)

plot(r, axes=FALSE, legend=FALSE)
add_box()
for (i in 1:4) add_mtext("margin text", i, cex=i, col=i, line=2-i)

# align the text along the margin
plot(r, axes=FALSE, legend=FALSE)
add_box()
for (i in 1:4) {
	add_mtext("start", i, adj=0, col="blue")
	add_mtext("middle", i, adj=0.5, cex=1.5)
	add_mtext("end", i, adj=1, col="red")
}

}

\keyword{methods}
\keyword{spatial}
