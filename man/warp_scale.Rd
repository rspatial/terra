\name{warp_scale}

\alias{warp_scale}

\title{Compute warp resampling scale}

\description{
Compute the ratio between source and destination pixel sizes for use with \code{\link{project}}. These values correspond to the GDAL warp options XSCALE and YSCALE. By default, GDAL computes them independently for each processing chunk, which can produce visible block-boundary artifacts when the ratio varies across the raster (particularly in projections with significant deformation). Setting fixed values via \code{project(..., xscale=, yscale=)} eliminates these artifacts.

The function samples a grid of points in the source raster, projects them to the destination CRS, and returns summary statistics of the local resampling ratios.
}

\usage{
warp_scale(x, y, n=21)
}

\arguments{
  \item{x}{SpatRaster. The source raster}
  \item{y}{SpatRaster or character. A template raster or a CRS description for the destination}
  \item{n}{integer. Number of sample points along each axis (the total number of sample points is \code{n*n})}
}

\value{
A list with two named numeric vectors (\code{xscale} and \code{yscale}), each containing the quantiles (0\%, 25\%, 50\%, 75\%, 100\%) of the local resampling ratio across the sampled points. The median (\code{"50\%"}) is typically the best choice for a global fixed scale. Values equal one for no resampling, below one for downsampling, and above one for upsampling.
}

\seealso{\code{\link{project}}}

\examples{
\dontrun{
a <- rast(ncols=360, nrows=180, xmin=-180, xmax=180, ymin=-90, ymax=90,
          crs="+proj=longlat +datum=WGS84")
values(a) <- 1:ncell(a)

sc <- warp_scale(a, "+proj=robin")
sc
# Use the median values
b <- project(a, "+proj=robin", xscale=sc$xscale["50\%"], yscale=sc$yscale["50\%"])
}
}

\keyword{spatial}
