\name{pitfiller}

\alias{pitfiller}
\alias{pitfiller,SpatRaster-method}

\title{Pit Filler}

\description{
Fills pits (depressions with no outlet) in a elevation raster using the PEM4PIT terrain erosion model.
}

\usage{
\S4method{pitfiller}{SpatRaster}(x, pit, flowdir, niter = 10, lambda = 0,
                 deviation_type = "lad", max_iters=10^6, U = 1, D = 300, beta = 0.9,
                 theta_exp = 0.5, filename = "", ...)
}

\arguments{
  \item{x}{SpatRaster with elevation values}
  \item{pit}{SpatRaster with pits identified with values > 0. See \code{\link{pitfinder}}}
  \item{flowdir}{SpatRaster with flow direction. See \code{\link{flowDir}}}
  \item{niter}{Number of iterations}
  \item{lambda}{Deviation parameter}
  \item{deviation_type}{Type of deviation. See \code{\link{flowDir}}}
  \item{max_iters}{maximum iterations for drainage path starting points detection.See \code{\link{flowDir}}}
  \item{U}{Uplift rate. See Details}
  \item{D}{Diffusion coefficient. See Details}
  \item{beta}{Erosion coefficient. See Details}
  \item{theta_exp}{Exponent for contributing area. See Details}
  \item{filename}{Character. Output filename}
  \item{...}{Additional arguments for writing files, as in \code{\link{writeRaster}}}
}

\value{
A \code{\link{SpatRaster}} with pits filled.
}

\details{
The PEM4PIT (Physical Erosion Model for PIT and flat areas correction) equation [Grimaldi at al, 2007] corrects DEMs by simulating erosion processes. It improves hydrologic modeling by enhancing flow direction and contributing area representation. Function is still experimental. 

The model is based on the steady-state sediment balance equation:

\deqn{0 = U - \beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) + \frac{4D}{\Delta x^2} \left( \bar{z} - z \right)}

Rearranging terms to isolate \eqn{z}:

\deqn{\beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) = U + \frac{4D}{\Delta x^2} \left( \bar{z} - z \right)}

Distributing terms:

\deqn{\beta A^{\theta} \frac{z}{\Delta l} - \beta A^{\theta} \frac{z_d}{\Delta l} = U + \frac{4D \bar{z}}{\Delta x^2} - \frac{4D z}{\Delta x^2}}

Combining like terms:

\deqn{\beta A^{\theta} \frac{z}{\Delta l} + \frac{4D z}{\Delta x^2} = U + \frac{4D \bar{z}}{\Delta x^2} + \beta A^{\theta} \frac{z_d}{\Delta l}}

Factoring out \eqn{z}:

\deqn{z \left( \beta A^{\theta} \frac{1}{\Delta l} + \frac{4D}{\Delta x^2} \right) = U + \frac{4D \bar{z}}{\Delta x^2} + \beta A^{\theta} \frac{z_d}{\Delta l}}

Solving for \eqn{z}:

\deqn{z = \frac{U \Delta l \Delta x^2 + 4D \Delta l \bar{z} + \beta A^{\theta} \Delta x^2 z_d}{\beta A^{\theta} \Delta x^2 + 4D \Delta l}}
}
\references{
Grimaldi, S., Nardi, F., Di Benedetto, F., Istanbulluoglu, E., & Bras, R. L. (2007). A physically-based method for removing pits in digital elevation models. \emph{Advances in Water Resources}, 30(10), 2151--2158. \doi{10.1016/j.advwatres.2006.11.016}

Tucker, G. E., & Bras, R. L. (1998). Hillslope processes, drainage density, and landscape morphology. \emph{Water Resources Research}, 34(10), 2751--2764. \doi{10.1029/98WR01474}

Niemann, J. D., Bras, R. L., & Veneziano, D. (2003). A physically based interpolation method for fluvially eroded topography. \emph{Water Resources Research}, 39, 1017. \doi{10.1029/2001WR001050}

Niemann, J. D., Bras, R. L., Veneziano, D., & Rinaldo, A. (2001). Impacts of surface elevation on the growth and scaling properties of simulated river networks. \emph{Geomorphology}, 40(1--2), 37--55. \doi{10.1016/S0169-555X(01)00036-8}
}

\author{
Emanuele Cordano
}

\seealso{
\code{\link{terrain}}, \code{\link{watershed}}, \code{\link{flowDir}}, \code{\link{pitfinder}}
}

\examples{
elev <- rast(system.file("ex/elev.tif", package = "terra"))
elev <- project(elev, "epsg:32632")

lambda <- 0.5 ## try also 0 (default)
flowdir <- flowDir(elev, lambda = lambda)
pits <- pitfinder(flowdir, pits_on_boundary = FALSE)
fill <- pitfiller(elev, pits, flowdir, lambda=lambda)

flowdir2 <- terrain(elev2, "flowdir")
flowdir2 <- flowDir(elev2, lambda = lambda)
pits2 <- pitfinder(flowdir2, pits_on_boundary = FALSE)
fill2 <- pitfiller(elev, pits2, flowdir2, lambda=lambda)
}
\keyword{spatial}
