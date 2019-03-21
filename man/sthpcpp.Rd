\name{sthpcpp}
\alias{sthpcpp}

\title{Spatio-temporal hot-spots cluster point process model}
\description{Generate a realisation of the hot-spots cluster process in a region \eqn{S\times T}{S x T}.}

\usage{
sthpcpp(lambp, r, mu, s.region, t.region)
}

\arguments{
\item{s.region}{Two-column matrix specifying polygonal region containing
  all data locations.If \code{s.region} is missing, the unit square is considered.}
  
\item{t.region}{Vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the interval \eqn{[0,1]}{[0,1]} is considered.}

\item{lambp}{Intensity of the Poisson process of cluster centres. A single positive number, a function, or a pixel image.}
  
\item{r}{Radius parameter of the clusters.}

\item{mu}{Average number of daughter per parent (a single positive number) or reference intensity for the cluster points (a function or a pixel image).}
}

\details{This function generates a realisation of spatio-temporal cluster process, which can be considered as generalisation of the classical Matern cluster process, inside the spatio-temporal window.

Consider a Poisson point process in the plane with intensity \eqn{\lambda_{p}}{\lambda_{p}} as cluster centres for all times 'parent', as well as a infinite cylinder of radius \eqn{R}{R} around of each Poisson point, orthogonal to the plane. The scatter uniformly in all cylinders of all points which are of the form \eqn{(x,y,z)}{(x,y,z)}, the number of points in each cluster being random with a Poisson (\eqn{\mu}{\mu}) distribution. The resulting point pattern is a spatio-temporal cluster point process with \eqn{t=z}{t=z}. This point process has intensity \eqn{\lambda_{p}\times\mu}{\lambda_{p} x \mu}.}

\value{The simulated spatio-temporal point pattern.}

\references{
Baddeley, A., Rubak, E., Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. CRC Press, Boca Raton.

Chiu, S. N., Stoyan, D., Kendall, W. S., and Mecke, J. (2013). Stochastic Geometry and its Applications. John Wiley & Sons.

Gabriel, E., Rowlingson, B., Diggle P J. (2013) \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.

Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Stoyan, D., Rodriguez-Cortes, F. J., Mateu, J., and Gille, W. (2017). Mark variograms for spatio-temporal point processes. Spatial Statistics. 20, 125-147.
}
\author{Francisco J. Rodriguez Cortes <frrodriguezc@unal.edu.co>}

\examples{
# First example
X <- sthpcpp(lambp=20, r=0.05, mu=100)
plot(X$xyt)

# Spatio-temporal 3D scatter plot
par(mfrow=c(1,1))
plot(X$xyt,type="scatter")

## Spatio-temporal hot-spots cluster point process model
data(northcumbria)
Northcumbria <- northcumbria/1000
Xo <- sthpcpp(lambp=0.0035, r=5, mu=200,
s.region=Northcumbria, t.region=c(28,198))
plot(Xo$xyt,s.region=Northcumbria)

# Spatio-temporal 3D scatter plot
par(mfrow=c(1,1))
plot(Xo$xyt,type="scatter",style="elegant",theta=45,phi=30,cex=0.05,
ticktype="detailed",col="black")
}