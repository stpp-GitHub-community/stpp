\name{stdcpp}
\alias{stdcpp}

\title{Generate double-cluster point pattern}

\description{Generate a realisation of the double-cluster process in a region \eqn{S\times T}{S x T}.}

\usage{stdcpp(lambp, a, b, c, mu, s.region, t.region)}

\arguments{
 \item{s.region}{Two-column matrix specifying polygonal region containing
  all data locations.If \code{s.region} is missing, the unit square is considered.}
  \item{t.region}{Vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the interval \eqn{[0,1]}{[0,1]} is considered.}
  \item{lambp}{Intensity of the parent process. Can be either a numeric value, a function, or a 3d-array (see \code{\link{rpp}}).}
  \item{a}{Length of the semi-axes \eqn{x}{x} of ellipsoid.}
  \item{b}{Length of the semi-axes \eqn{y}{y} of ellipsoid.}
  \item{c}{Length of the semi-axes \eqn{y}{y} of ellipsoid.}
  \item{mu}{Average number of daughter per parent. (a single positive number).}
}

\details{We consider the straightforward extension of the classical Matern cluster process on the \eqn{R^3} case (with ellipsoid or balls) by considering the \eqn{z}{z}-coordiantes as times.

Consider a Poisson point process in the plane with intensity \eqn{\lambda_p}{\lambda_p} as cluster centres for all times 'parent', as well as a ellipsoid (or ball) where the semi-axes are of lengths \eqn{a}{a}, \eqn{b}{b} and \eqn{c}{c}, around of each Poisson point under a random general rotation. The scatter uniformly in all ellipsoid (or ball)  of all points which are of the form \eqn{(x,y,z)}{(x,y,z)}, the number of points in each cluster being random with a Poisson (\eqn{\mu}{\mu}) distribution. The resulting point pattern is a spatio-temporal cluster point process with \eqn{t=z}{t=z}. This point process has intensity \eqn{\lambda_{p} \times \mu}{\lambda_{p} x \mu}.}

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
# Ellipsoid
Xe <- stdcpp(lambp=20,a=0.5,b=0.09,c=0.07,mu=100)
plot(Xe$xyt)

# Spatio-temporal 3D scatter plot
par(mfrow=c(1,1))
plot(Xe$xyt,type="scatter")

# Balls
Xb <- stdcpp(lambp=20,a=0.07,b=0.07,c=0.07,mu=100)
plot(Xb$xyt)

\donttest{
# Spatio-temporal 3D scatter plot
par(mfrow=c(1,1))
plot(Xb$xyt,type="mark",style="elegant")

# Northcumbria
data(northcumbria)
Northcumbria <- northcumbria/1000
X <- stdcpp(lambp=0.00004,a=10,b=10,c=10,mu=120,
s.region=Northcumbria,t.region=c(0,200))
plot(X$xyt,s.region=Northcumbria, cex=0.5)

# Spatio-temporal 3D scatter plot
par(mfrow=c(1,1))
plot(X$xyt,type="scatter",theta=45,phi=30,cex=0.1,
ticktype="detailed",col="black",style="elegant")
}}