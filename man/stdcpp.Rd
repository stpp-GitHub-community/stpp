\name{stdcpp}
\alias{stdcpp}

\title{Spatio-temporal double cluster point process model}

\description{Generate a random spatio-temporal point pattern, a simulated realisation of the cluster process.}

\usage{stdcpp(lambp, a, b, c, mu, s.region, t.region)}

\arguments{
  \item{lambp}{Intensity of the Poisson process of cluster centres. A single positive number or function.}
  \item{a}{Length of the semi-axes x of ellipsoid.}
  \item{b}{Length of the semi-axes y of ellipsoid.}
  \item{c}{Length of the semi-axes z of ellipsoid.}
  \item{mu}{Mean number of points per cluster (a single positive number).}
  \item{s.region}{A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.}
  \item{t.region}{A vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}
}

\details{We consider the straightforward extension of the classical Matern cluster process on the \eqn{R^3} case (with ellipsoid or balls) by considering the \eqn{z}-coordiantes as times.

Consider a Poisson point process in the plane with intensity \eqn{\lambda_{p}} as cluster centres for all times 'parent', as well as a ellipsoid (or ball) where the semi-axes are of lengths \eqn{a}, \eqn{b} and \eqn{c}, around of each Poisson point under a random general rotation. The scatter uniformly in all ellipsoid (or ball)  of all points which are of the form \eqn{(x,y,z)}, the number of points in each cluster being random with a Poisson (\eqn{\mu}) distribution. The resulting point pattern is a spatio-temporal cluster point process with \eqn{t=z}. This point process has intensity \eqn{\lambda_{p} * \mu}.}

\value{The simulated spatio-temporal point pattern.}

\references{
Baddeley, A., Rubak, E., Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. CRC Press, Boca Raton.

Chiu, S. N., Stoyan, D., Kendall, W. S., and Mecke, J. (2013). Stochastic Geometry and its Applications. John Wiley & Sons.

Gabriel, E., Rowlingson, B., Diggle P J. (2013) \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.

Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Stoyan, D., Rodriguez-Cortes, F. J., Mateu, J., and Gille, W. (2017). Mark variograms for spatio-temporal point processes. Spatial Statistics. 20, 125-147.
}

\author{Edith Gabriel, Francisco J. Rodriguez Cortes <cortesf@uji.es> }

\examples{
## Not run:
#################

# Ellipsoid
Xe <- stdcpp(lambp=20, a=0.12, b=0.09, c=0.07, mu=100)
plot(Xe$xyt)

# Spatio-temporal 3D scatter plot
par(mfrow=c(1,1))
plot(Xe$xyt,scatter=TRUE,theta=45,phi=30,main="Spatio-temporal point pattern",
xlab="\n x",ylab="\n y",zlab="\n t",ticktype="detailed",col="black")

# Balls
Xb <- stdcpp(lambp=20, a=0.05, b=0.05, c=0.05, mu=100)
plot(Xb$xyt)

# Spatio-temporal 3D scatter plot
par(mfrow=c(1,1))
plot(Xb$xyt,scatter=TRUE,theta=45,phi=30,
main="Spatio-temporal point pattern",xlab="\n x",ylab="\n y",
zlab="\n t",ticktype="detailed",col="black")

# Northcumbria
data(northcumbria)
Northcumbria <- northcumbria/1000
X <- stdcpp(lambp=0.00004, a=10,b=10,c=10, mu=120, s.region=Northcumbria, t.region=c(0,200))
plot(X$xyt,s.region=Northcumbria)

# Spatio-temporal 3D scatter plot
par(mfrow=c(1,1))
plot(X$xyt,scatter=TRUE,theta=45,phi=30,
main="Spatio-temporal point pattern",xlab="\n x",ylab="\n y",
zlab="\n t",ticktype="detailed",col="black")

## End(Not run)
}