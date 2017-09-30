\name{stdcpp}
\alias{stdcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatio-temporal double cluster point process model
%%  ~~function to do ... ~~
}
\description{Generate a random spatio-temporal point pattern, a simulated realisation of the cluster process.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stdcpp(lambp, a, b, c, mu, s.region, t.region)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambp}{Intensity of the Poisson process of cluster centres. A single positive number or function.
%%     ~~Describe \code{lambp} here~~
}
  \item{a}{Length of the semi-axes x of ellipsoid.
%%     ~~Describe \code{a} here~~
}
  \item{b}{Length of the semi-axes y of ellipsoid.
%%     ~~Describe \code{b} here~~
}
  \item{c}{Length of the semi-axes z of ellipsoid.
%%     ~~Describe \code{c} here~~
}
  \item{mu}{Mean number of points per cluster (a single positive number).
%%     ~~Describe \code{mu} here~~
}
  \item{s.region}{A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.
%%     ~~Describe \code{s.region} here~~
}
  \item{t.region}{A vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.
%%     ~~Describe \code{t.region} here~~
}
}
\details{We consider the straightforward extension of the classical Matern cluster process on the \eqn{R^3} case (with ellipsoid or balls) by considering the \eqn{z}-coordiantes as times.

Consider a Poisson point process in the plane with intensity \eqn{\lambda_{p}} as cluster centres for all times 'parent', as well as a ellipsoid (or ball) where the semi-axes are of lengths \eqn{a}, \eqn{b} and \eqn{c}, around of each Poisson point under a random general rotation. The scatter uniformly in all ellipsoid (or ball)  of all points which are of the form \eqn{(x,y,z)}, the number of points in each cluster being random with a Poisson (\eqn{\mu}) distribution. The resulting point pattern is a spatio-temporal cluster point process with \eqn{t=z}. This point process has intensity \eqn{\lambda_{p} * \mu}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{The simulated spatio-temporal point pattern.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Baddeley, A., Rubak, E., Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. CRC Press, Boca Raton.

Chiu, S. N., Stoyan, D., Kendall, W. S., and Mecke, J. (2013). Stochastic Geometry and its Applications. John Wiley & Sons.

Gabriel, E., Rowlingson, B., Diggle P J. (2013) \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software 53, 1-29.

Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Stoyan, D., Rodriguez-Cortes, F. J., Mateu, J., and Gille, W. (2017). Mark variograms for spatio-temporal point processes. Spatial Statistics. 20, 125-147.
}
\author{Edith Gabriel, Francisco J. Rodriguez Cortes <cortesf@uji.es> 
}
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Not run:
#################

# Ellipsoid
Xe <- stdcpp(lambp=20, a=0.12, b=0.09, c=0.07, mu=100)
plot(Xe$xyt)

# plot scatterplot3d
#par(mfrow=c(1,1))
#scatter3D(Xe$xyt[,1],Xe$xyt[,2],Xe$xyt[,3],theta=45,phi=30,
#main="Spatio-temporal point pattern",xlab="\n x",ylab="\n y",
#zlab="\n t",ticktype="detailed",col="black")

# Balls
Xb <- stdcpp(lambp=20, a=0.05, b=0.05, c=0.05, mu=100)
plot(Xb$xyt)

# plot scatterplot3d
#par(mfrow=c(1,1))
#scatter3D(Xb$xyt[,1],Xb$xyt[,2],Xb$xyt[,3],theta=45,phi=30,
#main="Spatio-temporal point pattern",xlab="\n x",ylab="\n y",
#zlab="\n t",ticktype="detailed",col="black")

# Northcumbria
data(northcumbria)
Northcumbria <- northcumbria/1000
X <- stdcpp(lambp=0.00004, a=10,b=10,c=10, mu=120, s.region=Northcumbria, t.region=c(0,200))
plot(X$xyt,s.region=Northcumbria)

# plot scatterplot3d
#par(mfrow=c(1,1))
#scatter3D(X$xyt[,1],X$xyt[,2],X$xyt[,3],theta=45,phi=30,
#main="Spatio-temporal point pattern",xlab="\n x",ylab="\n y",
#zlab="\n t",ticktype="detailed",col="black")

## End(Not run)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")%
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line