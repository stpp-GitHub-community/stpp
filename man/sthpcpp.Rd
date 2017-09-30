\name{sthpcpp}
\alias{sthpcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatio-temporal hot-spots cluster point process model
%%  ~~function to do ... ~~
}
\description{Generate a random spatio-temporal point pattern, a simulated realisation of the cluster process.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sthpcpp(lambp, r, mu, s.region, t.region)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambp}{Intensity of the Poisson process of cluster centres.
%%     ~~Describe \code{lambp} here~~
}
  \item{r}{Radius parameter of the clusters.
%%     ~~Describe \code{r} here~~
}
  \item{mu}{Mean number of points per cluster (a single positive number) or reference intensity for the cluster points (a function or a pixel image).
%%     ~~Describe \code{mu} here~~
}
  \item{s.region}{A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.
%%     ~~Describe \code{s.region} here~~
}
  \item{t.region}{A vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.
%%     ~~Describe \code{t.region} here~~
}
}
\details{This function generates a realisation of spatio-temporal cluster process, which can be considered as generalisation of the classical Matern cluster process, inside the spatio-temporal window.

Consider a Poisson point process in the plane with intensity \eqn{\lambda_{p}} as cluster centres for all times 'parent', as well as a infinite cylinder of radius \eqn{R} around of each Poisson point, orthogonal to the plane. The scatter uniformly in all cylinders of all points which are of the form \eqn{(x,y,z)}, the number of points in each cluster being random with a Poisson (\eqn{\mu}) distribution. The resulting point pattern is a spatio-temporal cluster point process with \eqn{t=z}. This point process has intensity \eqn{\lambda_{p} * \mu}.

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
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Not run:
#################

X <- sthpcpp(lambp=20, r=0.05, mu=100)
plot(X$xyt)

#par(mfrow=c(1,1))
#scatter3D(X$xyt[,1],X$xyt[,2],X$xyt[,3],theta=45,phi=30,
#main="Spatio-temporal point pattern",xlab="\n x",ylab="\n y",
#zlab="\n t",ticktype="detailed",col="black")

## Spatio-temporal hot-spots cluster point process model
data(northcumbria)
Northcumbria <- northcumbria/1000
Xo <- sthpcpp(lambp=0.0035, r=5, mu=200, s.region=Northcumbria, 
                   t.region=c(28,198))
plot(Xo$xyt,s.region=Northcumbria)

# plot scatterplot3d
#par(mfrow=c(1,1))
#scatter3D(Xo$xyt[,1],Xo$xyt[,2],Xo$xyt[,3],theta=45,phi=30,
#main="Spatio-temporal point pattern",xlab="\n x",ylab="\n y",
#zlab="\n t",ticktype="detailed",col="black")

## End(Not run)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")%
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line