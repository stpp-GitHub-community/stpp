\name{gsp}
\alias{gsp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial mark variogram function
%%  ~~function to do ... ~~
}
\description{Computes an estimator of the spatial mark variogram function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsp(xyt,s.region,s.lambda,ds,ks="epanech",hs,correction="none",approach="simplified")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyt}{Spatial coordinates and times \eqn{(x,y,t)} of the point pattern.}
  \item{s.region}{Two-column matrix specifying polygonal region containing all data locations. }
  \item{s.lambda}{Vector of values of the spatial intensity function evaluated at the points \eqn{(x,y)} in \eqn{W}. If \code{s.lambda} is missing, the estimate of the spatial mark correlation function is computed as for the homogeneous case, i.e. considering \eqn{n/|W|} as an estimate of the spatial intensity under the parameter \code{approach="standardised"}.}
  \item{ds}{A vector of distances \code{u} at which \code{gsp(u)} is computed.}
  \item{ks}{A kernel function for the spatial distances. The default is the \code{"epanech"} kernel. It can also be \code{"box"} for the uniform kernel, or \code{"biweight"}.}
  \item{hs}{A bandwidth of the kernel function \code{ks}.}
  \item{correction}{A character vector specifying the edge-correction(s) to be applied among \code{"isotropic"}, \code{"border"}, \code{"modified.border"}, \code{"translate"}, \code{"setcovf"} and \code{"none"}. The default is \code{"none"}.}
  \item{approach}{A character vector specifying the approach to use for the estimation to be applied among "simplified" or \code{"standardised"}. If approach is missing, \code{"simplified"} is considered by default.}
%%     ~~Describe \code{x} here~~
}
\details{By default, this command calculates an estimate of the spatial mark variogram function \eqn{\gamma_[sp](r)} for a spatio-temporal point pattern.}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{egsp}{ A vector containing the values of \eqn{\gamma_{sp}(r)} estimated}
\item{ds}{If \code{ds} is missing, a vector of distances \code{u} at which \code{gsp(u)} is computed from 0 to until quarter of the maximum distance between the points in the pattern.}
\item{kernel}{A vector of names and bandwidth of the spatial kernel.}
\item{gsptheo}{Value under the Poisson case is calculated considering \eqn{\tau}=\code{max(xyt[,3])-min(xyt[,3])}.}}
%% ...
\references{
%% ~put references to the literature/web site here ~
Baddeley, A., Rubak, E., Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. CRC Press, Boca Raton.

Chiu, S. N., Stoyan, D., Kendall, W. S., and Mecke, J. (2013). Stochastic Geometry and its Applications. John Wiley & Sons.

Gabriel, E., Rowlingson, B., Diggle P J. (2013) \code{stpp}: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software. 53, 1-29.

Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Stoyan, D., Rodriguez-Cortes, F. J., Mateu, J., and Gille, W. (2017). Mark variograms for spatio-temporal point processes. Spatial Statistics. 20, 125-147.
}
\author{Francisco J. Rodriguez Cortes <frrodriguezc@unal.edu.co> \url{https://fjrodriguezcortes.wordpress.com}
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Not run:
#################

# A realisation of spatio-temporal homogeneous Poisson point processes
hpp <- rpp(lambda = 100, replace = FALSE)$xyt

# R plot
plot(hpp)

# This function provides an kernel estimator of the spatial mark variogram function
out <- gsp(hpp)

# R plot - Spatial mark variogram function
par(mfrow=c(1,1))
xl <- c(0,0.25)
yl <- c(0,max(out$gsptheo,out$egsp))
plot(out$ds,out$egsp,type="l",xlab="r = distance",ylab=expression(gamma[sp](r)),
                 xlim=xl,ylim=yl,col=1,cex.lab=1.5,cex.axis=1.5)
lines(out$ds,rep(out$gsptheo,length(out$ds)),col=11)

## End(Not run)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")%
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
