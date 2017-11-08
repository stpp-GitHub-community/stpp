\name{LISTAhat}
\alias{LISTAhat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{Computes an edge-corrected kernel estimator of the space-time pair correlarion LISTA functions.
}
\usage{
LISTAhat(xyt, s.region, t.region, dist, times, lambda, ks = "box", hs, kt = "box", ht, correction = "isotropic")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xyt}{Spatial-temporal coordinates \eqn{(x,y,t)} of the point pattern.
}
  \item{s.region}{
%%     ~~Describe \code{s.region} here~~
}
  \item{t.region}{
%%     ~~Describe \code{t.region} here~~
}
  \item{dist}{
%%     ~~Describe \code{dist} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{ks}{
%%     ~~Describe \code{ks} here~~
}
  \item{hs}{
%%     ~~Describe \code{hs} here~~
}
  \item{kt}{
%%     ~~Describe \code{kt} here~~
}
  \item{ht}{
%%     ~~Describe \code{ht} here~~
}
  \item{correction}{
%%     ~~Describe \code{correction} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.

Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.

Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.

Siino, M., Rodriguez-Cortes, F. J., Mateu, J. and Adelfio, G. (2017). Testing for local structure in spatio-temporal point pattern data. Environmetrics. DOI: 10.1002/env.2463.

Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Ellipsoid
Xe <- stdcpp(lambp=20, a=0.12, b=0.09, c=0.07, mu=10)
plot(Xe$xyt)

# Spatio-temporal 3D scatter plot
par(mfrow=c(1,1))
plot(Xe$xyt,scatter=TRUE)

out <- LISTAhat(Xe$xyt,lambda=200,ks="gaussian",k="gaussian",correction="isotropic")

persp3D(out$dist,out$times,out$list.lista[,,190])
}
