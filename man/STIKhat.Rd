\name{STIKhat} \alias{STIKhat}
\title{Estimation of the Space-Time Inhomogeneous K-function}

\description{Compute an estimate of the Space-Time Inhomogeneous K-function.}

\usage{STIKhat(xyt, s.region, t.region, dist, times, lambda, 
correction="isotropic", infectious=FALSE) }

\arguments{
  \item{xyt}{Coordinates and times \eqn{(x,y,t)}{(x,y,t)} of the point pattern.}
  \item{s.region}{Two-column matrix specifying polygonal region containing all data locations. If \code{s.region} is missing, the bounding box of \code{xyt[,1:2]} is considered.}
  \item{t.region}{Vector containing the minimum and maximum values of the time interval. If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}
  \item{dist}{Vector of distances \eqn{u}{u} at which \eqn{K(u,v)}{K(u,v)} is computed. If missing, the maximum of \code{dist} is given by \eqn{\min(S_x,S_y)/4}{min(S_x, S_y)/4}, where \eqn{S_x}{S_x} and \eqn{S_y}{S_y} represent the maximum width and height of the bounding box of \code{s.region}.}
  \item{times}{Vector of times \eqn{v}{v} at which \eqn{K(u,v)}{K(u,v)} is computed. If missing, the maximum of \code{times} is given by \eqn{(T_{\max} - T_{\min})/4}{(T_max - T_min)/4}, where \eqn{T_{\min}}{T_min} and \eqn{T_{\max}}{T_max} are the minimum and maximum of the time interval \eqn{T}{T}.} 
  \item{lambda}{Vector of values of the space-time intensity function evaluated at the points \eqn{(x,y,t)}{(x,y,t)} in \eqn{S\times T}{S x T}. If \code{lambda} is missing, the estimate of the space-time K-function is computed as for the homogeneous case (Diggle et al., 1995), i.e. considering \eqn{n/|S \times T|}{n/|S x T|} as an estimate of the space-time intensity.}
 \item{correction}{A character vector specifying the edge correction(s) to be applied among \code{"isotropic"}, \code{"border"}, \code{"modified.border"}, \code{"translate"} and \code{"none"} (see Details). The default is \code{"isotropic"}.}
 \item{infectious}{Logical value. If \code{TRUE}, only future events are considered and the isotropic edge correction method is used. See Details.}}

\details{Gabriel (2014) proposes the following unbiased estimator for the STIK-function, based on data giving the locations of events \eqn{x_i: i=1,\ldots,n}{x_i: i = 1,...,n} on a spatio-temporal region \eqn{S\times T}{SxT}, where \eqn{S}{S} is an arbitrary polygon and \eqn{T}{T} is a time interval:
\deqn{\widehat{K}(u,v)=\sum_{i=1}^{n}\sum_{j\neq i}\frac{1}{w_{ij}}\frac{1}{\lambda(x_i)\lambda(x_j)}\mathbf{1}_{\lbrace \|s_i - s_j\| \leq u \ ; \ |t_i - t_j| \leq v \rbrace},}{K(u,v) = sum_{i = 1,...,n} sum_{j != i} 1/w_ij 1/(lambda(x_i) lambda(x_j)) 1{u_ij <= u ; t_j - t_i <= v},}
where \eqn{\lambda(x_i)}{lambda(x_i)} is the intensity at \eqn{x_i = (s_i,t_i)}{x_i = (s_i, t_i)} and \eqn{w_{ij}}{w_ij} is an edge correction factor to deal with spatial-temporal edge effects. The edge correction methods implemented are:

\code{isotropic}: \eqn{w_{ij} = |S \times T| w_{ij}^{(t)} w_{ij}^{(s)}}{w_ij = |S x T| w_ij^(s) w_ij^(t)}, where the temporal edge correction factor \eqn{w_{ij}^{(t)} = 1}{w_ij^(t) = 1} if both ends of the interval of length \eqn{2 |t_i - t_j|}{2|t_i - t_j|} centred at \eqn{t_i}{t_i} lie within \eqn{T}{T} and \eqn{w_{ij}^{(t)}=1/2}{w_ij^(t) = 1/2}  otherwise and \eqn{w_{ij}^{(s)}}{w_ij^(s)} is the proportion of the circumference of a circle centred at the location \eqn{s_i}{s_i} with radius \eqn{\|s_i -s_j\|}{||s_i - s_j||} lying in \eqn{S}{S} (also called Ripley's edge correction factor).

\code{border}: \eqn{w_{ij}=\frac{\sum_{j=1}^{n}\mathbf{1}\lbrace d(s_j,S)>u \ ; \ d(t_j,T) >v\rbrace/\lambda(x_j)}{\mathbf{1}_{\lbrace d(s_i,S) > u \ ; \ d(t_i,T) >v \rbrace}}}{w_ij = (sum_{j = 1,...,n} 1{d(s_j, S) > u ; d(t_j, T) > v}/
lambda(x_j)) / 1{d(s_i, S) > u ; d(t_i, T) > v}}, where \eqn{d(s_i,S)}{d(s_i, S)} denotes the distance between \eqn{s_i}{s_i}  and the boundary of \eqn{S}{S} and \eqn{d(t_i,T)}{d(t_i, T)} the distance between \eqn{t_i}{t_i} and the boundary of \eqn{T}{T}.

\code{modified.border}: \eqn{w_{ij} = \frac{|S_{\ominus u}|\times|T_{\ominus v}|}{\mathbf{1}_{\lbrace d(s_i,S) > u \ ; \ d(t_i,T) >v \rbrace}}}{w_ij = |S_(-u) x T_(-v)| / 1{d(s_i, S) > u ; d(t_i, T) > v}}, where \eqn{S_{\ominus u}}{S_(-u)} and \eqn{T_{\ominus v}}{T_(-v)} are the eroded spatial and temporal region respectively, obtained by trimming off a margin of width \eqn{u}{u} and \eqn{v}{v} from the border of the original region.

\code{translate}: \eqn{w_{ij} =|S \cap S_{s_i-s_j}| \times |T \cap T_{t_i-t_j}|}{w_ij = |S intersect S_(s_i - s_j)
 x T intersect T_(t_i - t_j)|}, where \eqn{S_{s_i-s_j}}{S_(s_i - s_j)} and \eqn{T_{t_i-t_j}}{T_(t_ i - t_j)} 
 are the translated spatial and temporal regions.

\code{none}: No edge correction is performed and \eqn{w_{ij}=|S \times T|}{w_ij = |S x T|}.

If parameter \code{infectious = TRUE}, ony future events are considered and the estimator is, using an isotropic edge correction factor (Gabriel and Diggle, 2009): 
\deqn{\widehat{K}(u,v)=\frac{1}{|S\times T|}\frac{n}{n_v}\sum_{i=1}^{n_v}\sum_{j=1; j > i}^{n_v} \frac{1}{w_{ij}} \frac{1}{\lambda(x_i) \lambda(x_j)}\mathbf{1}_{\left\lbrace u_{ij} \leq u\right\rbrace}\mathbf{1}_{\left\lbrace t_j - t_i \leq v \right\rbrace}.}{K(u,v) = 1/|SxT| n/n_v sum_{i = 1,...,n_v} sum_{j = 1,...,n_v; j > i} 1/w_ij 1/(lambda(x_i) lambda(x_j)) 1{u_ij <= u} 1{t_j - t_i <= v}.}

In this equation, the points \eqn{x_i=(s_i, t_i)}{x_i = (s_i, t_i)} are ordered so that \eqn{t_i < t_{i+1}}{t_i <
t_(i+1)}, with ties due to round-off error broken by randomly unrounding if necessary. To deal with temporal edge-effects, for each \eqn{v}{v}, \eqn{n_v}{n_v} denotes the number of events for which \eqn{t_i \leq T_1 -v}{t_i <= T_1 - v}, with \eqn{T=[T_0,T_1]}{T=[T_0, T_1]}. To deal with spatial edge-effects, we use Ripley's method.

If \code{lambda} is missing in argument, \code{STIKhat} computes an estimate of the space-time (homogeneous)
K-function: \deqn{\widehat{K}(u,v)=\frac{|S\times T|}{n_v(n-1)} \sum_{i=1}^{n_v}\sum_{j=1;j>i}^{n_v}\frac{1}{w_{ij}}\mathbf{1}_{\lbrace u_{ij}\leq u \rbrace}\mathbf{1}_{\lbrace t_j - t_i \leq v \rbrace}}{K(u,v) = |SxT|/(n_v (n-1)) sum_{i = 1,...,n_v} sum_{j = 1,...,n_v; j > i} 1/w_ij 1{u_ij <= u} 1{t_j - t_i <= v}}}

\value{ A list containing:

\item{Khat}{\code{ndist} x \code{ntimes} matrix containing values of \eqn{\hat{K}_{ST}(u,v).}{K(u,v)}.}
\item{Ktheo}{\code{ndist} x \code{ntimes} matrix containing theoretical values for a Poisson process; \eqn{\pi u^2 v}{pi u^2 v} for \eqn{K}{K} and \eqn{2 \pi u^2 v}{2 pi u^2 v} for \eqn{K^*}{K^*}.}
\item{dist, times, infectious}{Parameters passed in argument.}
\item{correction}{The name(s) of the edge correction method(s) passed in argument.}}

\references{
Baddeley A., Moller J. and Waagepetersen R. (2000). Non- and semi-parametric estimation of interaction in inhomogeneous point patterns. Statistica Neerlandica, 54, 329--350.

Baddeley, A., Rubak, E., Turner, R., (2015). Spatial Point Patterns: Methodology and Applications with R. CRC Press, Boca Raton.

Diggle P. , Chedwynd A., Haggkvist R. and Morris S. (1995). Second-order analysis of space-time clustering. Statistical Methods in Medical Research, 4, 124--136.

Gabriel E., Diggle P. (2009). Second-order analysis of inhomogeneous spatio-temporal point process data. Statistica Neerlandica, 63, 43--51.

Gabriel E., Rowlingson B., Diggle P. (2013). stpp: an R package for plotting, simulating and analyzing Spatio-Temporal Point Patterns. Journal of Statistical Software, 53(2), 1--29.

Gabriel E. (2014). Estimating second-order characteristics of inhomogeneous spatio-temporal point processes: influence of edge correction methods and intensity estimates. Methodology and computing in Applied Probabillity, 16(2), 411--431.
}

\author{Edith Gabriel <edith.gabriel@inrae.fr> }

\examples{\donttest{
# First example

data(fmd)
data(northcumbria)
FMD<-as.3dpoints(fmd[,1]/1000,fmd[,2]/1000,fmd[,3])
Northcumbria=northcumbria/1000

# estimation of the temporal intensity
Mt<-density(FMD[,3],n=1000)
mut<-Mt$y[findInterval(FMD[,3],Mt$x)]*dim(FMD)[1]

# estimation of the spatial intensity
h<-mse2d(as.points(FMD[,1:2]), Northcumbria, nsmse=50, range=4)
h<-h$h[which.min(h$mse)]
Ms<-kernel2d(as.points(FMD[,1:2]), Northcumbria, h, nx=5000, ny=5000)
atx<-findInterval(x=FMD[,1],vec=Ms$x)
aty<-findInterval(x=FMD[,2],vec=Ms$y)
mhat<-NULL
for(i in 1:length(atx)) mhat<-c(mhat,Ms$z[atx[i],aty[i]])

# estimation of the STIK function
u <- seq(0,10,by=1)
v <- seq(0,15,by=1)
stik1 <- STIKhat(xyt=FMD, s.region=northcumbria/1000,t.region=c(1,200), 
lambda=mhat*mut/dim(FMD)[1], dist=u, times=v, infectious=TRUE)

# plotting the estimation
plotK(stik1)
plotK(stik1,type="persp",theta=-65,phi=35)
 }
# Second example

xyt=rpp(lambda=200)
stik2=STIKhat(xyt$xyt,dist=seq(0,0.16,by=0.02),
times=seq(0,0.16,by=0.02),correction=c("border","translate"))
plotK(stik2,type="contour",legend=TRUE,which="translate")
}
