\name{rlgcp}
\alias{rlgcp}
\title{Generate log-Gaussian Cox point patterns}

\description{
Generate one (or several) realisation(s) of the log-Gaussian cox process 
in a region \eqn{S\times T}{S x T}.
}

\usage{
 rlgcp(s.region, t.region, replace=TRUE, npoints=NULL, nsim=1, 
 nx=100, ny=100, nt=100,separable=TRUE,model="exponential",
 param=c(1,1,1,1,1,1), scale=c(1,1),var.grf=1,mean.grf=0,
 lmax=NULL,discrete.time=FALSE,exact=FALSE,anisotropy=FALSE,ani.pars=NULL)
}

\arguments{
  \item{s.region}{Two-column matrix specifying polygonal region containing
  all data locations. If \code{s.region} is missing, the unit square is considered.} 
  \item{t.region}{Vector containing the minimum and maximum values of
  the time interval. If \code{t.region} is missing, the interval \eqn{[0,1]}{[0,1]} is considered.} 
  \item{npoints}{Number of points to simulate. If \code{NULL}, the
  number of points is from a Poisson distribution with mean the double integral of 	  
  \eqn{\lambda(s,t)}{\lambda(s,t)} over \code{s.region} and \code{t.region}. }
  \item{nsim}{number of simulations to generate. Default is 1.}
  \item{separable}{Logical. If \code{TRUE}, the covariance function of the
  Gaussian random field is separable. }
  \item{model}{Vector of length 1 or 2 specifying the model(s) of
  covariance of the Gaussian random field. If \code{separable=TRUE} and 
  \code{model} is of length 2, then the elements of \code{model} define the 
  spatial and temporal covariances respectively.
   If \code{separable=TRUE} and \code{model} is of length 1, then the
  spatial and temporal covariances belongs to the same class of covariances, 
  among "matern", "exponential", "stable", "cauchy" and "wave" (see Details). 
   If \code{separable=FALSE}, \code{model} must be of length 1 and is either
  "gneiting" or "cesare" (see Details).} 
  \item{param}{\eqn{(\alpha_1,\alpha_2,\alpha_3,\alpha_4,\alpha_5,\alpha_6)}{(\alpha_1, \alpha_2, \alpha_3, \alpha_4, \alpha_5, \alpha_6)}. Vector of parameters of the covariance function (see
  Details).} 
  \item{scale}{Vector of length 2 defining the spatial and temporal scale.}
  \item{var.grf}{Variance of the Gaussian random field.}
  \item{mean.grf}{Mean of the Gaussian random field.}
  \item{replace}{Logical allowing times repeat.} 
  \item{nx,ny,nt}{Define the size of the 3-D grid on which the intensity
  is evaluated.} 
  \item{lmax}{Upper bound for the value of \eqn{\lambda(x,y,t)}{\lambda(x,y,t)}.} 
  \item{discrete.time}{If \code{TRUE}, times belong to \eqn{N}{N},
  otherwise belong to \eqn{R^+}{R^+}.} 
  \item{exact}{logical allowing exact simulation of Gaussian random fields (see manual for details).}
  \item{anisotropy}{If \code{TRUE}, simulate an anisotropic point pattern. Currently only implemented for 
  separable covariance functions.}
  \item{ani.pars}{Vector of length 2, the anisotropy angle 
  and the anisotropy ratio, respectively (see details).}
}

\details{
We implemented stationary, isotropic spatio-temporal covariance functions.

\emph{Separable covariance functions}

\deqn{c(h,t) = c_s(\| h \|) \, c_t(|t|) , h \in S, t \in T}{c(h,t)
=c_s(||h|) c_t(|t|), h in S, t in T}

The purely spatial and purely temporal covariance functions can be:
\itemize{
\item Exponential: \eqn{c(r) = \exp(-r)}{c(r)=exp(-r)},
\item Stable: \eqn{c(r) = \exp(-r^\alpha)}{c(r)=exp(-r^\alpha)}, 
\eqn{\alpha \in [0,2]}{\alpha in [0,2]},
\item Cauchy: \eqn{c(r) = (1+r^2)^{-\alpha}}{c(r)=(1+r^2)^(-\alpha)},
\eqn{\alpha >0}{\alpha > 0},
\item Wave: \eqn{c(r) = \frac{\sin(r)}{r}}{c(r)=sin(r)/r} if \eqn{r>0}{r>0}, 
 \eqn{c(0)=1}{c(0)=1},
\item Matern: \eqn{c(r) = \frac{(\alpha r)^\nu}{2^{\nu-1}\Gamma(\nu)}
  {\cal K}_{\nu}(\alpha r)}{c(r)={(\alpha r)^\nu}/{2^(\nu-1)\Gamma(\nu)}
  K_\nu(\alpha r)}, \eqn{\nu > 0}{\nu>0} and \eqn{\alpha > 0}{\alpha>0}.

\eqn{{\cal K}_{\nu}}{K_\nu} is the modified Bessel function of second kind:
\deqn{{\cal K}_{\nu}(x) = \frac{\pi}{2} \frac{I_{-\nu}(x) -
    I_{\nu}(x)}{\sin(\pi \nu)},}{K_\nu(x) = (\pi/2) * (I_(-\alpha)(x) 
    - I_\nu(x))/sin(\pi \nu),} with
  \eqn{I_{\nu}(x) = \left( \frac{x}{2} \right)^{\nu} \sum_{k=0}^\infty
  \frac{1}{k! \Gamma(\nu+k+1)} \left( \frac{x}{2} \right)^{2k}}{I_\nu(x) 
  = (x/2)^\nu sum_{k=0}^{\infty} 1/(k! \Gamma(\nu+k+1)) (x/2)^(2k)}.
  
The parameters \eqn{\alpha_1}{\alpha_1} and \eqn{\alpha_2}{\alpha_2} correspond 
to the parameters of the spatial and temporal covariance respectively. For
the Matern model, the parameters \eqn{\alpha_1}{\alpha_1}, \eqn{\alpha_3}{\alpha_3} 
and \eqn{\alpha_2}{\alpha_2}, \eqn{\alpha_4}{\alpha_4} correspond 
to the parameters \eqn{\nu}{\nu}, \eqn{\alpha}{\alpha} of the spatial and temporal 
covariance.
}


\emph{Non-separable covariance functions}

The spatio-temporal covariance function can be:
\itemize{
\item Gneiting: \eqn{c(h,t) = \psi (t^2/\beta_2)^{-\alpha_6} \phi \left(\frac{h^2/
      \beta_1}{\psi (t^2/\beta_2)} \right)}{c(h,t)=\psi(t^2/\beta_2)^(-\alpha_6) 
      \phi( (h^2/\beta_1)/\psi(t^2/\beta_2) )}, \eqn{\beta_1, \beta_2 >0}{\beta_1, 
      \beta_2 >0} are spatial and temporal scales respectively,
      \itemize{
  \item If \eqn{\alpha_2=1}{\alpha_2=1}, \eqn{\phi(r)}{\phi(r)} is the Stable model.
  \item if \eqn{\alpha_2=2}{\alpha_2=2}, \eqn{\phi(r)}{\phi(r)} is the Cauchy model.
  \item If \eqn{\alpha_2=3}{\alpha_2=3}, \eqn{\phi(r)}{\phi(r)} is the Matern model.
  \item If \eqn{\alpha_5=1}{\alpha_5=1}, \eqn{\psi(r) = (r^{\alpha_3}+
  1)^{\alpha_4}}{\psi(r) = (r^\alpha_3+1)^\alpha_4},
  \item If \eqn{\alpha_5=2}{\alpha_5=2}, \eqn{\psi(r) = (\alpha_4^{-1} r^{\alpha_3}
    +1)/(r^{\alpha_3}+1)}{\psi(r) = (\alpha_4^(-1) r^\alpha_3 +1)/(r^\alpha_3+1)},
  \item If \eqn{\alpha_5=3}{\alpha_5=3}, \eqn{\psi(r) = \log(r^{\alpha_3} +
   \alpha_4)/\log {\alpha_4}}{\psi(r) = log(r^\alpha_3+\alpha_4)/
   log(\alpha_4)},
                }
        
  The parameter \eqn{\alpha_1}{\alpha_1} is the respective parameter for the model of
  \eqn{\phi(\cdot)}{\phi(.)}, \eqn{\alpha_3 \in (0,1]}{\alpha_3 in (0,1]}, 
  \eqn{\alpha_4 \in (0,1]}{\alpha_4 in (0,1]} and \eqn{\alpha_6 \geq 2}{\alpha_6 >= 2}.
\item cesare: \eqn{c(h,t) = \left( 1 + (h/\beta_1)^{\alpha_1} +
    (t/\beta_2)^{\alpha_2} \right)^{-\alpha_3}}{c(h,t) = (1 + (h/\beta_1)^\alpha_1 +
    (t/\beta_2)^\alpha_2)^(-\alpha_3)}, \eqn{\beta_1, \beta_2 >0}{\beta_1, 
    \beta_2 >0}, \eqn{\alpha_1, \alpha_2 \in [1,2]}{\alpha_1, \alpha_2 in [1,2]}
     and \eqn{\alpha_3 \geq 3/2}{\alpha_3 >= 3/2}.
	}

We also implemented anisotropic Log-Gaussian Cox processes.
We considered geometric spatial anisotropy (see Moller and Toftaker, 2014). In this case the covariance function
is elliptical and anisotropy is characterized by two parameters:
the anisotropy angle \eqn{0 \leq \theta < \pi}{\pi > \theta >=0} and the anisotropy ratio \eqn{0 < \delta \leq 1}{1 >= \delta >0} of the minor axis \eqn{2 \omega \delta}{2 \omega \delta}
and the major axis \eqn{2 \omega}{2 \omega}.

\deqn{C(h,t)=C_0\left( \sqrt{h \Sigma^{-1} h'},t \right), \ h \in R^2.}{C(h,t)=C_0(sqrt(h * solve(\Sigma) * h',t)).} 

}

\value{
A list containing:
\item{xyt}{Matrix (or list of matrices if \code{nsim}>1)
containing the points \eqn{(x,y,t)}{(x,y,t)} of the simulated point pattern.
\code{xyt} (or any element of the list if \code{nsim}>1) is an object 
of the class \code{stpp}.}
\item{s.region, t.region}{parameters passed in argument.}
\item{Lambda}{\eqn{nx \times ny \times nt}{nx x ny x nt} array (or list of array if \code{nsim}>1) 
of the intensity.}
}

\references{
Chan, G. and Wood A. (1997).  An algorithm for simulating stationary
Gaussian random fields. Applied Statistics, Algorithm Section, 46, 171--181.

Chan, G. and Wood A. (1999).  Simulation of stationary Gaussian vector fields.
Statistics and Computing, 9, 265--268.

Gneiting T. (2002). Nonseparable, stationary covariance functions
for space-time data. Journal of the American Statistical Association,
97, 590--600.

Moller J. and Toftaker H. (2014). Geometric anisotropic spatial point pattern analysis
and Cox processes. Scandinavian Journal of Statistics, 41, 414--435.
}

\author{
Edith Gabriel <edith.gabriel@inrae.fr>, Peter J Diggle.
}

\seealso{
 \code{\link{plot.stpp}}, \code{\link{animation}} and \code{\link{stan}} for plotting space-time point patterns.
 }

\examples{
\donttest{
# non separable covariance function: 
lgcp1 <- rlgcp(npoints=200, nx=50, ny=50, nt=50, separable=FALSE, 
model="gneiting", param=c(1,1,1,1,1,2), var.grf=1, mean.grf=0)
N <- lgcp1$Lambda[,,1];for(j in 2:(dim(lgcp1$Lambda)[3])){N <-
N+lgcp1$Lambda[,,j]}
image(N,col=grey((1000:1)/1000));box()
animation(lgcp1$xyt, cex=0.8, runtime=10, add=TRUE, prevalent="orange")

# separable covariance function: 
lgcp2 <- rlgcp(npoints=200, nx=50, ny=50, nt=50, separable=TRUE, 
model="exponential", param=c(1,1,1,1,1,2), var.grf=2, mean.grf=-0.5*2)
N <- lgcp2$Lambda[,,1];for(j in 2:(dim(lgcp2$Lambda)[3])){N <-
N+lgcp2$Lambda[,,j]}
image(N,col=grey((1000:1)/1000));box()
animation(lgcp2$xyt, cex=0.8, pch=20, runtime=10, add=TRUE,
prevalent="orange")

# anisotropic

sigma2=0.5
simlgcp <- rlgcp(npoints=500,nx=250, ny=200, nt=50,separable=TRUE,
s.region=matrix(c(0,2,2,0,0,0,0.5,0.5),byrow=FALSE,ncol=2), model="exponential", 
param=c(1,1,1,1,1,2), var.grf=sigma2, mean.grf=-0.5*sigma2,anisotropy = TRUE,
ani.pars = c(pi/4,0.1))

N <- simlgcp$Lambda[,,1];for(j in 2:dim(simlgcp$Lambda)[3]){N <- N+simlgcp$Lambda[,,j]}
image(x=simlgcp$grid[[1]]$x,y=simlgcp$grid[[1]]$y,z=N,col=grey((1000:1)/1000));box()
points(simlgcp$xyt[,1:2],pch=19,cex=0.25,col=2)
}
\dontshow{
lgcp4 <- rlgcp(npoints=50, nx=50, ny=50, nt=50, separable=TRUE, 
model="exponential", param=c(1,1,1,1,1,2), var.grf=2, mean.grf=-0.5*2)
}
}
