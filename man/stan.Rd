\name{stan}
\alias{stan}
\title{(3D) space-time data animation}
\description{
 Displays \eqn{(x,y,t)}{(x,y,t)} point data and enables dynamic highlighting of time slices.
}
\usage{
stan(xyt,tlim=range(xyt[,3],na.rm=TRUE),twid=diff(tlim)/20,
persist=FALSE,states,bgpoly,bgframe=TRUE,bgimage,
bgcol=gray(seq(0,1,len=12)),axes=TRUE)
}
\arguments{
  \item{xyt}{A 3-column matrix of \eqn{x}{x}, \eqn{y}{y}, \eqn{t}{t} coordinates}
  \item{tlim}{A two-element vector of upper and lower time limits}
  \item{twid}{The initial time window width}
  \item{persist}{Whether to display points before time window}
  \item{states}{How to display points - see Details}
  \item{bgpoly}{A polygon to draw on the background plane}
  \item{bgframe}{Whether to extend the bgpoly to the front plane}
  \item{bgimage}{An list with \eqn{x}{x}, \eqn{y}{y} vectors and \eqn{z}{z} matrix to display on the background plane}
  \item{bgcol}{A colour palette vector with which to draw the bgimage}
  \item{axes}{Logical value indicating whether labels should be added.}
}
\details{
This function requires the \code{rpanel} and \code{rgl} packages. It uses \code{rpanel} for the sliders to control the graphics, and \code{rgl} for its ability to do flicker-free graphics.

The sliders set the position and width of the temporal highlight window. For 'time' slider set to time \eqn{T}{T} and 'width' slider set to \eqn{S}{S}, highlighted points are those with time coordinate \eqn{t}{t} such that \eqn{T-S< t< T}{T-S< t< T}.

How points are shown is configured with the states parameter. This is a list of length 3 specifying how points before the time window, inside the time window, and after the time window are displayed. Each element is a list of parameters as would be passed to \code{material3d()} together with a radius element. Points are drawn as spheres with the corresponding
material and radius as a fraction of the spatial span of the data.

By default the third state is invisible, and the first two states are different. By calling with the default for states and \code{persist=TRUE}, then the first state is set to the same as the second state. This has the effect of showing all points at time \eqn{< T}{< T} with the same sphere type.

If the user specifies the states parameter, then persist is ignored. The user can emulate the persist behaviour by specifying a states list with identical parameters for states 1 and 2.

Note that each state element should specify all \code{material3d} parameters used in any of the state elements. This is to make sure the parameters are reset for each of the sets of points.

The background polygon must be a simple 2-column vector of \eqn{x}{x} and \eqn{y}{y} coordinates. When used with \code{bgframe=TRUE}, the polygon is also drawn on the front plane, and the convex hull points are connected front to back in order to visualise the space-time prism that the data are contained in.

A raster image can be displayed on the back plane by setting the bgimage parameter. This must be a list with \eqn{x}{x}, \eqn{y}{y} and \eqn{z}{z} components as needed by the image function. Note that \eqn{x}{x} and \eqn{y}{y} define the center of cells and so must be the same length as the dimensions of \eqn{z}{z} - the image function can accept \eqn{x}{x} and \eqn{y}{y} values that are one longer than the dimensions of \eqn{z}{z} to define the edges, but bgimage does not allow that.

}
\value{
A list of the slider parameters when the dialog is quitted.
}
\author{Barry Rowlingson <b.rowlingson@lancaster.ac.uk>, Edith Gabriel}

\keyword{hplot}

