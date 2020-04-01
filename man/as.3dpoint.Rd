\name{as.3dpoints}
\alias{as.3dpoints}
\title{Create data in spatio-temporal point format}

\description{
Create data in spatio-temporal point format.}

\usage{
as.3dpoints(...)
}

\arguments{
  \item{...}{Any object(s), such as \eqn{x}{x}, \eqn{y}{y} and \eqn{t}{t} vectors of the same length, or a list or data frame containing \eqn{x}{x}, \eqn{y}{y} and \eqn{t}{t} vectors. Valid options for \code{...} are: a points object: returns it unaltered; a list with \eqn{x}{x}, \eqn{y}{y} and \eqn{t}{t} elements of the same length: returns a points object with the \eqn{x}{x}, \eqn{y}{y} and \eqn{t}{t} elements as the coordinates of the points; three vectors of equal length: returns a points object with the first vector as the \eqn{x}{x} coordinates, the second vector as the \eqn{y}{y}-coordinates and the third vector as the \eqn{t}{t}-coordinates.}
}

\value{
The output is an object of the class \code{stpp}.
as.3dpoints tries to return the argument(s) as a spatio-temporal points object.
}

\author{Edith Gabriel <edith.gabriel@inrae.fr>, Peter Diggle, Barry Rowlingson.}
