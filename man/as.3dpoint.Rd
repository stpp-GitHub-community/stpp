\name{as.3dpoints}
\alias{as.3dpoints}
\title{Create data in spatio-temporal point format}

\description{
Create data in spatio-temporal point format.}

\usage{
as.3dpoints(...)
}

\arguments{
  \item{...}{any object(s), such as x, y and t vectors of the same length, 
  or a list or data frame containing x, y and t vectors. Valid options for 
  ... are: a points object: returns it unaltered; a list with x, y and t 
  elements of the same length: returns a points object with the x, y and t 
  elements as the coordinates of the points; three vectors of equal length:
  returns a points object with the first vector as the x coordinates, the 
  second vector as the y-coordinates and the third vector as the t-coordinates. }
 }

\value{
The output is an object of the class \code{stpp}.
as.3dpoints tries to return the argument(s) as a spatio-temporal points object.
}

\author{Edith Gabriel <edith.gabriel@univ-avignon.fr>, Peter Diggle, Barry Rowlingson.}
