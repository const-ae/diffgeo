% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grassmann.R
\name{grassmann_map}
\alias{grassmann_map}
\title{Exponential map on a Grassmann manifold}
\usage{
grassmann_map(v, base_point)
}
\arguments{
\item{v}{a tangent vector}

\item{base_point}{a point on a Grassmann manifold}
}
\description{
Go from the \code{base_point} in the direction \code{v}.
}
\details{
The exponential map on a Grassmann manifold is
\deqn{\exp_p(v) = p V \text{diag}(\cos(d)) V^T + U \text{diag}(\sin(d)) V^T,}
where \eqn{U \text{diag}(d) V^T = v} is the SVD of the tangent vector \eqn{v} and \eqn{p}
is the \code{base_point}.
}
