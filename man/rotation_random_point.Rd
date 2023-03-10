% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rotation.R
\name{rotation_random_point}
\alias{rotation_random_point}
\title{Generate a random point on a rotation manifold (SO(n))}
\usage{
rotation_random_point(n)
}
\arguments{
\item{n}{the dimension of the rotation manifold.}
}
\value{
an orthogonal matrix with \code{Matrix::det(x) == 1}
}
\description{
The elements of a rotation manifold (also called special orthogonal) are
\eqn{\text{SO}(n) := \{x \in \mathbb{R}^{n\times n} | x^Tx = xx^T = I_n, \text{det}(x) = 1\}}
}
