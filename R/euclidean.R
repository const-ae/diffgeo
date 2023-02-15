
#' Check if a point is from (the tangent space of) a Euclidean manifold
#'
#' @param x point for which we check if it is on a Euclidean manifold
#' @param v,base_point tangent vector for we which we check if is in the
#'   tangent space of a Euclidean manifold at `base_point`
#' @param error a flag to indicate if the function should
#'   throw an error if the check fails
#'
#' @return
#'  `TRUE` if the check succeeds, `FALSE` if the check fails (
#'  and `error = FALSE`).
#'
#' @export
euclidean_check_point <- function(x, error = TRUE){
  if(! is.numeric(x)){
    if(error) stop("The elements of the Euclidean manifold must be numbers (i.e., scalars, vectors or matrices)")
    FALSE
  }
  TRUE
}

#' @rdname euclidean_check_point
#' @export
euclidean_check_tangent <- function(v, base_point, error = TRUE){
  if(! is.numeric(v)){
    if(error) stop("The elements of the Euclidean manifold must be numbers (i.e., scalars, vectors or matrices)")
    FALSE
  }
  if(! all(DIM(v) == DIM(base_point))){
    if(error) stop("The dimensions of 'v' and 'base_point' must match")
    FALSE
  }
  TRUE
}

#' @rdname sphere_injectivity_radius
#' @export
euclidean_injectivity_radius <- function(){
  Inf
}


#' Generate a random point on a Euclidean manifold
#'
#' The elements of a Euclidean manifold (also called special orthogonal) are
#' \eqn{\text{Eu}(n,\cdot) := \{x \in \mathbb{R}^{n\times \cdot}\}}
#'
#' @param n the dimension of the Euclidean manifold (can be a vector)
#' @param ... additional arguments passed to [`rnorm`]
#'
#' @return an array
#'
#' @export
euclidean_random_point <- function(n, ...){
  array(rnorm(prod(n), ...), dim = n)
}

#' Generate a random tangent vector to a point on a Euclidean manifold
#'
#' The elements of the tangent space for a \eqn{p \in \text{SO}(n)}
#' are skew-symmetric matrices (\eqn{\mathcal{T}_p\text{Eu}(n,\cdot) :=  \{v \in \mathbb{R}^{n\times \cdot}\}})
#'
#' @param base_point the point from a euclidean manifold
#'
#' @return an array with the same dimensions as `base_point`
#'
#' @export
euclidean_random_tangent <- function(base_point, ...){
  array(rnorm(prod(DIM(base_point)), ...), dim = DIM(base_point))
}


#' Project a point to (the tangent space of ) a Euclidean manifold
#'
#' @param x,v,base_point arrays
#'
#' @return the projected arrays
#'
#' @export
euclidean_project_point <- function(x){
  x
}

#' @rdname rotation_project_point
#' @export
euclidean_project_tangent <- function(v, base_point){
  v
}


#' Exponential map on a Euclidean manifold
#'
#' Go from the `base_point` in the direction `v`.
#'
#' @param v a tangent vector
#' @param base_point a point on a Euclidean manifold
#'
#' @details
#'   The exponential map on a Euclidean manifold is
#'     \deqn{\exp_p(v) = p + v}
#'   where \eqn{p} is the `base_point`
#'
#' @export
euclidean_map <- function(v, base_point){
  base_point + v
}

#' Find the tangent vector that connects two points on a Euclidean manifold
#'
#' @param base_point,target_point two points on a Euclidean manifold
#'
#' @details
#'   The inverse of the exponential map on a Euclidean manifold is
#'     \deqn{\log(p, q) = q - p}
#'   where \eqn{p} is the `base_point` and \eqn{q} is the `target_point`.
#'
#' @export
euclidean_log <- function(base_point, target_point){
  target_point - base_point
}




