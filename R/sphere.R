
#' Check if a point is from (the tangent space of) a sphere
#'
#' @param x point for which we check if it is on a sphere
#' @param v,base_point tangent vector for we which we check if is in the
#'   tangent space of a sphere at `base_point`
#' @param tol the tolerance for the test
#' @param error a flag to indicate if the function should
#'   throw an error if the check fails
#'
#' @return
#'  `TRUE` if the check succeeds, `FALSE` if the check fails (
#'  and `error = FALSE`).
#'
#' @export
sphere_check_point <- function(x, tol = 1e-12, error = TRUE){
  if(! is.matrix(x) || ncol(x) != 1){
    if(error) stop("A sphere must be specified by a 1-column matrix")
    FALSE
  }
  if(abs(1 - sum(x^2)) > tol){
    if(error) stop("sum(x^2) must equal 1.")
    FALSE
  }
  TRUE
}

#' @rdname sphere_check_point
sphere_check_tangent <- function(v, base_point, tol = 1e-12, error = TRUE){
  if(! is.matrix(v) || ncol(v) != 1){
    if(error) stop("A sphere tangent must be specified by a 1-column matrix")
    FALSE
  }
  if(abs(sum(base_point * v)) > tol){
    if(error) stop("'sum(base_point * v)' must be zero.")
    FALSE
  }
  TRUE
}

#' The injectivity radius of a manifold
#'
#' @param n the dimension of the manifold
#'
#' @details
#'   The injectivity radius is the size of the
#'   largest tangent vector for which the exponential map
#'   is a diffeomorphism, i.e., inverting the exponential map
#'   would return the original tangent vector.
#'
#'   On a sphere going from a point \eqn{p} in direction \eqn{v} beyond
#'   its pole to \eqn{q} means that the \eqn{\log(p, q) \ne v}.
#'
#' @return
#'   injectivity radii:
#'   * Sphere: \eqn{\pi}
#'   * Rotation: \eqn{\pi/\sqrt{2}}
#'   * Grassmann: \eqn{\pi/2}
#'   * Stiefel: \eqn{0.89\,\pi}
#'
#' @references
#'  Rentmeesters, Quentin. “Algorithms for Data Fitting on Some Common Homogeneous Spaces,” (2013)
#'
#' @export
sphere_injectivity_radius <- function(n){
  pi
}

#' Generate a random point on the n-sphere
#'
#' The elements of a sphere are \eqn{\text{S}(n) := \{x \in \mathbb{R}^{n+1} | \sum_i x_i^2 = 1\}}
#'
#' @param n the dimension of the sphere-manifold. The
#'   familiar sphere in 3D dimensions for `n = 2`.
#'
#' @return a matrix with one column and `n+1` rows
#'
#' @export
sphere_random_point <- function(n){
  x <- randn(n+1, 1)
  x / sqrt(sum(x^2))
}

#' Generate a random tangent vector to a point on the n-sphere
#'
#' The elements of the tangent space for a \eqn{p \in \text{S}(n)}
#' are \eqn{\mathcal{T}_p\text{S}(n) :=  \{v \in \mathbb{R}^n | \sum_i v_i p_i = 0\}}
#'
#' @param base_point the point from the n-sphere. A matrix with one column.
#'
#' @return a matrix with one column and `nrow(base_point)` rows
#'
#' @export
sphere_random_tangent <- function(base_point, ...){
  v <- randn(nrow(base_point) - 1, 1, ...)
  # Get the nullspace of the base_point
  ns <- qr.Q(qr(base_point), complete = TRUE)[,-1,drop=FALSE]
  ns %*% v
}

#' Project a point to (the tangent space of ) a sphere
#'
#' @param x,v,base_point matrix with one column
#'
#' @return the projected matrix
#'
#' @export
sphere_project_point <- function(x){
  matrix(x, ncol = 1) / sqrt(sum(x^2))
}

#' @rdname sphere_project_point
sphere_project_tangent <- function(v, base_point){
  base_point <- base_point / sqrt(sum(base_point^2))
  v - drop(t(base_point) %*% v) * base_point
}


#' Exponential map on a n-sphere
#'
#' Go from the `base_point` in the direction `v`.
#'
#' @param v a tangent vector
#' @param base_point a point on a sphere
#'
#' @details
#'   The exponential map on the sphere is
#'     \deqn{\exp_p(v) = \cos(||v||) p + \sin(||v||)\frac{v}{||v||},}
#'   where \eqn{||v||} is the norm of the tangent vector and \eqn{p}
#'   is the `base_point`.
#'
#' @export
sphere_map <- function(v, base_point){
  norm_v <- sqrt(sum(v^2))
  if(abs(norm_v < 1e-18)){
    cos(norm_v) * base_point + v
  }else{
    cos(norm_v) * base_point + sin(norm_v) / norm_v * v
  }
}

#' Find the tangent vector that connects two points on a sphere
#'
#' @param base_point,target_point two points on a sphere
#'
#' @details
#'   The inverse of the exponential map on the sphere is
#'     \deqn{\log(p, q) = \frac{\alpha}{\sin(\alpha)}(q - (p^T q) p),}
#'   where \eqn{\alpha = \cos^{-1}(p^T q)}, \eqn{p} is the `base_point`
#'   and \eqn{q} is the `target_point`.
#'
#' @export
sphere_log <- function(base_point, target_point){
  tol <- 1e-12
  sphere_check_point(base_point)
  sphere_check_point(target_point)

  # Following https://github.com/JuliaManifolds/Manifolds.jl/blob/9cdc063df740dd0579f12469e3e663935de2df0e/src/manifolds/Sphere.jl#L296-L309
  cosAngle <- drop(pmin(pmax(t(base_point) %*% target_point, -1), 1))
  res <- if(abs(cosAngle + 1) < tol){ # cosAngle \approx -1
    res <- matrix(0, nrow = length(base_point))
    if(abs(base_point[1] - 1) < tol && length(base_point) > 1){
      res[2] <- 1
    }else{
      res[1] <- 1
    }
    res <- res - drop(t(base_point) %*% res) * base_point
    res * pi / sqrt(sum(res^2))
  }else if(abs(cosAngle - 1) < tol){ # cosAngle \approx 1
    target_point
  }else{
    angle <- acos(cosAngle)
    (target_point - cosAngle * base_point) * angle / sin(angle)
  }
  sphere_project_tangent(res, base_point)
}


