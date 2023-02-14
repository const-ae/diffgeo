
#' Check if a point is from (the tangent space of) a Grassmann manifold
#'
#' @param x point for which we check if it is on a Grassmann manifold
#' @param v,base_point tangent vector for we which we check if is in the
#'   tangent space of a Grassmann manifold at `base_point`
#' @param tol the tolerance for the test
#' @param error a flag to indicate if the function should
#'   throw an error if the check fails
#'
#' @return
#'  `TRUE` if the check succeeds, `FALSE` if the check fails (
#'  and `error = FALSE`).
#'
#' @export
grassmann_check_point <- function(x, tol = 1e-12, error = TRUE){
  stiefel_check_point(x, tol = tol, error = error)
}

#' @rdname grassmann_check_point
grassmann_check_tangent <- function(v, base_point, tol = 1e-12, error = TRUE){
  stiefel_check_tangent(v, base_point, tol = tol, error = error)
}

#' @rdname sphere_injective_radius
grassmann_injective_radius <- function(n){
  pi / 2
}

#' Generate a random point on a Grassmann manifold
#'
#' The elements of a Grassmann manifold are represented
#' by matrices with orthogonal columns
#' (\eqn{\text{Gr}(n, k) := \{\text{span}(x): x \in \mathbb{R}^{n\times k} | x^T x = I\}}).
#'
#' @param n,k the dimensions of the Grassmann manifold.
#'
#' @return a matrix with `n` rows and `k` orthogonal columns.
#'
#' @export
grassmann_random_point <- function(n, k, ...){
  x <- randn(n, k, ...)
  grassmann_project_point(x)
}

#' Generate a random tangent vector to a point on a Grassmann manifold
#'
#' The elements of the tangent space for a \eqn{p \in \text{Gr}(n,k)}
#' are \eqn{\mathcal{T}_p\text{Gr}(n,k) :=  \{v \in \mathbb{R}^n | p^Tv + v^Tp = 0\}}
#'
#' @param base_point the point from a Grassmann manifold. A matrix with orthogonal columns.
#'
#' @return a matrix with `n` rows and `k` columns.
#'
#' @export
grassmann_random_tangent <- function(base_point, ...){
  n <- nrow(base_point)
  k <- ncol(base_point)
  v <- randn(n, k, ...)
  grassmann_project_tangent(v, base_point)
}

#' Project a point to (the tangent space of ) a Grassmann manifold
#'
#' @param x,v,base_point matrix
#'
#' @return the projected matrix
#'
#' @export
grassmann_project_point <- function(x){
  project_stiefel(x)
}

#' @rdname grassmann_project_point
grassmann_project_tangent <- function(x, base_point){
  x - base_point %*% t(base_point) %*% x
}


#' Exponential map on a Grassmann manifold
#'
#' Go from the `base_point` in the direction `v`.
#'
#' @param v a tangent vector
#' @param base_point a point on a Grassmann manifold
#'
#' @details
#'   The exponential map on a Grassmann manifold is
#'     \deqn{\exp_p(v) = p V \text{diag}(\cos(d)) V^T + U \text{diag}(\sin(d)) V^T,}
#'   where \eqn{U \text{diag}(d) V^T = v} is the SVD of the tangent vector \eqn{v} and \eqn{p}
#'   is the `base_point`.
#'
#' @export
grassmann_map <- function(x, base_point){
  # Adapted from https://github.com/JuliaManifolds/Manifolds.jl/blob/master/src/manifolds/GrassmannStiefel.jl#L93
  if(ncol(base_point) == 0 || nrow(base_point) == 0){
    base_point
  }else if(any(is.na(x))){
    matrix(NA, nrow = nrow(x), ncol = ncol(x))
  }else{
    svd <- svd(x)
    z <- base_point %*% svd$v %*% diag(cos(svd$d), nrow = length(svd$d)) %*% t(svd$v) +
      svd$u %*% diag(sin(svd$d), nrow = length(svd$d)) %*% t(svd$v)
    z
  }
}


#' Find the tangent vector that connects two points on a Grassmann manifold
#'
#' @param base_point,target_point two points on a Grassmann manifold
#'
#' @details
#'   The inverse of the exponential map on a Grassmann manifold is
#'     \deqn{\log(p, q) = V \text{diag}(\tan^{-1}(d))U^T}
#'   where \eqn{U\text{diag}(d)V^T = (q^Tp)^{-1} (q - q^Tpp^T)} is a SVD,
#'   \eqn{p} is the `base_point` and \eqn{q} is the `target_point`.
#'
#' @export
grassmann_log <- function(base_point, target_point){
  # Adapted from https://github.com/JuliaManifolds/Manifolds.jl/blob/master/src/manifolds/GrassmannStiefel.jl#L174
  # The Grassmann manifold handbook proposes an alternative algorithm in section 5.2
  n <- nrow(base_point)
  k <- ncol(base_point)
  stopifnot(nrow(target_point) == n, ncol(target_point) == k)
  if(n == 0 || k == 0){
    base_point
  }else{
    z <- t(target_point) %*% base_point
    At <- t(target_point) - z %*% t(base_point)
    Bt <- lm.fit(z, At)$coefficients
    svd <- svd(t(Bt), k, k)
    svd$u %*% diag(atan(svd$d), nrow = k) %*% t(svd$v)
  }
}


#' Find the rotation angle for between elements of Grassmann manifolds
#'
#' @param v tanget vector
#' @param base_point,target_point two points from a Grassmann manifold
#'
#' @references
#'   Bendokat, Thomas, Ralf Zimmermann, and P-A. Absil. "A Grassmann manifold handbook: Basic geometry and computational aspects." arXiv preprint arXiv:2011.13699 (2020).
#'
#' @export
grassmann_angle_from_tangent <- function(v, normalized = TRUE){
  # Conversion of tangent vector to angle taken from Proposition 5.1 of the Grassmann manifold handbook
  thetas <- (svd(v)$d / pi * 180)
  if(normalized){
    thetas <- thetas %% 180
    max(pmin(thetas, 180 - thetas))
  }else{
    thetas[1]
  }
}

#' @rdname grassmann_angle_from_tangent
grassmann_angle_from_points <- function(base_point, target_point){
  grassmann_angle_from_tangent(grassmann_log(base_point, target_point))
}

