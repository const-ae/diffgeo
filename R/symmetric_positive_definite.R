
#' Check if a point is from (the tangent space of) a symmetric positive definite (SPD) manifold
#'
#' @param x point for which we check if it is on an SPD manifold
#' @param v,base_point tangent vector for we which we check if is in the
#'   tangent space of an SPD manifold at `base_point`
#' @param tol the tolerance for the test
#' @param error a flag to indicate if the function should
#'   throw an error if the check fails
#'
#' @return
#'  `TRUE` if the check succeeds, `FALSE` if the check fails (
#'  and `error = FALSE`).
#'
#' @export
spd_check_point <- function(x, tol = 1e-12, error = TRUE){
  if(! is.matrix(x) || ncol(x) != nrow(x)){
    if(error) stop("A SPD point must be specified by a square matrix")
    FALSE
  }
  if(! all(eigen(x)$values >= -tol)){
    if(error) stop("The eigenvalues of a SPD point must be all positive")
    FALSE
  }
  TRUE
}

#' @rdname spd_check_point
spd_check_tangent <- function(v, base_point, tol = 1e-12, error = TRUE){
  if(! is.matrix(v) || ncol(v) != nrow(v)){
    if(error) stop("A SPD tangent must be specified by a square matrix")
    FALSE
  }
  if(sum((v - t(v))^2) > tol){
    if(error) stop("'The SPD vector must be symmetric ('v - t(v) == 0')")
    FALSE
  }
  TRUE
}

#' @rdname sphere_injectivity_radius
spd_injectivity_radius <- function(){
  Inf
}


#' Generate a random point on an SPD manifold
#'
#' The elements of an SPD manifold are
#' \eqn{\text{SPD}(n) := \{x \in \mathbb{R}^{n\times n} | x = x^T,\, a^T x a > 0 \text{ for all } a \in \mathbb{R}^n\}}
#'
#' @param n the dimension of the SPD manifold.
#'
#' @return a symmetric positive definite matrix
#'
#' @export
spd_random_point <- function(n, ...){
  vec <- qr.Q(qr(randn(n, n)))
  val <- rlnorm(n, meanlog = 0, sd = 0.1)
  vec %*% diag(val, nrow = n) %*% t(vec)
}


#' Generate a random tangent vector to a point on a rotation manifold
#'
#' The elements of the tangent space for a \eqn{p \in \text{SPD}(n)}
#' are skew-symmetric matrices (\eqn{\mathcal{T}_p\text{SPD}(n) :=  \{v \in \mathbb{R}^{n\times n} | v = v^T\}})
#'
#' @param base_point the point from a rotation manifold. A square matrix.
#'
#' @return a skew-symmetruc matrix.
#'
#' @export
spd_random_tangent <- function(base_point, ...){
  n <- nrow(base_point)
  stopifnot(n == ncol(base_point))
  v <- randn(n, n, ...)
  spd_project_tangent(v, base_point)
}

#' Project a point to (the tangent space of ) an SPD manifold
#'
#' @param x,v,base_point square matrix
#' @param min_eigen_value the minimum eigenvalue of the matrix.
#'   Note that making this value too small leads to problems
#'   when calculating the inverse.
#'
#' @return the projected matrix
#'
#' @references
#'   Hickham blog post (2021), https://nhigham.com/2021/01/26/what-is-the-nearest-positive-semidefinite-matrix/
#'
#' @export
spd_project_point <- function(x, min_eigen_value = 0.01){
  if(nrow(x) == 0){
    x
  }else{
    eigen <- eigen(sym(x))
    eigen$vectors %*% diag(pmax(min_eigen_value, eigen$values), nrow = nrow(x)) %*% t(eigen$vectors)
  }
}

spd_project_tangent <- function(v, base_point = NULL){
  sym(v)
}


#' Exponential map on an SPD manifold
#'
#' Go from the `base_point` in the direction `v`.
#'
#' @param v a tangent vector
#' @param base_point a point on a SPD manifold
#'
#' @details
#'   The exponential map on a rotation manifold is
#'     \deqn{\exp_p(v) = p^{1/2} \exp(p^{-1/2} v p^{-1/2}) p^{1/2}}
#'   where \eqn{\exp(v)} is the matrix exponential and \eqn{p} is
#'   the `base_point`
#'
#' @export
spd_map <- function(v, base_point){
  ps <- spd_sqrt(base_point)
  psi <- spd_sqrt_inv(base_point)

  ps %*% expm::expm(psi %*% v %*% psi) %*% ps
}

#' Find the tangent vector that connects two points on a rotation manifold
#'
#' @param base_point,target_point two points on a rotation manifold
#'
#' @details
#'   The inverse of the exponential map on a rotation manifold is
#'     \deqn{\log(p, q) = p^{1/2} \log(p^{-1/2} q p^{-1/2}) p^{1/2}}
#'   where \eqn{\log(x)} is the matrix logarithm,
#'   \eqn{p} is the `base_point` and \eqn{q} is the `target_point`.
#'
#' @export
spd_log <- function(base_point, target_point){
  ps <- spd_sqrt(base_point)
  psi <- spd_sqrt_inv(base_point)

  ps %*% expm::logm(psi %*% target_point %*% psi) %*% ps
}

spd_sqrt <- function(x){
  stopifnot(nrow(x) == ncol(x))
  if(nrow(x) == 0){
    x
  }else{
    e <- eigen(sym(x))
    sqrt_s <- pmax(sqrt(e$values), .Machine$double.xmin)
    sym(e$vectors %*% (sqrt_s * t(e$vectors)))
  }
}

spd_sqrt_inv <- function(x){
  stopifnot(nrow(x) == ncol(x))
  if(nrow(x) == 0){
    x
  }else{
    e <- eigen(sym(x))
    sqrt_s <- 1/pmax(sqrt(e$values), .Machine$double.xmin)
    sym(e$vectors %*% (sqrt_s * t(e$vectors)))
  }
}

