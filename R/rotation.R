
#' Check if a point is from (the tangent space of) a Rotation Manifold
#'
#' @param x point for which we check if it is on a rotation manifold
#' @param v,base_point tangent vector for we which we check if is in the
#'   tangent space of a rotation manifold at `base_point`
#' @param tol the tolerance for the test
#' @param error a flag to indicate if the function should
#'   throw an error if the check fails
#'
#' @return
#'  `TRUE` if the check succeeds, `FALSE` if the check fails (
#'  and `error = FALSE`).
#'
#' @export
rotation_check_point <- function(x, tol = 1e-12, error = TRUE){
  if(! is.matrix(x) || ncol(x) != nrow(x)){
    if(error) stop("A rotation point must be specified by a square matrix")
    FALSE
  }
  xtx <- t(x) %*% x
  if(sum((xtx - diag(nrow = nrow(x)))^2) > tol){
    if(error) stop("The crossproduct 't(x) %*% x' must equal the identity matrix")
    FALSE
  }
  det <- Matrix::determinant(x)
  if(abs(det$modulus) > tol && det$sign != 1){
    if(error) stop("The determinant of 'x' must be '1'.")
    FALSE
  }
  TRUE
}

#' @rdname rotation_check_point
rotation_check_tangent <- function(v, base_point, tol = 1e-12, error = TRUE){
  if(! is.matrix(v) || ncol(v) != nrow(v)){
    if(error) stop("A rotation tangent must be specified by a square matrix")
    FALSE
  }
  if(sum((v + t(v))^2) > tol){
    if(error) stop("'The tangent vector must be skew-symmetric ('v + t(v) == 0')")
    FALSE
  }
  TRUE
}

#' @rdname sphere_injectivity_radius
rotation_injectivity_radius <- function(){
  pi / sqrt(2)
}


#' Generate a random point on a rotation manifold (SO(n))
#'
#' The elements of a rotation manifold (also called special orthogonal) are
#' \eqn{\text{SO}(n) := \{x \in \mathbb{R}^{n\times n} | x^Tx = xx^T = I_n, \text{det}(x) = 1\}}
#'
#' @param n the dimension of the rotation manifold.
#'
#' @return an orthogonal matrix with `Matrix::det(x) == 1`
#'
#' @export
rotation_random_point <- function(n, ...){
  x <- randn(n, n, ...)
  rotation_project_point(x)
}

#' Generate a random tangent vector to a point on a rotation manifold
#'
#' The elements of the tangent space for a \eqn{p \in \text{SO}(n)}
#' are skew-symmetric matrices (\eqn{\mathcal{T}_p\text{SO}(n) :=  \{v \in \mathbb{R}^{n\times n} | w = pv, v + v^T = 0\}})
#'
#' @param base_point the point from a rotation manifold. A square matrix.
#'
#' @return a skew-symmetruc matrix.
#'
#' @export
rotation_random_tangent <- function(base_point, ...){
  n <- nrow(base_point)
  stopifnot(n == ncol(base_point))
  Z <- randn(n, n, ...)
  rotation_project_tangent(Z, base_point)
}


#' Project a point to (the tangent space of ) a rotation manifold
#'
#' @param x,v,base_point square matrix
#'
#' @return the projected matrix
#'
#' @export
rotation_project_point <- function(x){
  if(nrow(x) == 0){
    x
  }else{
    svd <- svd(x)
    diag_elem <- c(rep(1, times = ncol(x) - 1), Matrix::det(svd$u %*% t(svd$v)))
    svd$u %*% diag(diag_elem, nrow = length(diag_elem)) %*% t(svd$v)
  }
}

#' @rdname rotation_project_point
rotation_project_tangent <- function(v, base_point){
  skew(v)
}


#' Exponential map on a rotation manifold
#'
#' Go from the `base_point` in the direction `v`.
#'
#' @param v a tangent vector
#' @param base_point a point on a rotation manifold
#'
#' @details
#'   The exponential map on a rotation manifold is
#'     \deqn{\exp_p(v) = p \exp(v)}
#'   where \eqn{\exp(v)} is the matrix exponential and \eqn{p} is
#'   the `base_point`
#'
#' @export
rotation_map <- function(v, base_point){
  base_point %*% expm::expm(v)
}

#' Find the tangent vector that connects two points on a rotation manifold
#'
#' @param base_point,target_point two points on a rotation manifold
#'
#' @details
#'   The inverse of the exponential map on a rotation manifold is
#'     \deqn{\log(p, q) = \log(p^Tq)}
#'   where \eqn{\log(x)} is the matrix logarithm,
#'   \eqn{p} is the `base_point` and \eqn{q} is the `target_point`.
#'
#' @export
rotation_log <- function(base_point, target_point){
  if(nrow(base_point) == 0){
    base_point
  }else{
    suppressWarnings({
      logm <- expm::logm(t(base_point) %*% target_point)
    })
    if(all(is.na(logm))){
      # The Higham08 algorithm failed. Try Eigen
      logm <- tryCatch({
        expm::logm(t(base_point) %*% target_point, method = "Eigen")
      }, error = function(error){
        for(idx in 1:10){
          tmp <- tryCatch({
            expm::logm(t(base_point) %*% target_point + randn(nrow(base_point), ncol(base_point), sd = 1e-12), method = "Eigen")
          }, error = function(error2) NULL)
          if(is.null(tmp)){
            # try once more
          }else{
            break
          }
        }
        tmp
      })
      if(is.null(logm)){
        stop("Cannot calculate the matrix logarithm. It may not exist")
      }
    }

    skew(logm)
  }
}




