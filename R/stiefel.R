

#' Check if a point is from (the tangent space of) a Stiefel manifold
#'
#' @param x point for which we check if it is on a Stiefel manifold
#' @param v,base_point tangent vector for we which we check if is in the
#'   tangent space of a Stiefel manifold at `base_point`
#' @param tol the tolerance for the test
#' @param error a flag to indicate if the function should
#'   throw an error if the check fails
#'
#' @return
#'  `TRUE` if the check succeeds, `FALSE` if the check fails (
#'  and `error = FALSE`).
#'
#' @export
stiefel_check_point <- function(x, tol = 1e-12, error = TRUE){
  if(! is.matrix(x) ){
    if(error) stop("A point on a Stiefel manifold must be specified by a matrix")
    FALSE
  }

  if(sum((t(x) %*% x - diag(nrow = ncol(x)))^2) > tol){
    if(error) stop("'t(x) %*% x' must equal the identity matrix.")
    FALSE
  }
  TRUE
}

#' @rdname stiefel_check_point
stiefel_check_tangent <- function(v, base_point, tol = 1e-12, error = TRUE){
  if(! is.matrix(v) ){
    if(error) stop("A Stiefel tangent must be specified by a matrix")
    FALSE
  }
  error <- t(base_point) %*% v + t(v) %*% base_point
  if(sum(error^2) > tol){
    if(error) stop("'sum(base_point * v)' must be zero.")
    FALSE
  }
  TRUE
}
