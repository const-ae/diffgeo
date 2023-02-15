

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
#' @export
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

#' @rdname sphere_injectivity_radius
#' @export
stiefel_injectivity_radius <- function(){
  0.89 * pi
}

#' Generate a random point on a Stiefel manifold
#'
#' The elements of a Grassmann manifold are represented
#' by matrices with orthogonal columns
#' (\eqn{\text{St}(n, k) := \{x \in \mathbb{R}^{n\times k} | x^T x = I\}}).
#'
#' @param n,k the dimensions of the Grassmann manifold.
#' @param ... additional parameters passed to [`rnorm`]
#'
#' @return a matrix with `n` rows and `k` orthogonal columns.
#'
#' @export
stiefel_random_point <- function(n, k, ...){
  x <- randn(n, k, ...)
  stiefel_project_point(x)
}

#' Generate a random tangent vector to a point on a Stiefel manifold
#'
#' The elements of the tangent space for a \eqn{p \in \text{St}(n,k)}
#' are \eqn{\mathcal{T}_p\text{St}(n,k) :=  \{v \in \mathbb{R}^n | p^Tv + v^Tp = 0\}}
#'
#' @param base_point the point from a Grassmann manifold. A matrix with orthogonal columns.
#' @param ... additional parameters passed to [`rnorm`]
#'
#' @return a matrix with `n` rows and `k` columns.
#'
#' @export
stiefel_random_tangent <- function(base_point, ...){
  n <- nrow(base_point)
  k <- ncol(base_point)
  v <- randn(n, k, ...)
  stiefel_project_tangent(v, base_point)
}


#' Project a point to (the tangent space of ) a Stiefel manifold
#'
#' @param x,v,base_point matrix
#'
#' @details The function uses the \eqn{Q} matrix of a QR-decomposition.
#'
#' @return the projected matrix
#'
#' @export
stiefel_project_point <- function(x){
  qr.Q(qr(x))
}

#' @rdname stiefel_project_point
#' @export
stiefel_project_tangent <- function(v, base_point){
  v - base_point %*% sym(t(base_point) %*% v)
}

#' Exponential map on a Stiefel manifold
#'
#' Go from the `base_point` in the direction `v`.
#'
#' @param v a tangent vector
#' @param base_point a point on a Stiefel manifold
#'
#' @details
#'   The exponential map on a Stiefel manifold for the canonical
#'   metric. For details please see Edelman et al. (1998).
#'
#' @references
#'   A. Edelman, T. A. Arias, and S. T. Smith. The geometry of algorithms with orthogonality constraints. SIAM Journal on Matrix Analysis and Applications, 20(2):303–353, 1998.
#'
#' @export
stiefel_map <- function(v, base_point){
  # Implementation based on https://github.com/JuliaManifolds/Manifolds.jl/blob/bf1859a860617eae62a209cf6ed03042cbaaa8fe/src/manifolds/StiefelCanonicalMetric.jl#L75
  # Note that the documentation of Julia Manifolds has a typo in the formula.
  k <- ncol(base_point)
  zeros <- matrix(0, nrow = k, ncol = k)
  ptv <- t(base_point) %*% v
  qr <- qr(v - base_point %*% ptv)
  block <- rbind(cbind(ptv, -t(qr.R(qr))),
                 cbind(qr.R(qr), zeros))
  BC_ext <- expm::expm(block)
  B <- BC_ext[seq_len(k),seq_len(k),drop=FALSE]
  C <- BC_ext[k+seq_len(k),seq_len(k),drop=FALSE]
  base_point %*% B + qr.Q(qr) %*% C
}

#' Find the tangent vector that connects two points on a Stiefel manifold
#'
#' @param base_point,target_point two points on a Stiefel manifold
#' @param tolerance,max_iter parameters to control the numerical approximation
#'   of the tangent vector calculation.
#'
#' @details
#'   For details see algorithm 4 from Zimmermann et al. (2022).
#'
#' @references
#'   Zimmermann, Ralf, and Knut Hüper. “Computing the Riemannian Logarithm on the Stiefel
#'   Manifold: Metrics, Methods and Performance.” arXiv, February 8, 2022. http://arxiv.org/abs/2103.12046.
#'
#' @export
stiefel_log <- function(base_point, target_point, tolerance = 1e-8, max_iter = 100){
# See https://github.com/RalfZimmermannSDU/RiemannStiefelLog/blob/main/Stiefel_log_general_metric/Matlab/Stiefel_Log.m
  n <- nrow(base_point)
  p <- ncol(base_point)
  zero <- matrix(0, nrow = p, ncol = p)

  M <- t(base_point) %*% target_point
  qn <- qr(target_point - base_point %*% M)
  # Orthogonal completion of [M N]' (step 3)
  Vk <- qr.Q(qr(rbind(M, qr.R(qn))), complete = TRUE)
  Vk[,seq_len(p)] <- rbind(M, qr.R(qn))
  if(Matrix::determinant(Vk)$sign != 1){
    # det == -1 --> flip one unimportant axis
    Vk[,p+1] <- -1 * Vk[,p+1]
  }

  for(k in seq(max_iter)){
    Vk <- rotation_project_point(Vk)
    log_Vk <- my_matrix_log(Vk)
    C <- log_Vk[p+seq_len(p),p+seq_len(p),drop=FALSE]
    if(sqrt(sum(C^2)) < tolerance){
      break
    }
    Phi_k <- expm::expm(-C)
    Wk <- rbind(cbind(diag(nrow = p), zero),
                cbind(zero, Phi_k))
    Vk <- Vk %*% Wk
  }
  base_point %*% log_Vk[seq_len(p),seq_len(p),drop=FALSE] + qr.Q(qn) %*% log_Vk[p+seq_len(p), seq_len(p),drop=FALSE]
}
