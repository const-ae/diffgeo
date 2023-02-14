
randn <- function(n, m, ...){
  matrix(rnorm(n * m, ...), nrow = n, ncol = m)
}

sym <- function(M){
  0.5 * (M + t(M))
}

skew <- function(M){
  0.5 * (M - t(M))
}
