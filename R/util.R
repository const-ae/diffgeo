
randn <- function(n, m, ...){
  matrix(rnorm(n * m, ...), nrow = n, ncol = m)
}
