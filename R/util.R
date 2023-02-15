
randn <- function(n, m, ...){
  matrix(rnorm(n * m, ...), nrow = n, ncol = m)
}

sym <- function(M){
  0.5 * (M + t(M))
}

skew <- function(M){
  0.5 * (M - t(M))
}

DIM <- function(x){
  if(is.array(x)) dim(x)
  else length(x)
}


my_matrix_log <- function(M){
  suppressWarnings({
    logm <- expm::logm(M)
  })
  if(all(is.na(logm))){
    # The Higham08 algorithm failed. Try Eigen
    logm <- tryCatch({
      expm::logm(M, method = "Eigen")
    }, error = function(error){
      for(idx in 1:10){
        tmp <- tryCatch({
          expm::logm(M + randn(nrow(M), ncol(M), sd = 1e-12), method = "Eigen")
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
  logm
}
