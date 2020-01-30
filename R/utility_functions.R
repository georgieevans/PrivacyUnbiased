#' Bias Correct OLS
#'
#' This function takes OLS coefficients as inputs and bias corrects according to the noise matrix, S.
#' @param coef Vector of OLS coefficients
#' @param S Matrix of error terms. Dimensions must be length(coef) x length(coef)
#' @param X Design matrix that produced OLS coefficients
#' @param N Number of rows in X
#' @return
#' Returns a list:
#' \item{bias_correct}{Vector of bias corrected coefficients}
#' \item{x_prime_x}{X'X matrix to be used in variance estimationo}
#' @export
betaTilde <- function(coef, S, X, N){
  Q <- t(X)%*%X
  omega <- (1/N)*Q - S
  C <- solve(Q/N)%*%omega
  bias_correct  <- as.numeric(solve(C)%*%coef)
  names(bias_correct) <- names(coef)
  return(list(bias_correct = bias_correct, x_prime_x = Q))
}


#' Bias Correct OLS model with interaction term
#'
#' This function takes OLS coefficients from a moel including an
#' interaction term as inputs and bias corrects according to the noise matrix, S.
#' @param coef Vector of OLS coefficients
#' @param S Matrix of error terms. Dimensions must be length(coef) x length(coef)
#' @param X Design matrix that produced OLS coefficients
#' @param N Number of rows in X
#' @param int1 Column index of first interaction variable in X
#' @param int2 Column index of second interaction variable in X
#' @return
#' Returns a list:
#' \item{bias_correct}{Vector of bias corrected coefficients}
#' \item{x_prime_x}{X'X matrix to be used in variance estimationo}
#' @export
betaTildeInt <- function(coef, S, X, N, int1, int2){
  # Construct Omega
  Q <- t(X)%*%X
  omega <- (1/N)*Q - S
  int_col <- ncol(X)
  omega[int1, int_col] <- omega[int1, int_col] - S[int1, int1]*mean(X[, int2])
  omega[int_col, int1] <- omega[int1, int_col]
  omega[int2, int_col] <- omega[int2, int_col] - S[int2, int2]*mean(X[, int1])
  omega[int_col, int2] <- omega[int2, int_col]
  mu1 <- mean(X[, int1]^2) - S[int1, int1]
  mu2 <- mean(X[, int2]^2) - S[int2, int2]
  omega[int_col, int_col] <- omega[int_col, int_col] - S[int2, int2]*S[int1, int1] - S[int2, int2]*mean(mu1) - S[int1, int1]*mean(mu2)

  # Bias correct
  C <- solve(Q/N)%*%omega
  bias_correct  <- as.numeric(solve(C)%*%coef)
  names(bias_correct) <- names(coef)
  return(list(bias_correct = bias_correct, x_prime_x = Q))
}


#' Bias Correct OLS model with squared variable
#'
#' This function takes OLS coefficients from a model including
#' a squared variable as inputs and bias corrects according to the noise matrix, S
#' @param coef Vector of OLS coefficients
#' @param S Matrix of error terms. Dimensions must be length(coef) x length(coef)
#' @param X Design matrix that produced OLS coefficients
#' @param N Number of rows in X
#' @param index_sq Column index of X covariate that is squared
#' @return
#' Returns a list:
#' \item{bias_correct}{Vector of bias corrected coefficients}
#' \item{x_prime_x}{X'X matrix to be used in variance estimation}
#' @export
betaTildeSq <- function(coef, S, X, N, index_sq){

  # Construct X'Z
  Q <- t(X)%*%X
  X_Z <- Q - N*S
  K <- ncol(X)
  sq_col <- ncol(X)
  X_Z[sq_col, index_sq] <- X_Z[sq_col, index_sq] - 2*N*S[index_sq, index_sq]*mean(X[, index_sq])
  non_sqcol <- c(1:(index_sq-1), (index_sq + 1):(K-1))
  non_sqcol <- non_sqcol[-which(non_sqcol == sq_col|non_sqcol == index_sq)]
  for(i in non_sqcol){
    X_Z[i, sq_col] <- X_Z[i, sq_col] - N*S[index_sq, index_sq]*mean(X[, i])
  }
  X_Z[index_sq, sq_col] <- X_Z[index_sq, sq_col] - 3*N*S[index_sq, index_sq]*mean(X[, index_sq])
  mu2 <- mean(X[, index_sq]^2) - S[index_sq, index_sq]
  X_Z[sq_col, sq_col] <- X_Z[sq_col, sq_col] - 5*N*S[index_sq, index_sq]*mu2 - 3*N*(S[index_sq, index_sq])^2
  omega <- X_Z/N

  # Bias correct
  C <- solve(Q/N)%*%omega
  bias_correct  <- as.numeric(solve(C)%*%coef)
  names(bias_correct) <- names(coef)
  return(list(bias_correct = bias_correct, x_prime_x = Q))
}


#' Estimates \eqn{\sigma^2} from regression
#'
#' This function takes OLS coefficients and a Y vector to estimate \eqn{\sigma^2}
#' for the model.
#' @param Y Dependent variable from OLS model
#' @param bias_corrrect Bias corrected OLS coefficients
#' @param S Matrix of error terms. Dimensions must be length(coef) x length(coef)
#' @param X Design matrix that produced OLS coefficients
#' @return Estimate of \eqn{\sigma^2}
#' @export
sigmaSq <- function(Y, bias_correct, S, X){
  var(Y - as.numeric(bias_correct%*%t(X))) - t(bias_correct^2)%*%diag(S)
}

#' Estimates \eqn{Cov(X_k'y, X_j'X_m)}
#'
#' This function estimates components of variance-covariance matrix to subsequently
#' draw random variable of \eqn{\tilde{\beta}} from a multivariate normal.
#' @param k column index of first X.
#' @param j column index of second X.
#' @param m column index of third X.
#' @param S Matrix of error terms. Dimensions must be ncol(X) x ncol(X)
#' @param X Design matrix that produced OLS coefficients
#' @param Y Dependent variable from OLS model
#' @param N Number of rows in X
#' @return Estimate of Cov(X_k'y, X_j'X_m)
#' @export
covXyXX <- function(k, j, m, S, X, Y, N){
  cov_est <-  Y%*%X[, m]%*%S[k,j] + Y%*%X[, j]%*%S[k,m]
  return(cov_est)
}


#' Estimates \eqn{Cov(X_k'y, X_j'y)}
#'
#' This function estimates components of variance-covariance matrix to subsequently
#' draw random variable of \eqn{\tilde{\beta}} from a multivariate normal.
#' @param k column index of first X.
#' @param j column index of second X.
#' @param S Matrix of error terms. Dimensions must be ncol(X) x ncol(X)
#' @param X Design matrix that produced OLS coefficients
#' @param Y Dependent variable from OLS model
#' @param sigma_sq estimate of \Sigma^2
#' @param x_prime_x X'X
#' @param N Number of rows in X
#' @return Estimate of Cov(X_k'y, X_j'y)
#' @export
covXyXy <- function(k, j, S, X, Y, sigma_sq, x_prime_x, N){
  Z_prime_Z <- x_prime_x - N*S
  cov_est <- sigma_sq*Z_prime_Z[k,j] +  S[k, j]*(t(Y)%*%Y)
  return(cov_est)
}

#' Estimates \eqn{Cov(X_k'X_j, X_m'X_l)}
#'
#' This function estimates components of variance-covariance matrix to subsequently
#' draw random variable of \eqn{\tilde{\beta}} from a multivariate normal.
#' @param k column index of first X.
#' @param j column index of second X.
#' @param m column index of third X.
#' @param l column index of fourth X.
#' @param S Matrix of error terms. Dimensions must be ncol(X) x ncol(X)
#' @param X Design matrix that produced OLS coefficients
#' @param x_prime_x X'X
#' @param N Number of rows in X
#' @return Estimate of Cov(X_k'X_j, X_m'X_l)
#' @export
covXX <- function(k, j, l, m, S, X, x_prime_x, N){
  omega <- x_prime_x - N*S
  cov_est <- omega[k,l]*S[j, m] + omega[k, m]*S[j,l] + omega[j, l]*S[k,m] + omega[j, m]*S[k,l]  + 2*N*(S[k,l])*(S[j,m])
  return(cov_est)
}

#' Draw re-sampled estimates
#'
#' This function draws random variables vec(X'X, X'y) from a multivariate normal
#' and constructs a simulated \eqn{\tilde{\beta}} and b.
#' @param X Design matrix that produced OLS coefficients
#' @param mu Mean vector for multivariate normal
#' @param Sigma Variance covariance matrix for mjultivariate normal
#' @param S Matrix of error terms. Dimensions must be ncol(X) x ncol(X)
#' @param N Number of rows in X
#' @return One draw of vector of re-sample b and \eqn{\tilde{\beta}}
#' @export
simulateEst <- function(X, mu, Sigma, S, N){

  Sigma <- as.matrix(Matrix::nearPD(Sigma[-1, -1])$mat)
  draws <- c(mu[1], mvtnorm::rmvnorm(1, mean = mu[-1], sigma = Sigma)) # removing the column of 1's
  K <- ncol(X)

  # Construct X'X
  Q_star <- matrix(0, K, K)
  to <- length(draws) - ncol(X)
  draws_x_x <- draws[c(1:to)]
  upperTriangle(Q_star, diag = TRUE, byrow = TRUE) <- draws_x_x
  Q_star <- as.matrix(Matrix::forceSymmetric(Q_star, uplo="U"))

  # Construct X'y
  from <- length(draws) - ncol(X) + 1
  to <- length(draws)
  X_y <- draws[c(from:to)]

  # b star
  b_star <- as.numeric(solve(Q_star)%*%as.matrix(X_y))

  # beta tilde star
  omega_star <- (1/N)*Q_star - S
  C_star <- solve(Q_star/N)%*%omega_star
  bias_correct  <- as.numeric(solve(C_star)%*%b_star)

  # Bias correct
  bc_star <- as.numeric(solve(C_star)%*%b_star)
  bc_star <- c(b_star, bc_star)

  return(bc_star)
}

#' Simulate variance
#'
#' This function takes repeated draws from a MVN to simulate the approximate
#' distribution of \eqn{\tilde{\beta}} and b.
#' @param X Design matrix that produced OLS coefficients
#' @param Y Dependent variable from OLS model
#' @param sigma_sq Estimate of \eqn{\sigma^2}
#' @param nsims Number of times to draw from a MVN
#' @param coef Vector of OLS coefficientss
#' @param S Matrix of error terms. Dimensions must be length(coef) x length(coef)
#' @param N Number of rows in X
#' @param x_prime_x X'X
#' @param x_y X'y
#' @return
#' Returns a list:
#' \item{sims}{Draws of b and \eqn{\tilde{\beta}}}
#' \item{Sigma}{Estimate of variance-covariance that was used for simulation}
#' @export
varianceMVN <- function(X, Y, sigma_sq, nsims, coef, S, N, x_prime_x, x_y){

  # Extract upper triangle
  x_x <- unlist(lapply(1:nrow(x_prime_x), FUN = function(i) as.numeric(x_prime_x[i, c(i:ncol(x_prime_x))])))

  # Mean vector
  mu <- c(x_x, x_y)

  # Initialise matrix
  dim <- ncol(X)
  size <- (1/2)*dim*(dim+1) + dim
  threshold <- (1/2)*dim*(dim+1) + 1
  Sigma <- matrix(0, nrow = size, ncol = size)

  # Generate associated k's, j's, l's and m's (look-up table)
  k <- c(unlist(lapply(1:dim, FUN = function(i) rep(i, dim - i + 1))), 1:dim) #plus  X'y
  j <- unlist(lapply(1:dim, FUN = function(i) c(i:dim)))

  # m is first x, j is second x (for columns)
  m <- c(unlist(lapply(1:dim, FUN = function(i) rep(i, dim - i + 1))), 1:dim)
  l <- unlist(lapply(1:dim, FUN = function(i) c(i:dim)))


  # Build the covariance matrix
  for(r in 1:size){
    for(c in 1:size){
      if(r >= threshold & c < threshold){ # then we are in Cov(X'y, X'X) world
        Sigma[r, c] <-  as.numeric(covXyXX(k = k[r], j = m[c], m = l[c], S = S, X = X, Y = Y, N = N))
      }else{
        if(r < threshold & c >= threshold){ # then we are in Cov(X'X, X'y) world
          Sigma[r, c] <-  as.numeric(covXyXX(k = m[c], j = k[r], m = j[r], S = S, X = X, Y = Y, N = N))
        }else{
          if(r < threshold & c < threshold){ # then we are in cov(X'X, X'X) world
            Sigma[r, c] <- as.numeric(covXX(k = k[r], j =  j[r], m = m[c], l = l[c], S = S, X = X, x_prime_x, N = N))
          }else{ # then we are in cov(X'y, X'y) world
            Sigma[r, c] <- as.numeric(covXyXy(k = k[r], j = m[c], S = S, X = X, Y = Y, sigma_sq = sigma_sq,  x_prime_x, N = N))
          }
        }
      }
    }
  }
  # Now we have the covariance matrix and the mean vector; run simulations
  sims <- do.call(rbind, lapply(1:nsims, FUN = function(i) simulateEst(X = X, mu = mu, Sigma = Sigma, S = S, N = N)))

  return(list(sims = sims, Sigma = Sigma))
}


#' Bootstrap variance
#'
#' This function takes repeated samples from the data to estimate a bootstrap
#' distribution of \eqn{\tilde{\beta}} and b.
#' @param Y Dependent variable from OLS model
#' @param X Design matrix that produced OLS coefficients
#' @param S Matrix of error terms. Dimensions must be length(coef) x length(coef)
#' @param nsims Number of times to draw from a MVN
#' @param int_vars Vector of interaction variables
#' @param sq_vars Vector of squared variables
#' @param int1 Column index of first interaction variable
#' @param int2 Column index of second interaction variable
#' @param index_sq Column index of second interaction variable
#' @return Bootstrap samples of b and \eqn{\tilde{\beta}}
varianceBoot <- function(Y, X, S, nsims, N, int_vars, sq_vars, int1, int2, index_sq){

  boot_dist <- do.call(rbind, lapply(1:nsims, FUN = function(i){

    # Resample data
    sample_i <- sample(1:N, N, replace = TRUE)
    X_star <- X[sample_i, ]
    Y_star <- Y[sample_i]

    # Caculate estimates
    b_star <- as.numeric(coef(lm(Y_star ~ X_star[, -1])))

    if(length(int_vars) == 0 & length(sq_vars) == 0){
      bias_correct <- betaTilde(coef = b_star, S = S, X = X_star, N = N)
      bc_star <- bias_correct$bias_correct
    }else{
      if(length(int1) > 0){
        bias_correct <- betaTildeInt(coef = b_star, S = S, X = X_star, N = N, int1 = int1, int2 = int2)
        bc_star <- bias_correct$bias_correct
      }else{
        bias_correct <- betaTildeSq(coef = b_star, S = S, X = X_star, N = N, index_sq = index_sq)
        bc_star <- bias_correct$bias_correct
      }
    }

    return(c(b_star, bc_star))
  }))

  return(boot_dist)
}

