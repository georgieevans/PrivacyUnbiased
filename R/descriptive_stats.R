#' Estimate the Rth raw moment
#'
#' This function estimates the Rth moment of a column of the private data according to the DP data and the noise matrix, S.
#' @param X Column from DP data
#' @param S Scalar value of DP error.
#' @param r Moment to be estimate
#' @return Estimate of Rth moment
#' @export
momentR <- function(X, S, r){
  fn <- mpoly::hermite(r)
  herm <- function(x){
    suppressMessages((as.function(fn)(x)))
  }
  herm <- herm(X/S)
  herm_poly <- (S^r)*herm
  return(c(mean(herm_poly), mean(herm^2)))
}

#' Estimate lagged value of column
#'
#' This function estimates first R moments of a column of the private data according to the DP data and the noise matrix, S.
#' @param x vector to lag
#' @param k length of lag
#' @return vector of lags
#' @export
lagHerm <- function(x, k){
  if (k>0){
    return (c(rep(NA, k), x)[1:length(x)] )
  }
  else{
    return(c(x[(-k+1):length(x)], rep(NA, -k)))
  }
}


#' Estimate the first R raw moments
#'
#' This function estimates first R moments of a column of the private data according to the DP data and the noise matrix, S.
#' @param X Column from DP data
#' @param S Scalar value of DP error.
#' @param R First R moments to be estimated
#' @return Vector of estimates of first R moments
#' @export
Rmoments <- function(X, S, R, n){
  df <- as.data.frame(do.call(rbind, lapply(1:R, function(r){
    return(momentR(X, S, r))
  })))
  rows <- 1:nrow(df)
  df[, 2] <- ((4*rows^2*S^(2*rows))/n)*lagHerm(df[ ,2], 1)
  df[1, 2] <- S^2/n^2
  df[,2] <- sqrt(df[,2])
  names(df) <- c('raw_moment', 'se')
  return(round(df, 6))
}


#' Estimate the first R raw moments
#'
#' This function estimates descriptive statistics of private data from column in DP data
#' @param variable Variable in data to estimate descriptive statistics of
#' @param data Data
#' @return Table containing descriptive stats of private data column (estimated) and DP data column (observed)
#' @export
descriptiveDP <- function(variable, data){
  var <- eval(substitute(variable),data, parent.frame())
  S <- var[1]
  X <- var[-1]
  moments <- Rmoments(X, S, 4)[, 1]
  m1 <- moments[1]
  m2 <- moments[2]
  m3 <- moments[3]
  m4 <- moments[4]
  var_z <- -m1^2 + m2
  sd_z <- sqrt(var_z)
  skew_z <- (m3 - 3*m1*sd_z^2 - m1^3)/(sd_z^3)
  kurt_z <- (-3*m1^4 + 6*m1^2*m2 - 4*m1*m3 + m4)/(sd_z^4)

  Z_df <- data.frame(Mean = m1, `Std Dev` = sqrt(var_z), Skewness = skew_z, Kurtosis = kurt_z)
  X_df <- data.frame(Mean = mean(X), `Std Dev` = sd(X), Skewness = moments::skewness(X),
                     Kurtosis = moments::kurtosis(X))
  df <- round(rbind(Z_df, X_df), 4)
  rownames(df) <- c('Before DP noise (estimated)', 'After DP noise (observed)')
  return(df)
}



piZNB <- function(moments){
  (moments[1]^2*moments[2] + moments[1]*(moments[2] + moments[3]) -
     2*moments[2]^2 - moments[1]^3)/(moments[1]*(moments[2] + moments[3]) - 2*moments[2]^2)
}

pZNB <- function(moments){
  (moments[1]*(moments[2] - moments[3]) + moments[2]^2 - moments[1]^2)/
    (moments[2]^2 - moments[1]*moments[3])
}

rZNB <- function(moments){
  (2*moments[2]^2 - moments[1]*(moments[2] + moments[3]))/
    (moments[1]^2 + moments[1]*(moments[3] - moments[2]) - moments[2]^2)
}

estZNB <- function(moments){
  pi <- piZNB(moments)
  p <- pZNB(moments)
  r <- rZNB(moments)
  return(c(pi = pi, p = p, r = r))
}

estZIP <- function(moments){
  pi <- 1 - moments[2]/moments[1]
  lambda <- (moments[2] - moments[1])/moments[1]
  return(c(pi = pi, lambda = lambda))
}

# Solve for params ZIP
paramsZIP <- function(X, S, R = 10, n){
  moments_df <- Rmoments(X, S, R, n = n)
  moments <- moments_df[,1]
  est <- as.numeric(estZIP(moments))
  if(est[1] < 0| est[1] > 1){
     stop(message('Estimates outside logical bounds'))
    return(list(params = est))
  }
  ind <- sample(c(0, 1), n, prob = c(est[1], 1 - est[1]), replace = TRUE)
  Z_est <- ind*rpois(n, est[2])
  implied_moments <- sapply(1:length(moments), FUN = function(r) mean(Z_est^r))

  return(list(params = est,
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2)))

}


# Solve for params ZNB
paramsZINB <- function(X, S, R = 10, n){
  moments_df <- Rmoments(X, S, R, n = n)
  moments <- moments_df[,1]
  est <- as.numeric(estZNB(moments))
  if(est[2] < 0| est[2] > 1| est[1] < 0| est[1] > 1){
    stop(message('Estimates outside logical bounds'))
  }
  ind <- sample(c(0, 1), 2*n, prob = c(est[1], 1 - est[1]), replace = TRUE)
  Z_est <- ind*rnbinom(2*n, mu = est[2]*est[3]/(1-est[2]), size = est[3])
  implied_moments <- sapply(1:length(moments), FUN = function(r) mean(Z_est^r))

  return(list(params = est,
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2)))

}


paramsNormal <- function(X, S, R = 10, n){
  moments_df <- Rmoments(X, S, R, n = n)
  moments <- moments_df[, 1]
  mu <- mean(X)
  var <- var(X) - S^2
  Z_est <- rnorm(n, mean = mu, sd = sqrt(var))
  implied_moments <- sapply(1:nrow(moments_df), FUN = function(r) mean(Z_est^r))

  return(list(params = c(mu, var),
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2)))
}


paramsPoisson <- function(X, S, R = 10, n){
  moments_df <- Rmoments(X, S, R, n = n)
  moments <- moments_df[, 1]
  lambda <- mean(X)
  Z_est <- rpois(n, lambda)
  implied_moments <- sapply(1:nrow(moments_df), FUN = function(r) mean(Z_est^r))

  return(list(params = c(lambda),
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2)))

}

paramsNB <- function(X, S, R = 10, n){
  moments_df <- Rmoments(X, S, R, n = n)
  moments <- moments_df[, 1]
  var <- var(X) - S^2
  p <- 1 - (mean(X)/var)
  r <- (-mean(X)^2)/(mean(X) - var)
  Z_est <- rnbinom(n, size = r, mu = p*r/(1-p))
  implied_moments <- sapply(1:nrow(moments_df), FUN = function(r) mean(Z_est^r))
  return(list(params = c(p, r),
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2)))

}


