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


#' Estimate the Scaled Moments
#'
#' This function estimates mean, variance, skew and kurtosis of the private data
#' @param moments Vector of estimated raw moments
#' @return Vector of scaled moments
#' @export
scaledMoments <- function(moments){
  m1 <- moments[1]
  m2 <- moments[2]
  m3 <- moments[3]
  m4 <- moments[4]
  var_z <- -m1^2 + m2
  sd_z <- sqrt(var_z)
  skew_z <- (m3 - 3*m1*sd_z^2 - m1^3)/(sd_z^3)
  kurt_z <- (-3*m1^4 + 6*m1^2*m2 - 4*m1*m3 + m4)/(sd_z^4)
  return(c(m1, var_z, skew_z, kurt_z))
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

