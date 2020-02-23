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
  herm_poly <- (S^r)*herm(X/S)
  return(mean(herm_poly))
}

#' Estimate the first R raw moments
#'
#' This function estimates first R moments of a column of the private data according to the DP data and the noise matrix, S.
#' @param X Column from DP data
#' @param S Scalar value of DP error.
#' @param R First R moments to be estimated
#' @return Vector of estimates of first R moments
#' @export
Rmoments <- function(X, S, R){
   sapply(1:R, function(r){
    return(momentR(X, S, r))
  })
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
  moments <- Rmoments(X, S, 4)
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


