#' DP regression diagnostics
#'
#' Run diagnostic check on DP linear regression
#' @param lmdp_obj output from \code{\link{lmdp}}
#' @return Summary tesst of normality and homoskedasticity of errors
#' @export
#'
diagnosticsDP <- function(lmdp_obj)
{
  Y <- lmdp_obj$Y[-1]
  S_y <- lmdp_obj$Y[1]
  X <- lmdp_obj$X[-1, ]
  S <- lmdp_obj$S
  beta_tilde <- lmdp_obj$beta_tilde
  dp_res <- Y - X%*%beta_tilde
  e_hat <- (dp_res)^2 - as.numeric(t(beta_tilde)%*%S%*%beta_tilde) - S_y^2
  err_vec <- c(0, lmdp_obj$X[1, -1])
  dat <- as.data.frame(rbind(err_vec, cbind(e_hat = e_hat, X[, -1])))
  names(dat)[1] <- 'e_hat'
  formula <- as.formula(paste0('e_hat ~ ', strsplit(as.character(lmdp_obj$formula), '~')[[3]]))

  # Regression: Heteroskedasticity test
  het_reg <- lmdp(formula, data = dat)

  # Testing normality of errors
  S_e <- sqrt(as.numeric(S_y^2 + t(beta_tilde)%*%S%*%beta_tilde))
  e_moments <- Rmoments(dp_res, S_e, 4, n = nrow(X))
  normality <- as.data.frame(cbind(`Est. error` = scaledMoments(e_moments[, 1])[3:4], `Normal dist.`= c(0, 3)))
  rownames(normality) <- c('Skewness', 'Kurtosis')

  cat("\n Heteroskedsaticity test via Variance Regression (Bias Corrected Residuals):\n\n")
  print.default(summary(het_reg),
                print.gap = 2L, quote = FALSE)
  cat("\n Error Normality Test:\n\n")
  if(sum(e_moments[3:4, 1]/e_moments[3:4, 2] < 2) > 0){
    cat(' Warning: Error moment estimates are imprecise - unable to accurately test normality of errors')
  }else{
    print(normality)
  }
  output <- list(error_moments = e_moments,
                 heteroskedasticity_reg = summary(het_reg),
                 normality_test = normality)

  return(invisible(output))

}




