
herteroskedasticityDP <- function(lmdp_obj)
  {
  Y <- lmdp_obj$Y[-1]
  S_y <- lmdp_obj$Y[1]
  X <- lmdp_obj$X[-1, ]
  S <- lmdp_obj$S
  beta_tilde <- lmdp_obj$beta_tilde
  e_hat <- (Y - X%*%beta_tilde)^2 - as.numeric(t(beta_tilde)%*%S%*%beta_tilde) - S_y
  err_vec <- c(0, lmdp_obj$X[1, -1])
  dat <- as.data.frame(rbind(err_vec, cbind(e_hat = e_hat, X[, -1])))
  names(dat)[1] <- 'e_hat'
  formula <- as.formula(paste0('e_hat ~ ', strsplit(as.character(lmdp_obj$formula), '~')[[3]]))
  het_reg <- lmdp(formula, data = dat)
  return(list(heteroskedasticity_test = summary(het_reg)))
  }








