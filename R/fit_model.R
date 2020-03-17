#' OLS differential privacy
#'
#' This function works similarly to \code{\link{lm}}. It takes a formula and data and
#' returns an lmdp object containing bias corrected OLS coefficients and standard errors. The output can be summarised
#' by inputting it as an argument to \code{\link{summary.lmdp}}. See an overview at \url{bit.ly/PrivUexample}.
#'
#' @param formula An \code{\link{lm}} style formula.
#' @param data The data to estimate the model on. The first row of data should contain the DP
#' standard error associated with that column unless noise argument is not NULL.
#' @param bootstrap_var If FALSE, then the variance is estimated via simulation.
#' If TRUE then the variance covariance matrix is estimated via bootstrap methods. Default is
#' FALSE unless model contains interaction terms/squared terms or fewer than 10000 obsesrvations
#' @param nsims_var Number of bootstrap samples/simulations. Default is 500
#' @param noise Set a default differentially private standard error for every column of the data matrix
#' @return
#' Returns an object of class lmdp containing:
#' \item{b}{Inconsistent OLS coefficient estimate}
#' \item{b_vcov}{Estimate of variance covariance of b}
#' \item{beta_tilde}{Consistent estimates of coefficients, \eqn{\tilde{\beta}}}
#' \item{beta_tilde_vcov}{Estimate of variance covariance of \eqn{\tilde{\beta}}}
#' \item{var_sims}{Full set of simulated/bootstrap estimates of b and \eqn{\tilde{\beta}}}
#' \item{Sigma_sq_hat}{Estimate of \eqn{\sigma^2}}
#' \item{vc_pos_def}{Indicator variable = 1 if covariance estimate was PD. NA if bootstrap used}
#' \item{boot}{Indicator variable = 1 if bootstrap was used to estimate variance}
#' \item{est_vc}{Variance-covariance matrix used in variance simulation. NA if bootstrap used}
#' \item{X}{Matrix of covariates}
#' \item{Y}{Dependent variable vector}
#' \item{formula}{Model formula}
#' @export
#' @examples
#' \dontrun{data(dp_data)}
#' \dontrun{lmdp_test <- lmdp(Y ~ X1 + X2 + X3, data = dp_data)}
#' \dontrun{summary(lmdp_test)}

lmdp <- function(formula, data, bootstrap_var = FALSE, nsims_var = 500, noise = NULL)
  {

  # Construct S vector and remove error row if included
  if(is.null(noise)){
  S_vec <- as.numeric(data[1, ])^2
  data <- data[-1, ]
  }else{
    S_vec <- rep(noise^2, ncol(data))
  }

  N <- nrow(data)

  # Run OLS
  reg <- lm(formula, data = data)
  b <- coef(reg)

  # Extract X matrix and Y vector
  X <- as.matrix(model.matrix(reg))
  Y <- as.numeric(model.frame(reg)[, 1])

  # Save inconsistent OLS variance
  names(b) <- colnames(X)
  b_var <- vcov(reg)

  # Testing for interactions/ squared variables
  int_vars <- grep(':', colnames(X), value=TRUE)
  sq_vars <-  grep("\\^2", colnames(X), value = TRUE)
  int1 <- NA
  int2 <- NA
  index_sq <- NA

  if(length(int_vars) > 1){
    stop("Model cannot contain more than one interaction variable")
  }

  if(length(sq_vars) > 1){
    stop("Model cannot contain more than one squared variable")
  }

  if(length(int_vars) > 0 & length(sq_vars) > 0){
    stop("Model cannot contain interaction and squared variable")
  }

  boot <- TRUE
  vc_pos_def <- NA
  est_vc <- NA

  # No squared or interaction terms

  if(length(int_vars) == 0 & length(sq_vars) == 0){

    # Construct S matrix
    S <- matrix(0, nrow = ncol(X), ncol = ncol(X))

    # Find position of relevant columns
    err_ind <- match(colnames(X)[-1], colnames(data))
    diag(S) <- c(0, S_vec[err_ind])

    # Bias correct

    bias_correct <- betaTilde(coef = b, S = S, X = X, N = N)
    beta_tilde <- bias_correct$bias_correct
    x_prime_x <- bias_correct$x_prime_x

    # Estimate sigma^2
    sigma_sq <- sigmaSq(Y = Y, bias_correct = beta_tilde, S = S, X = X)

    if(sigma_sq < 0){
      warning("Sigma squared estimate negative")
    }

    # Estimate variance
    b_index <- 1:length(b)
    beta_tilde_index <- (length(b) + 1):(2*length(b))

    if(bootstrap_var == FALSE & N > 9999){

      x_y <- t(X)%*%Y

      var_sims <-  varianceMVN(X = X, Y = Y, sigma_sq = sigma_sq,
                               nsims = nsims_var, coef = b, S = S, N = N,
                               x_prime_x = x_prime_x, x_y = x_y)

      # Was vc est PD
      vc_pos_def <-  matrixcalc::is.positive.definite(var_sims$Sigma[-1, -1])

      if(vc_pos_def == FALSE){
        warning("VC matrix not positive definite")
      }

      # b vcov
      b_vcov <- cov(var_sims$sims[, b_index])

      # Beta tilde vcov
      beta_tilde_vcov <- cov(var_sims$sims[, beta_tilde_index])
      est_vc <- var_sims$Sigma
      var_sims <- var_sims$sims
      boot <- FALSE

    }else{

      var_sims <- varianceBoot(Y = Y, X = X, S = S, nsims = nsims_var, N = N, int_vars = int_vars,
                               sq_vars = sq_vars, int1 = int1, int2 = int2, index_sq = index_sq)

      # b vcov
      b_vcov <- cov(var_sims[, b_index])

      # Beta tilde vcov
      beta_tilde_vcov <- cov(var_sims[, beta_tilde_index])
    }

  }else{

    if(length(int_vars) > 0){
      # Interaction model
      int_names <- unlist(strsplit(int_vars, split = ":"))
      int1 <- which(colnames(X) == int_names[1])
      int2 <- which(colnames(X) == int_names[2])

      # Construct S matrix
      S <- matrix(0, nrow = ncol(X), ncol = ncol(X))

      # Find position of relevant columns
      err_ind <- match(colnames(X)[-c(1, ncol(X))], colnames(data))
      diag(S) <- c(0, S_vec[err_ind], 0)

      # Bias correct
      bias_correct <- betaTildeInt(coef = b, S = S, X = X, N = N, int1 = int1, int2 = int2)
      beta_tilde <- bias_correct$bias_correct
      x_prime_x <- bias_correct$x_prime_x

      # Estimate sigma^2
      sigma_sq <- sigmaSq(Y = Y, bias_correct = beta_tilde, S = S, X = X)

      # Variance
      b_index <- 1:length(b)
      beta_tilde_index <- (length(b) + 1):(2*length(b))
      var_sims <- varianceBoot(Y = Y, X = X, S = S, nsims = nsims_var, N = N, int_vars = int_vars,
                               sq_vars = sq_vars, int1 = int1, int2 = int2, index_sq = index_sq)
      b_vcov <- cov(var_sims[, b_index])
      beta_tilde_vcov <- cov(var_sims[, beta_tilde_index])

    } else{
      # Squared model
      sq_var <- gsub("I\\(", "", sq_vars)
      sq_var <- gsub("\\^2\\)", "", sq_var)
      index_sq <- which(colnames(X) == sq_var)

      # Construct S matrix
      S <- matrix(0, nrow = ncol(X), ncol = ncol(X))

      # Find position of relevant columns
      err_ind <- match(colnames(X)[-c(1, ncol(X))], colnames(data))
      diag(S) <- c(0, S_vec[err_ind], 0)

      # Bias correct
      bias_correct <- betaTildeSq(coef = b, S = S, X = X, N = N, index_sq = index_sq)
      beta_tilde <- bias_correct$bias_correct
      x_prime_x <- bias_correct$x_prime_x

      # Estimate sigma^2
      sigma_sq <- sigmaSq(Y = Y, bias_correct = beta_tilde, S = S, X = X)

      # Variance
      b_index <- 1:length(b)
      beta_tilde_index <- (length(b) + 1):(2*length(b))
      var_sims <- varianceBoot(Y = Y, X = X, S = S, nsims = nsims_var, N = N, int_vars = int_vars,
                               sq_vars = sq_vars, int1 = int1, int2 = int2, index_sq = index_sq)
      b_vcov <- cov(var_sims[, b_index])
      beta_tilde_vcov <- cov(var_sims[, beta_tilde_index])

    }
  }

  output <- list(
    b = b,
    b_vcov = b_vcov,
    beta_tilde = beta_tilde,
    beta_tilde_vcov = beta_tilde_vcov,
    var_sims = var_sims,
    Sigma_sq_hat = sigma_sq,
    vc_pos_def = vc_pos_def,
    boot = boot,
    est_vc = est_vc,
    Y = c(S_vec[match(colnames(model.frame(reg))[1], colnames(data))], Y),
    X = rbind(c(0, S_vec[match(colnames(X)[-1], colnames(data))]), X),
    S = S,
    formula = formula
  )

  class(output) <- "lmdp"

  return(invisible(output))
  }


#' Summary lmdp
#'
#' This function takes an lmdp object as input and produces a summary table comparable to
#' a summary of an \code{\link{lm}} object.
#'
#' @param lmdp_object Output from \code{\link{lmdp}}
#' @return Table of consistent coefficient estimates and corrresponding standard errors, t-stats
#' and p-values.
#' @export
summary.lmdp <- function(lmdp_object)
  {
  beta_tilde <- lmdp_object$beta_tilde
  beta_tilde_se <- sqrt(diag(lmdp_object$beta_tilde_vcov))
  t_val <- beta_tilde/beta_tilde_se
  p_val <- 2*(1 - pnorm(abs(t_val)))
  mod_output <- cbind(beta_tilde, beta_tilde_se, t_val, p_val)
  rownames(mod_output) <- names(beta_tilde)
  colnames(mod_output) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
  return(round(mod_output, 4))
  }
