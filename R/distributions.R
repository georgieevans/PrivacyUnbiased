plotDist <- function(X, Z_est, plot_dp){
  if(plot_dp == TRUE){
    plot <- ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = X, fill = 'X'), alpha = .7, bins = 35)  +
      ggplot2::geom_histogram(ggplot2::aes(x = Z_est, fill = 'Est Z'), alpha = .7, bins = 40) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = c(0.9, 0.9),
            legend.justification =  c(0.9, 0.9),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(x = '', title = 'Estimated vs. observed data', fill = '') +
      ggplot2::scale_fill_manual(values = c('navyblue', 'orange'))
  }else{
    plot <- ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = Z_est, fill = 'Est Z'), alpha = .7, bins = 40) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = c(0.9, 0.9),
                     legend.justification =  c(0.9, 0.9),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(x = '', title = 'Estimated data', fill = '') +
      ggplot2::scale_fill_manual(values = c('navyblue'))
  }
  return(plot)

}

paramsNormal <- function(X, S, R, moments_df, plot_dp){
  moments <- moments_df[, 1]
  mu <- mean(X)
  var <- var(X) - S^2
  Z_est <- rnorm(n, mean = mu, sd = sqrt(var))
  implied_moments <- sapply(1:nrow(moments_df), FUN = function(r) mean(Z_est^r))
  return(list(params = c(mu, var),
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2),
              plot = plotDist(X, Z_est, plot_dp)))
}


paramsPoisson <- function(X, S, R, moments_df, plot_dp){
  moments <- moments_df[, 1]
  lambda <- mean(X)
  Z_est <- rpois(n, lambda)
  implied_moments <- sapply(1:nrow(moments_df), FUN = function(r) mean(Z_est^r))

  return(list(params = c(lambda),
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2),
              plot = plotDist(X, Z_est, plot_dp)))

}

paramsNB <- function(X, S, R, moments_df, plot_dp){
  moments <- moments_df[, 1]
  var <- var(X) - S^2
  p <- 1 - (mean(X)/var)
  r <- (-mean(X)^2)/(mean(X) - var)
  Z_est <- rnbinom(n, size = r, mu = p*r/(1-p))
  implied_moments <- sapply(1:nrow(moments_df), FUN = function(r) mean(Z_est^r))
  return(list(params = c(p, r),
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2),
              plot = plotDist(X, Z_est, plot_dp)))

}

estZIP <- function(moments){
  pi <- 1 - moments[2]/moments[1]
  lambda <- (moments[2] - moments[1])/moments[1]
  return(c(pi = pi, lambda = lambda))
}


paramsZIP <- function(X, S, R, moments_df, plot_dp){
  moments <- moments_df[,1]
  est <- as.numeric(estZIP(moments))
  if(est[1] < 0| est[1] > 1){
    warning('ZIP estimates outside logical bounds')
    return(list(params = est))
  }
  ind <- sample(c(0, 1), n, prob = c(est[1], 1 - est[1]), replace = TRUE)
  Z_est <- ind*rpois(n, est[2])
  implied_moments <- sapply(1:length(moments), FUN = function(r) mean(Z_est^r))

  return(list(params = est,
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2),
              plot = plotDist(X, Z_est, plot_dp)))
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

paramsZINB <- function(X, S, R, moments_df, plot_dp){
  moments <- moments_df[,1]
  est <- as.numeric(estZNB(moments))
  if(est[2] < 0| est[2] > 1| est[1] < 0| est[1] > 1){
    warning('Estimates outside logical bounds')
    return(list(params = est))
  }
  ind <- sample(c(0, 1), 2*n, prob = c(est[1], 1 - est[1]), replace = TRUE)
  Z_est <- ind*rnbinom(2*n, mu = est[2]*est[3]/(1-est[2]), size = est[3])
  implied_moments <- sapply(1:length(moments), FUN = function(r) mean(Z_est^r))

  return(list(params = est,
              moment_df = moments_df,
              moment_fit = round(moments/implied_moments, 2),
              moment_precision = round(moments_df[, 1]/moments_df[, 2], 2),
              plot = plotDist(X, Z_est, plot_dp)))
}


distributionDP <- function(variable, data, distributions, moments_fit = 6, plot = TRUE, plot_dp = TRUE){

  var <- eval(substitute(variable),data, parent.frame())
  S <- var[1]
  X <- var[-1]
  n <- length(X)
  moments_df <- Rmoments(X, S, R, n = n)
  moments <- moments_df[,1]

  if(!is.na(match('Normal', distributions))){
    params_normal <- paramsNormal(X = X, S = S, R = moments_fit, plot_dp = plot_dp)
  }

  if(!is.na(match('Poisson', distributions))){
    params_poisson <- paramsPoisson(X = X, S = S, R = moments_fit, plot_dp = plot_dp)
  }

  if(!is.na(match('ZIP', distributions))){
    params_ZIP <- paramsZIP(X = X, S = S, R = moments_fit, plot_dp = plot_dp)
  }


  if(!is.na(match('NB', distributions))){
    params_NB <- paramsNB(X = X, S = S, R = moments_fit, plot_dp = plot_dp)
  }

  if(!is.na(match('ZNB', distributions))){
    params_ZINB <- paramsZINB(X = X, S = S, R = moments_fit, plot_dp = plot_dp)
  }

}
