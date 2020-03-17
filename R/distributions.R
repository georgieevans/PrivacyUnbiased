#' Distribution plots
#'
#'
#' @param X dp vector
#' @param Z_est Estimate of private vector
#' @return Histogram of dp vector and estimated private vector
#' @export
#'
plotDist <- function(X, Z_est, plot_dp){
  if(plot_dp == TRUE){
    plot <- ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = X, fill = 'X'), alpha = .7, bins = length(X)/2000)  +
      ggplot2::geom_histogram(ggplot2::aes(x = Z_est, fill = 'Est Z'), alpha = .7, bins = length(X)/2000) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = c(0.9, 0.9),
            legend.justification =  c(0.9, 0.9),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(x = '', title = 'Estimated vs. observed data', fill = '') +
      ggplot2::scale_fill_manual(values = c('navyblue', 'orange'))
  }else{
    plot <- ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = Z_est, fill = 'Est Z'), alpha = .7, bins = length(X)/2000) +
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

#' Normal fit
#'
#'  Estimate via method of moments fit and paramaters of a normal distribution from DP data
#' @param X dp vector
#' @param S dp noise added to obtain X
#' @param R numbers of implied moments to fit
#' @param moments_df estimated moments of Z
#' @param plot_dp TRUE/FALSE value denoting whether to generate plot
#' @return
#' Returns list of estimates:
#' \item{params}{Method of moment estimates of mean and variance of normal distribution}
#' \item{moment_df}{Moments used to estimate params}
#' \item{moment_fit}{Ratio of estimated moments to implied moments from paramaterized distribution}
#' \item{moment_precision}{Ratio of moment mean to moment standard error}
#' \item{plot}{Plot of fitted distribution}
#' @export
#'
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

#' Poisson fit
#'
#'  Estimate via method of moments fit and paramaters of a Poisson distribution from DP data
#' @param X dp vector
#' @param S dp noise added to obtain X
#' @param R numbers of implied moments to fit
#' @param moments_df estimated moments of Z
#' @param plot_dp TRUE/FALSE value denoting whether to generate plot
#' @return
#' Returns list of estimates:
#' \item{params}{Method of moment estimate of mean/variance of lambda distribution}
#' \item{moment_df}{Moments used to estimate params}
#' \item{moment_fit}{Ratio of estimated moments to implied moments from paramaterized distribution}
#' \item{moment_precision}{Ratio of moment mean to moment standard error}
#' \item{plot}{Plot of fitted distribution}
#' @export
#'
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

#' NegBin fit
#'
#'  Estimate via method of moments fit and paramaters of a negative binomial distribution from DP data
#' @param X dp vector
#' @param S dp noise added to obtain X
#' @param R numbers of implied moments to fit
#' @param moments_df estimated moments of Z
#' @param plot_dp TRUE/FALSE value denoting whether to generate plot
#' @return
#' Returns list of estimates:
#' \item{params}{Method of moment estimate of paramaters of a NB distribution}
#' \item{moment_df}{Moments used to estimate params}
#' \item{moment_fit}{Ratio of estimated moments to implied moments from paramaterized distribution}
#' \item{moment_precision}{Ratio of moment mean to moment standard error}
#' \item{plot}{Plot of fitted distribution}
#' @export
#'
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

#' Zero Inflated Poisson paramaters
#'
#' Estimates paramaters of a ZIP using moments
#' @param moments vector of moments
#' @return Vector of ZIP paramaters
#' @export
#'
estZIP <- function(moments){
  pi <- 1 - moments[2]/moments[1]
  lambda <- (moments[2] - moments[1])/moments[1]
  return(c(pi = pi, lambda = lambda))
}

#' Zero-Inflated Poisson fit
#'
#'  Estimate via method of moments fit and paramaters of a ZIP distribution from DP data
#' @param X dp vector
#' @param S dp noise added to obtain X
#' @param R numbers of implied moments to fit
#' @param moments_df estimated moments of Z
#' @param plot_dp TRUE/FALSE value denoting whether to generate plot
#' @return
#' Returns list of estimates:
#' \item{params}{Method of moment estimate of mean/variance of ZIP distribution}
#' \item{moment_df}{Moments used to estimate params}
#' \item{moment_fit}{Ratio of estimated moments to implied moments from paramaterized distribution}
#' \item{moment_precision}{Ratio of moment mean to moment standard error}
#' \item{plot}{Plot of fitted distribution}
#' @export
#'
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


#' Zero-Inflated NegBin paramaters
#'
#' Estimates paramaters of a ZINB using estimated moments
#' @param moments vector of moments
#' @return Vector of ZNB paramaters
#' @export
#'
estZNB <- function(moments){
  pi <-  (moments[1]^2*moments[2] + moments[1]*(moments[2] + moments[3]) -
            2*moments[2]^2 - moments[1]^3)/(moments[1]*(moments[2] + moments[3]) - 2*moments[2]^2)
  p <- (moments[1]*(moments[2] - moments[3]) + moments[2]^2 - moments[1]^2)/
    (moments[2]^2 - moments[1]*moments[3])
  r <- (2*moments[2]^2 - moments[1]*(moments[2] + moments[3]))/
    (moments[1]^2 + moments[1]*(moments[3] - moments[2]) - moments[2]^2)
  return(c(pi = pi, p = p, r = r))
}

#' Zero-Inflated NegBin fit
#'
#'  Estimate via method of moments fit and paramaters of a ZINB distribution from DP data
#' @param X dp vector
#' @param S dp noise added to obtain X
#' @param R numbers of implied moments to fit
#' @param moments_df estimated moments of Z
#' @param plot_dp TRUE/FALSE value denoting whether to generate plot
#' @return
#' Returns list of estimates:
#' \item{params}{Method of moment estimate of mean/variance of ZINB distribution}
#' \item{moment_df}{Moments used to estimate params}
#' \item{moment_fit}{Ratio of estimated moments to implied moments from paramaterized distribution}
#' \item{moment_precision}{Ratio of moment mean to moment standard error}
#' \item{plot}{Plot of fitted distribution}
#' @export
#'
paramsZINB <- function(X, S, R, moments_df, plot_dp){
  moments <- moments_df[,1]
  est <- as.numeric(estZNB(moments))
  if(est[2] < 0| est[2] > 1| est[1] < 0| est[1] > 1){
    warning('ZINB estimates outside logical bounds')
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


#' Paramaterize distrributions
#'
#'  Paramaterize distributions for private data based on moments estimated from DP data
#' @param variable Variable name in dataa
#' @param data DP data (first row should contain sd of noise added to `variable``)
#' @param distributions Vector of distributions to fit. Options include 'Normal', 'Poisson', 'NB', 'ZIP', 'ZINB'
#' @param moments_fit Number of moments to fit implied distribution to
#' @param Plot Plot histogram of distributions TRUE/FALSE
#' @param plot_dp Include DP in plot
#' @return
#' Returns list of paramater estimates:
#' @export
#'
distributionDP <- function(variable, data, distributions, moments_fit = 6, plot = TRUE, plot_dp = TRUE){

  var <- eval(substitute(variable),data, parent.frame())
  S <- var[1]
  X <- var[-1]
  n <- length(X)
  moments_df <- Rmoments(X, S, R, n = n)
  moments <- moments_df[,1]

  output_list <- list()
  params_normal <- NULL
  params_poisson <- NULL
  params_ZIP <- NULL
  params_NB <- NULL
  params_ZINB <- NULL

  if(!is.na(match('Normal', distributions))){
    params_normal <- paramsNormal(X = X, S = S, R = moments_fit, moments_df = moments_df, plot_dp = plot_dp)
    if(!is.null(params_normal$moment_fit)){
      output_list[[length(output_list) + 1]] <- list(dist = 'Normal', fit = params_normal$moment_fit)
    }
  }

  if(!is.na(match('Poisson', distributions))){
    params_poisson <- paramsPoisson(X = X, S = S, R = moments_fit, moments_df = moments_df, plot_dp = plot_dp)
    if(!is.null(params_poisson$moment_fit)){
      output_list[[length(output_list) + 1]] <- list(dist = 'Poisson', fit = params_poisson$moment_fit)
    }
  }

  if(!is.na(match('ZIP', distributions))){
    params_ZIP <- paramsZIP(X = X, S = S, R = moments_fit, moments_df = moments_df, plot_dp = plot_dp)
    if(!is.null(params_ZIP$moment_fit)){
      output_list[[length(output_list) + 1]] <- list(dist = 'ZIP', fit = params_ZIP$moment_fit)
    }
  }


  if(!is.na(match('NB', distributions))){
    params_NB <- paramsNB(X = X, S = S, R = moments_fit, moments_df = moments_df, plot_dp = plot_dp)
    if(!is.null(params_NB$moment_fit)){
     output_list[[length(output_list) + 1]] <- list(dist = 'NB', fit = params_NB$moment_fit)
    }
  }

  if(!is.na(match('ZINB', distributions))){
    params_ZINB <- paramsZINB(X = X, S = S, R = moments_fit, moments_df = moments_df, plot_dp = plot_dp)
    if(!is.null(params_ZINB$moment_fit)){
      output_list[[length(output_list) + 1]] <- list(dist = 'ZINB', fit = params_ZINB$moment_fit)
    }
  }

  fit_table <- as.data.frame(do.call(rbind, lapply(1:length(output_list), FUN = function(i) output_list[[i]]$fit)))
  rownames(fit_table) <- sapply(1:length(output_list), FUN = function(i) output_list[[i]]$dist)
  colnames(fit_table) <- paste0('mu', c(1:moments_fit))

  cat("\n Estimated Raw Moments/Implied Distribution Moments:\n\n")
  print(fit_table)

  output <- list(
    Normal = params_normal,
    Poisson = params_poisson,
    ZIP = params_ZIP,
    NB = params_NB,
    ZINB = params_ZINB
  )

  return(invisible(output))
}
