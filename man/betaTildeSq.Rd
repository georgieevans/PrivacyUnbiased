% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{betaTildeSq}
\alias{betaTildeSq}
\title{Bias Correct OLS model with squared variable}
\usage{
betaTildeSq(coef, S, X, N, index_sq)
}
\arguments{
\item{coef}{Vector of OLS coefficients}

\item{S}{Matrix of error terms. Dimensions must be length(coef) x length(coef)}

\item{X}{Design matrix that produced OLS coefficients}

\item{N}{Number of rows in X}

\item{index_sq}{Column index of X covariate that is squared}
}
\value{
Returns a list:
\item{bias_correct}{Vector of bias corrected coefficients}
\item{x_prime_x}{X'X matrix to be used in variance estimation}
}
\description{
This function takes OLS coefficients from a model including
a squared variable as inputs and bias corrects according to the noise matrix, S
}
