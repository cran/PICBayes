\name{clusterIC_trt}
\alias{clusterIC_trt}
\title{PH model with random intercept and random treatment for clustered general interval-censored data}
\description{Fit a Bayesian semiparametric PH model with random intercept 
and random treatment for clustered general interval-censored data. 
Each random effect follows a normal distribution \code{N(0, tau^{-1})}.}
\usage{
clusterIC_trt(L, R, y, xcov, IC, scale.designX, scaled, xtrt, area, binary, I, 
order, knots, grids, a_eta, b_eta, a_ga, b_ga, a_tau, b_tau, a_tau_trt, 
b_tau_trt, beta_iter, phi_iter, beta_cand, phi_cand, beta_sig0, x_user, 
total, burnin, thin, conf.int, seed)
}
\arguments{
  \item{L}{The vector of left endpoints of the observed time intervals.}
  \item{R}{The vector of right endponts of the observed time intervals.}
  \item{y}{The vector of censoring indicator: 0=left-censored, 1=interval-censored, 
	2=right-censored.}
  \item{xcov}{The covariate matrix for the p predictors.}
  \item{IC}{The vector of general interval-censored indicator: 
	1=general interval-censored, 0=exact.}
  \item{scale.designX}{The TRUE or FALSE indicator of whether or not to scale the design matrix X.}
  \item{scaled}{The vector indicating whether each covariate is to be scaled: 1=to be scaled, 0=not.}
  \item{xtrt}{The covariate that has a random effect.}
  \item{area}{The vector of cluster ID.}
  \item{binary}{The vector indicating whether each covariate is binary.}
  \item{I}{The number of clusters.}
  \item{order}{The degree of basis I-splines: 1=linear, 2=quadratic, 3=cubic, etc.}
  \item{knots}{A sequence of knots to define the basis I-splines.}
  \item{grids}{A sequence of points at which baseline survival function is to be estimated.}
  \item{a_eta}{The shape parameter of Gamma prior for \code{gamma_l}.}
  \item{b_eta}{The rate parameter of Gamma prior for \code{gamma_l}.}
  \item{a_ga}{The shape parameter of Gamma prior for \code{e^{beta_r}}.}
  \item{b_ga}{The rate parameter of Gamma prior for \code{e^{beta_r}}.}
  \item{a_tau}{The shape parameter of Gamma prior for random intercept precision \code{tau}.}
  \item{b_tau}{The rate parameter of Gamma prior for random intercept precision \code{tau}.}
  \item{a_tau_trt}{The shape parameter of Gamma prior for random treatment precision \code{tau_trt}.}
  \item{b_tau_trt}{The rate parameter of Gamma prior for random treatment precision \code{tau_trt}.}
  \item{beta_iter}{The number of initial iterations in the Metropolis-Hastings sampling for \code{beta_r}.}
  \item{phi_iter}{The number of initial iterations in the Metropolis-Hastings sampling for \code{phi_i}.}
  \item{beta_cand}{The sd of the proposal normal distribution in the initial MH sampling for \code{beta_r}.}
  \item{phi_cand}{The sd of the proposal normal distribution in the initial MH sampling for \code{phi_i}.}
  \item{beta_sig0}{The sd of the prior normal distribution for \code{beta_r}.}
  \item{x_user}{The user-specified covariate vector at which to estimate survival function(s).}
  \item{total}{The number of total iterations.}
  \item{burnin}{The number of burnin.}
  \item{thin}{The frequency of thinning.}
  \item{conf.int}{The confidence level of the CI for \code{beta_r}.}
  \item{seed}{A user-specified random seed.}
}
\details{
The baseline cumulative hazard is approximated by a linear combination of I-splines:

                    \code{sum_{l=1}^{K}(gamma_l*b_l(t))}.

For a binary prdictor, we sample \code{e^{beta_r}}, with Gamma prior.

The regression coefficient \code{beta_r} for a continuous predictor, random intercept \code{phi_i}, 
and random treatment \code{phi_trt_i} are sampled using MH algorithm. 
During the initial \code{beta_iter} iterations, sd of the proposal distribution is \code{beta_cand}. 
Afterwards, proposal sd is set to be the sd of available MCMC draws. 
Same method for \code{phi_i} and \code{phi_trt_i}. 
}

\value{
a list containing the following elements:
 \item{N}{The sample size.}
 \item{parbeta}{A \code{total} by \code{p} matrix of MCMC draws of \code{beta_r}, r=1, ..., p.}
 \item{parsurv0}{A \code{total} by \code{length(grids)} matrix, each row contains the baseline survival at \code{grids} from one iteration.}
 \item{parsurv}{A \code{total} by \code{length(grids)*G} matrix, each row contains the survival at \code{grids} from one iteration. 
	G is the number of sets of user-specified covariate values.}
 \item{parphi}{A \code{total} by \code{I} matrix of MCMC draws of \code{phi_i}, i=1,...,I.}
 \item{parphi_trt}{A \code{total} by \code{I} matrix of MCMC draws of \code{phi_trt_i}, i=1,...,I.}
 \item{partau}{A \code{total} by 1 vector of MCMC draws of \code{tau}.}
 \item{partau_trt}{A \code{total} by 1 vector of MCMC draws of \code{tau_trt}.}
 \item{coef}{A vector of regression coefficient estimates.}
 \item{coef_ssd}{A vector of sample standard deviations of regression coefficient estimates.}
 \item{coef_ci}{The credible intervals for the regression coefficients.}
 \item{S0_m}{The estimated baseline survival at \code{grids}.}
 \item{S_m}{The estimated survival at \code{grids} with user-specified covariate values \code{x_user}.}
 \item{grids}{The sequance of points where baseline survival function is estimated.}
 \item{DIC}{Deviance information criterion.}
 \item{NLLK}{Negative log pseudo-marginal likelihood.}
}
\author{Chun Pan}

