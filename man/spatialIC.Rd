\name{spatialIC}
\alias{spatialIC}
\title{PH model for spatial general interval-censored data}
\description{Fit a Bayesian semiparametric PH model with spatial frailty for 
spatially dependent general interval-censored data.}
\usage{
spatialIC(L, R, y, xcov, IC, scale.designX, scaled, area, binary, I, C, nn, 
order, knots, grids, a_eta, b_eta, a_ga, b_ga, a_lamb, b_lamb, beta_iter, 
phi_iter, beta_cand, beta_sig0, x_user, total, burnin, thin, conf.int, seed)
}
\arguments{
  \item{L}{The vector of left endpoints of the observed time intervals.}
  \item{R}{The vector of right endponts of the observed time intervals.}
  \item{y}{The vector of censoring indicator: 0=left-censored, 1=interval-censored, 
	2=right-censored, 3=exact.}
  \item{xcov}{The covariate matrix for the p predictors.}
  \item{IC}{The vector of general interval-censored indicator: 
	1=general interval-censored, 0=exact.}
  \item{scale.designX}{The TRUE or FALSE indicator of whether or not to scale the design matrix X.}
  \item{scaled}{The vector indicating whether each covariate is to be scaled: 1=to be scaled, 0=not.}
  \item{area}{The vector of area ID.}
  \item{I}{The number of areas.}
  \item{C}{The adjacency matrix.}
  \item{nn}{The vector of number of neighbors for each area.}
  \item{binary}{The vector indicating whether each covariate is binary.}
  \item{order}{The degree of basis I-splines: 1=linear, 2=quadratic, 3=cubic, etc.}
  \item{knots}{A sequence of knots to define the basis I-splines.}
  \item{grids}{A sequence of points at which baseline survival function is to be estimated.}
  \item{a_eta}{The shape parameter of Gamma prior for \code{gamma_l}.}
  \item{b_eta}{The rate parameter of Gamma prior for \code{gamma_l}.}
  \item{a_ga}{The shape parameter of Gamma prior for \code{e^{beta_r}}.}
  \item{b_ga}{The rate parameter of Gamma prior for \code{e^{beta_r}}.}
  \item{a_lamb}{The shape parameter of Gamma prior for spatial precision \code{lambda}.}
  \item{b_lamb}{The rate parameter of Gamma prior for spatial precision \code{lambda}.}
  \item{beta_iter}{The number of initial iterations in the Metropolis-Hastings sampling for \code{beta_r}.}
  \item{phi_iter}{The number of initial iterations in the Metropolis-Hastings sampling for \code{phi_i}.}
  \item{beta_cand}{The sd of the proposal normal distribution in the MH sampling for \code{beta_r}.}
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

The regression coefficient \code{beta_r} for a continuous predictor is sampled 
using MH algorithm. During the initial \code{beta_iter} iterations, sd of the 
proposal distribution is \code{beta_cand}. Afterwards, proposal sd is set to be 
the sd of available MCMC draws. 
}

\value{
a list containing the following elements:
 \item{N}{The sample size.}
 \item{parbeta}{A \code{total} by \code{p} matrix of MCMC draws of \code{beta_r}, r=1, ..., p.}
 \item{parsurv0}{A \code{total} by \code{length(grids)} matrix, each row contains the baseline survival at \code{grids} from one iteration.}
 \item{parsurv}{A \code{total} by \code{length(grids)*G} matrix, each row contains the survival at \code{grids} from one iteration. 
	G is the number of sets of user-specified covariate values.}
 \item{parphi}{A \code{total} by \code{I} matrix of MCMC draws of \code{phi_i}, i=1,...,I.}
 \item{parlamb}{A \code{total} by 1 matrix of MCMC draws of \code{lambda}.}
 \item{coef}{A vector of regression coefficient estimates.}
 \item{coef_ssd}{A vector of sample standard deviations of regression coefficient estimates.}
 \item{coef_ci}{The credible intervals for the regression coefficients.}
 \item{S0_m}{The estimated baseline survival at \code{grids}.}
 \item{S_m}{The estimated survival at \code{grids} with user-specified covariate values \code{x_user}.}
 \item{grids}{The sequance of points where baseline survival functions is estimated.}
 \item{DIC}{Deviance information criterion.}
 \item{NLLK}{Negative log pseudo-marginal likelihood.}
}
\references{
Pan, C. and Cai, B. (2020). A Bayesian model for spatial partly interval-censored data.  
\emph{Communications in Statistics - Simulation and Computation}, 
DOI: 10.1080/03610918.2020.1839497.
}
\author{Chun Pan}
