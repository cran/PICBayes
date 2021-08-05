\name{PICBayes}
\alias{PICBayes}
\alias{PICBayes.formula}
\alias{PICBayes.default}
\title{Bayesian models for partly interval-censored data and general interval-censored data}
\description{Calls one of the 16 functions to fit the correspoinding model.}
\usage{
PICBayes(L, ...)

\method{PICBayes}{default}(L,R,y,xcov,IC,model,scale.designX,scaled,xtrt,zcov,
area,binary,I,C,nn,order=3,knots,grids,a_eta=1,b_eta=1,a_ga=1,b_ga=1,a_lamb=1,
b_lamb=1,a_tau=1,b_tau=1,a_tau_trt=1,b_tau_trt=1,a_alpha=1,b_alpha=1,H=5,
a_tau_star=1,b_tau_star=1,a_alpha_trt=1,b_alpha_trt=1,H_trt=5,
a_tau_trt_star=1,b_tau_trt_star=1,beta_iter=1001,phi_iter=1001,
beta_cand,phi_cand,beta_sig0=10,x_user=NULL,
total=6000,burnin=1000,thin=1,conf.int=0.95,seed=1,...)

\method{PICBayes}{formula}(formula, data, ...)
}

\arguments{
  \item{L}{The vector of left endpoints of the observed time intervals.}
  \item{R}{The vector of right endponts of the observed time intervals.}
  \item{y}{The vector of censoring indicator: 0=left-censored, 1=interval-censored, 
	2=right-censored, 3=exact.}
  \item{xcov}{The covariate matrix for the p predictors.}
  \item{IC}{The vector of general interval-censored indicator: 
	1=general interval-censored, 0=exact.}
  \item{model}{A character string specifying the type of model. See details.}
  \item{scale.designX}{The TRUE or FALSE indicator of whether or not to scale the design matrix X.}
  \item{scaled}{The vector indicating whether each covariate is to be scaled: 1=to be scaled, 0=not.}
  \item{xtrt}{The covariate that has a random effect.}
  \item{zcov}{The design matrix for the q random effects.}
  \item{area}{The vector of cluster ID.}
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
  \item{a_tau}{The shape parameter of Gamma prior for random intercept precision \code{tau}.}
  \item{b_tau}{The rate parameter of Gamma prior for random intercept precision \code{tau}.}
  \item{a_tau_trt}{The shape parameter of Gamma prior for random treatment precision \code{tau_trt}.}
  \item{b_tau_trt}{The rate parameter of Gamma prior for random treatment precision \code{tau_trt}.}
  \item{a_alpha}{The shape parameter of Gamma prior for \code{alpha}.}
  \item{b_alpha}{The rate parameter of Gamma prior for \code{alpha}.}
  \item{H}{The number of distinct components in DP mixture prior under blocked Gibbs sampler.}
  \item{a_tau_star}{The shape parameter of \code{G_0} in DP mixture prior.}
  \item{b_tau_star}{The rate parameter of \code{G_0} in DP mixture prior.}
  \item{a_alpha_trt}{The shape parameter of Gamma prior for \code{alpha_trt}.}
  \item{b_alpha_trt}{The rate parameter of Gamma prior for \code{alpha_trt}.}
  \item{H_trt}{The number of distinct components in DP mixture prior under blocked Gibbs sampler for random treatment.}
  \item{a_tau_trt_star}{The shape parameter of \code{G_0} in DP mixture prior for random treatment.}
  \item{b_tau_trt_star}{The rate parameter of \code{G_0} in DP mixture prior for random treatment.}
  \item{beta_iter}{The number of initial iterations in the Metropolis-Hastings sampling for \code{beta_r}.}
  \item{phi_iter}{The number of initial iterations in the Metropolis-Hastings sampling for \code{phi_i}.}
  \item{beta_cand}{The sd of the proposal normal distribution in the MH sampling for \code{beta_r}.}
  \item{phi_cand}{The sd of the proposal normal distribution in the initial MH sampling for \code{phi_i}.}
  \item{beta_sig0}{The sd of the prior normal distribution for \code{beta_r}.}
  \item{x_user}{The user-specified covariate vector at which to estimate survival function(s).}
  \item{total}{The number of total iterations.}
  \item{burnin}{The number of burnin.}
  \item{thin}{The frequency of thinning.}
  \item{conf.int}{The confidence level of the CI for \code{beta_r}.}
  \item{seed}{A user-specified random seed.}
  \item{formula}{A formula expression with the response returned by the Surv function in the survival package.}
  \item{data}{A data frame that contains the variables named in the formula argument.}
   \item{...}{Other arguments if any.}
}
\details{
Possible values are "PIC", "spatialPIC", "clusterPIC_int", "clusterPIC_int_DP", 
"clusterPIC_trt", "clusterPIC_trt_DP", "clusterPIC_Z", and "clusterPIC_Z_DP" 
for partly interval-censored data; and "IC", "spatialIC", "clusterIC_int", 
"clusterIC_int_DP", "clusterIC_trt", "clusterIC_trt_DP", "clusterIC_Z", and 
"clusterIC_Z_DP" for general interval-censored data.
}
\value{
An object of class \code{PICBayes}. Refere to each specific function for its specific values.
}
\author{Chun Pan}

















