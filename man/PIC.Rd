\name{PIC}
\alias{PIC}
\title{PH model for partly interval-censored data}
\description{Fit a Bayesian semiparametric PH model to partly interval-censored data.}
\usage{
PIC(L, R, y, xcov, IC, scale.designX, scaled, binary, order, knots, grids, 
a_eta, b_eta, a_ga, b_ga, beta_iter, beta_cand, beta_sig0, x_user, 
total, burnin, thin, conf.int, seed)
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
  \item{binary}{The vector indicating whether each covariate is binary: 1=binary, 0=not.}
  \item{order}{The degree of basis I-splines: 1=linear, 2=quadratic, 3=cubic, etc.}
  \item{knots}{A sequence of knots to define the basis I-splines.}
  \item{grids}{A sequence of points at which baseline survival function is to be estimated.}
  \item{a_eta}{The shape parameter of Gamma prior for \code{gamma_l}.}
  \item{b_eta}{The rate parameter of Gamma prior for \code{gamma_l}.}
  \item{a_ga}{The shape parameter of Gamma prior for \code{e^{beta_r}}.}
  \item{b_ga}{The rate parameter of Gamma prior for \code{e^{beta_r}}.}
  \item{beta_iter}{The number of initial iterations in the Metropolis-Hastings sampling for \code{beta_r}.}
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

The baseline hazard is approximated by a linear combination of basis M-splines:

                    \code{sum_{l=1}^{K}(gamma_l*M_l(t))}.

For a binary prdictor, we sample \code{e^{beta_r}}, with Gamma prior.

The regression coefficient \code{beta_r} for a continuous predictor is sampled 
using MH algorithm. During the initial \code{beta_iter} iterations, sd of the 
proposal distribution is \code{beta_cand}. Afterwards, proposal sd is set to be 
the sd of available MCMC draws. 
}


\value{
a list containing the following elements:
 \item{N}{The sample size.}
 \item{parbeta}{A \code{total} by \code{p} matrix of MCMC draws of \code{beta_r}, r=1, ..., p. }
 \item{parsurv0}{A \code{total} by \code{length(grids)} matrix, each row contains the baseline survival at \code{grids} from one iteration.}
 \item{parsurv}{A \code{total} by \code{length(grids)*G} matrix, each row contains the survival at \code{grids} from one iteration. 
	G is the number of sets of user-specified covariate values.}
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
Pan, C., Cai, B., and Wang, L. (2020). A Bayesian approach for analyzing partly 
interval-censored data under the proportional hazards model.  
\emph{Statistical Methods in Medical Research}, 

DOI: 10.1177/0962280220921552.
}
\author{Chun Pan}

\examples{
# Number of iterations set to very small for CRAN automatic testing
data(da1)
try1<-PICBayes(formula=Surv(L,R,type='interval2')~x1+x2,data=data.frame(da1),
model='PIC',IC=da1[,6],scale.designX=TRUE,scale=c(1,0),binary=c(0,1),
order=3,knots=c(0,2,6,max(da1[,1:2],na.rm=TRUE)+1),grids=seq(0.1,10.1,by=0.1),
a_eta=1,b_eta=1,a_ga=1,b_ga=1,beta_iter=11,beta_cand=1,beta_sig0=10,
x_user=NULL,total=60,burnin=10,thin=1,conf.int=0.95,seed=1)
}























