\name{clusterPIC_Z}
\alias{clusterPIC_Z}
\title{Mixed effects PH model for clustered partly interval-censored data}
\description{Fit a Bayesian semiparametric mixed effects PH model for clustered partly 
interval-censored data with random effects for one or more predictors. 
Each random effect follows a normal distribution \code{N(0, tau^{-1})}.}
\usage{
clusterPIC_Z(L, R, y, xcov, IC, scale.designX, scaled, zcov, area, binary, I, 
order, knots, grids, a_eta, b_eta, a_ga, b_ga, a_tau, b_tau, beta_iter, 
phi_iter, beta_cand, phi_cand, beta_sig0, x_user, 
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
  \item{zcov}{The design matrix for the q random effects.}
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
The mixed effects PH model is:

\code{h(t_{ij}|x_{ij},z_{i})=h_{0}(t_{ij})exp(beta'x_{ij}+phi_{i}'z_{i})},

for the jth subject in the ith cluster.

Each of the q random effects is sampled using MH algorithm separately.
}

\value{
a list containing the following elements:
 \item{N}{The sample size.}
 \item{parbeta}{A \code{total} by \code{p} matrix of MCMC draws of \code{beta_r}, r=1, ..., p.}
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
\author{Chun Pan}

\examples{
# Number of iterations set to very small for CRAN automatic testing
data(da4)
J=rep(1,nrow(da4))
zcov=cbind(J,da4[,4]) # The 4th column of da4 is x1.
try7<-PICBayes(formula=Surv(L,R,type='interval2')~x1+x2,data=data.frame(da4),
model='clusterPIC_Z',IC=da4[,7],scale.designX=TRUE,scaled=c(1,0),zcov=zcov,
area=da4[,6],binary=c(0,1),I=25,order=3,knots=c(0,2,6,max(da4[,1:2],na.rm=TRUE)+1),
grids=seq(0.1,10.1,by=0.1),a_eta=1,b_eta=1,a_ga=1,b_ga=1,a_tau=1,b_tau=1,
beta_iter=11,phi_iter=11,beta_cand=c(1,1),phi_cand=1,beta_sig0=10,
x_user=NULL,total=30,burnin=10,thin=1,conf.int=0.95,seed=1)
}















