\name{da4}
\alias{da4}
\docType{data}
\title{Clustered partly interva-censored data}
\description{A simulated clsutered partly interval-censored data set based on 
PH model with random intercept and random effect for x2:

	\code{lambda(t|x)=lambda_{0}(t)exp(x1+x2+phi+phi_trt*x2)}.
}
\usage{data(da4)}
\format{
  \tabular{ll}{
    L:\tab Left endpoints of observed time intervals. \cr
    R: \tab Right endpoints of observed time intervals. \cr
    y: \tab Censoring indicator: 0=left-censored, 1=interval-censored, 
	2=right-censored, 3=exact. \cr
    X1:\tab Covariate 1. \cr
    X2:\tab Covariate 2. \cr
    area:\tab Cluster ID. \cr
    IC:\tab General interval-censored indicator: 
	1=general interval-censored, 0=exact. \cr
    ID:\tab Subject ID. \cr
  }
}