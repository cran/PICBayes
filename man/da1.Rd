\name{da1}
\alias{da1}
\docType{data}
\title{Partly interva-censored data}
\description{A simulated partly interval-censored data set based on:

	\code{lambda(t|x)=lambda_{0}(t)exp(x1+x2)}.
}
\usage{data(da1)}
\format{
  \tabular{ll}{
    L:\tab Left endpoints of observed time intervals. \cr
    R: \tab Right endpoints of observed time intervals. \cr
    y: \tab Censoring indicator: 0=left-censored, 1=interval-censored, 
	2=right-censored, 3=exact. \cr
    X1:\tab Covariate 1. \cr
    X2:\tab Covariate 2. \cr
    IC:\tab General interval-censored indicator: 
	1=general interval-censored, 0=exact. \cr
    ID:\tab Subject ID. \cr
  }
}

