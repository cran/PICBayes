\name{mCRC}
\alias{mCRC}
\docType{data}
\title{Colorectal cancer data}
\description{A progression-free survival data set derived by the author from 
a phase 3 metastatic colorectal cancer clinical trial.}
\usage{data(mCRC)}
\format{
  \tabular{ll}{
    L:\tab Left endpoints of observed time intervals. \cr
    R: \tab Right endpoints of observed time intervals. \cr
    y: \tab Censoring indicator: 0=left-censored, 1=interval-censored, 
	2=right-censored, 3=exact. \cr
    TRT_C:\tab Treatment arm: 0 = FOLFIRI alone, 1 = Panitumumab + FOLFIRI. \cr
    KRAS_C:\tab Tumor KRAS mutation status: 0 = wild-type, 1 = mutant. \cr
    SITE:\tab Clinical site where a patient is treated. \cr
    IC:\tab General interval-censored indicator: 1=general interval-censored, 0=exact. \cr
    ID:\tab Subject ID. \cr
  }
}
