% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Intercept-class.R
\docType{methods}
\name{setOldSig<-,Intercept-method}
\alias{setOldSig<-,Intercept-method}
\title{Set P-values for Previous Iterations for \linkS4class{Intercept}}
\usage{
\S4method{setOldSig}{Intercept}(object) <- value
}
\arguments{
\item{object}{Object of class \linkS4class{Intercept}.}

\item{value}{4-D array with dimension of 2 x (numViews) x (number of rounds) x (number of glm iterations) containing -ddG p-values for 'Intercept' offsets.}
}
\description{
Sets I.oldSig slot for object of class \linkS4class{Intercept}.
}
