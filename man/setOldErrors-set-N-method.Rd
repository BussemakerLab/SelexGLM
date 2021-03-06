% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/N-class.R
\docType{methods}
\name{setOldErrors<-,N-method}
\alias{setOldErrors<-,N-method}
\title{Set Errors for Previous Iterations for \linkS4class{N}}
\usage{
\S4method{setOldErrors}{N}(object) <- value
}
\arguments{
\item{object}{Object of class \linkS4class{N}.}

\item{value}{Array of 4 x fpLen x (iterations-1) containing beta errors for mono-nucleotides in -ddG units}
}
\description{
Sets N.oldErrors slot from object of class \linkS4class{N}.
}
