% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Shape-class.R
\docType{methods}
\name{setErrors<-,Shape-method}
\alias{setErrors<-,Shape-method}
\title{Set Errors for \linkS4class{Shape}}
\usage{
\S4method{setErrors}{Shape}(object) <- value
}
\arguments{
\item{object}{Object of class \linkS4class{Shape}.}

\item{value}{Matrix of length(shapeParamsUsed[[1]]) x fpLen containing beta errors for shape coefficients in -ddG/angstrom units.}
}
\description{
Sets Shape.errors slot from object of class \linkS4class{Shape}.
}
