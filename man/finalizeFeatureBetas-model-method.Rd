% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-class.R
\docType{methods}
\name{finalizeFeatureBetas,model-method}
\alias{finalizeFeatureBetas,model-method}
\title{Finalize Beta Values}
\usage{
\S4method{finalizeFeatureBetas}{model}(object)
}
\arguments{
\item{object}{Object of class \linkS4class{model}.}
}
\description{
Finalizes beta values in \linkS4class{model} object by restoring NA's that were removed to allow for the iteration.
}
