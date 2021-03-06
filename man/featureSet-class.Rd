% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featureSet-class.R
\docType{class}
\name{featureSet-class}
\alias{featureSet-class}
\alias{featureSet}
\title{Class 'featureSet' Definition}
\description{
Defines an S4 class to represent features,
}
\section{Slots}{

\describe{
\item{\code{seedLen}}{number of base pairs in the seeding model.}

\item{\code{upFootprintExtend}}{maximum number of upstream positions to be fit beyond the footprint for any feature.}

\item{\code{downFootprintExtend}}{maximum number of downstream positions to be fit beyond the footprint for any feature.}

\item{\code{numViews}}{number of views scored on each DNA strand.}

\item{\code{rounds}}{rounds of SELEX data used to evaluate 'model' to which 'Intercept' belongs.}

\item{\code{rcSymmetric}}{logical: indicating whether reverse complement symmetric model fit should be used.}

\item{\code{shapeParamsUsed}}{list of shape features to be included in 'model'.}

\item{\code{N}}{object of class 'N'}

\item{\code{Intercept}}{object of class 'Intercept'}

\item{\code{Shape}}{object of class 'Shape'}

\item{\code{glmFits}}{List of glm fit summaries.}

\item{\code{designMatrixSummary}}{List of design matrices from getDesignMatrix corresponding to the glmFits.}
}}

