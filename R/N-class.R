#' @include N-validity.R
#' @include generic-functions-definitions.R
#' @include plottingTools.R
NULL

#' Class 'N' for Mononucleotide Features
#' 
#' Defines an S4 class to represent mononuclueotide features.
#' 
#' @slot seedLen number of base pairs in the seeding model.
#' @slot N.upFootprintExtend number of upstream nucleotides to fit beyond the seed.
#' @slot N.downFootprintExtend number of downstream nucleotides to fit beyond the seed.
#' @slot fS.upFootprintExtend maximum number of upstream positions to be fit beyond the footprint for any feature.
#' @slot fS.downFootprintExtend maximum number of downstream positions to be fit beyond the footprint for any feature. 
#' @slot fpLen footprint length for full set of features.
#' @slot N.set set of positions for which mononucleotides are fit where positions are given relative to the full feature set footprint.
#' @slot N.equivMat matrix encoding reverse complement symmetries between positions with -1 and identical beta values shared between positions with 1.  
#' @slot N.values beta values for positions in N.set.
#' @slot N.errors beta value errors for positions in N.set.
#' @slot N.z beta value z scores for positions in N.set.
#' @slot N.sig beta value significance for positions in N.set.
#' @slot N.oldValues previous iteration beta values for positions in N.set.
#' @slot N.oldErrors previous iteration beta values' errors for positions in N.set.
#' @slot N.oldZ previous iteration beta values' z scores for positions in N.set.
#' @slot N.oldSig previous iteration beta values' significance for positions in N.set.
#' @export N
N <- setClass("N",
              representation(
                seedLen = "numeric", 
                N.upFootprintExtend = "numeric", 
                N.downFootprintExtend = "numeric",
                fS.upFootprintExtend = "numeric",
                fS.downFootprintExtend = "numeric",
                fpLen = "numeric",
                N.set = "numeric", 
                N.equivMat = "matrix", 
                N.values = "matrix",
                N.errors = "matrix", 
                N.z = "matrix", 
                N.sig = "matrix",
                N.oldValues = "array", 
                N.oldErrors = "array", 
                N.oldZ = "array",
                N.oldSig = "array"
              ), 
              validity = N.validity
)


#' Constructor for S4 class 'N'
#' 
#' Initializes an S4 class to represent mononuclueotide features
#' 
#' @export N
setMethod("initialize",
          "N", 
          function(.Object,
                   seedLen, 
                   N.upFootprintExtend, 
                   N.downFootprintExtend,
                   fS.upFootprintExtend,
                   fS.downFootprintExtend,
                   fpLen,
                   N.set, 
                   N.equivMat, 
                   N.values,
                   N.errors, 
                   N.z, 
                   N.sig,
                   N.oldValues, 
                   N.oldErrors, 
                   N.oldZ,
                   N.oldSig,
                   rcSymmetric,
                   ...) {
            
            
            if (!missing(seedLen)) {
              .Object@seedLen <- seedLen
            } else {
              stop('seedLen must be specified')
            }
            
            
            if (!missing(fS.upFootprintExtend)) {
              .Object@fS.upFootprintExtend <- fS.upFootprintExtend
            } else {
              stop('fS.upFootprintExtend must be specified')
            }
            
            if (!missing(fS.downFootprintExtend)) {
              .Object@fS.downFootprintExtend <- fS.downFootprintExtend
            } else {
              .Object@fS.downFootprintExtend = .Object@fS.upFootprintExtend
            }
            
            if (!missing(fpLen)) {
              message("fpLen will be overwritten if inconsistent with fpLen = seedLen+fS.upFootprintExtend+fS.downFootprintExtend")
            } 
            .Object@fpLen = .Object@seedLen+.Object@fS.upFootprintExtend+.Object@fS.downFootprintExtend
            
            if (missing(rcSymmetric)) {rcSymmetric = FALSE}
            if ((!missing(N.downFootprintExtend)) & (missing(fS.downFootprintExtend))) {
              if (N.downFootprintExtend > fS.upFootprintExtend) {
                stop("Invalid combination of slot values.'fS.downFootprintExtend' (not specified) is assigned to 'fS.upFootprintExtend'. 'N.downFootprintExtend' must be <= 'fS.upFootprintExtend' = 'fS.downFootprintExtend'.")
              }
            }
            
            if ((!missing(N.set)) & (missing(fS.downFootprintExtend)) & (missing(N.downFootprintExtend))) {
              if (max(N.set) > 2*fS.upFootprintExtend+seedLen) {
                stop("Invalid combination of slot values. 'fS.downFootprintExtend' (not specified) is assigned to 'fS.upFootprintExtend'. max('N.set') must be <= 'fS.upFootprintExtend' + 'fS.downFootprintExtend' + 'seedLen'.")
              }
            }
            
            if (missing(N.upFootprintExtend)) {
              if (missing(N.downFootprintExtend)) {
                if (missing(N.set)) {
                  .Object@N.upFootprintExtend = .Object@fS.upFootprintExtend
                  .Object@N.downFootprintExtend = .Object@fS.downFootprintExtend
                  minSet = 1
                  maxSet = .Object@fpLen
                  .Object@N.set = c(minSet:maxSet)
                } else if (!((length(N.set) == 1) & (N.set[1] == 0))){
                  .Object@N.set = N.set
                  minSet = min(N.set)
                  maxSet = max(N.set)
                  if (minSet <= .Object@fS.upFootprintExtend) {
                    .Object@N.upFootprintExtend = .Object@fS.upFootprintExtend-(minSet-1)
                  } else {
                    .Object@N.upFootprintExtend = 0
                  }
                  if (maxSet > .Object@fS.upFootprintExtend+.Object@seedLen) {
                    .Object@N.downFootprintExtend = maxSet-(.Object@fS.upFootprintExtend+.Object@seedLen)
                  } else {
                    .Object@N.downFootprintExtend = 0
                  }
                } else {
                  .Object@N.set = N.set
                  .Object@N.upFootprintExtend = 0
                  .Object@N.downFootprintExtend = 0
                } 
              } else {
                .Object@N.downFootprintExtend = N.downFootprintExtend
                if (missing(N.set)) {
                  .Object@N.upFootprintExtend = .Object@fS.upFootprintExtend
                  .Object@N.set = c(1:(.Object@fS.upFootprintExtend+.Object@seedLen+.Object@N.downFootprintExtend))
                } else if (!((length(N.set) == 1) & (N.set[1] == 0))){
                  .Object@N.set = N.set
                  minSet = min(N.set)
                  if (minSet <= .Object@fS.upFootprintExtend) {
                    .Object@N.upFootprintExtend = .Object@fS.upFootprintExtend-(minSet-1)
                  } else {
                    .Object@N.upFootprintExtend = 0
                  }
                } else {
                  .Object@N.set = N.set
                  .Object@N.upFootprintExtend = 0
                }
              }
            } else {
              .Object@N.upFootprintExtend = N.upFootprintExtend
              if (missing(N.downFootprintExtend)) {
                if (missing(N.set)) {
                  .Object@N.downFootprintExtend = .Object@fS.downFootprintExtend
                  minSet = .Object@fS.upFootprintExtend-.Object@N.upFootprintExtend+1
                  maxSet = .Object@fS.upFootprintExtend+.Object@seedLen+.Object@N.downFootprintExtend
                  .Object@N.set = c(minSet:maxSet)
                } else if (!((length(N.set) == 1) & (N.set[1] == 0))) {
                  minSet = min(N.set)
                  maxSet = max(N.set)
                  if (maxSet > .Object@fS.upFootprintExtend+.Object@seedLen) {
                    .Object@N.downFootprintExtend = maxSet-(.Object@fS.upFootprintExtend+.Object@seedLen)
                  } else {
                    .Object@N.downFootprintExtend = 0
                  }
                  .Object@N.set = N.set
                } else {
                  .Object@N.set = N.set
                  .Object@N.downFootprintExtend = 0
                }  
              } else {
                .Object@N.downFootprintExtend = N.downFootprintExtend
                if (missing(N.set)) {
                  minSet = .Object@fS.upFootprintExtend-.Object@N.upFootprintExtend+1
                  maxSet = .Object@fS.upFootprintExtend+.Object@seedLen+.Object@N.downFootprintExtend
                  .Object@N.set = c(minSet:maxSet)
                } else {
                  .Object@N.set = N.set
                }
              }
            }
            if (rcSymmetric == TRUE) {
              if (!missing(N.equivMat)) {
                message('rcSymmetric == TRUE & N.equivMat specified. N.equivMat input will not be used.')
              } 
              .Object@N.equivMat = buildSymmetricEquivalenceMatrix(.Object)
            } else {
              if (!missing(N.equivMat)) {
                colnames(N.equivMat) = c(1:ncol(N.equivMat))
                rownames(N.equivMat) = c(1:nrow(N.equivMat))
                .Object@N.equivMat = N.equivMat
              } else {
                .Object@N.equivMat = buildNullEquivalenceMatrix(.Object)
              }
            }
            
            
            
            
            .Object@N.set = unique(.Object@N.set)
            .Object@N.set = .Object@N.set[order(.Object@N.set)]
            
            
            # initialize N.values slot
            if (missing(N.values)) {
              .Object@N.values = matrix(0, nrow = 4, ncol = .Object@fpLen)
              rownames(.Object@N.values) = c("N.A", "N.C", "N.G", "N.T")
              colnames(.Object@N.values) = 1:.Object@fpLen
            } else {
              .Object@N.values = N.values
              rownames(.Object@N.values) = c("N.A", "N.C", "N.G", "N.T")
              colnames(.Object@N.values) = 1:.Object@fpLen
            } 
            
            # initialize N.errors slot
            if (missing(N.errors)) {
              .Object@N.errors = matrix(0, nrow = 4, ncol = .Object@fpLen)
              rownames(.Object@N.errors) =  c("N.A", "N.C", "N.G", "N.T")
              colnames(.Object@N.errors) = 1:.Object@fpLen
            } else {
              .Object@N.errors = N.errors
              rownames(.Object@N.errors) =  c("N.A", "N.C", "N.G", "N.T")
              colnames(.Object@N.errors) = 1:.Object@fpLen
            } 
            
            
            
            # initialize N.z slot
            if (missing(N.z)) {
              .Object@N.z = matrix(0, nrow = 4, ncol = .Object@fpLen)
              rownames(.Object@N.z) =  c("N.A", "N.C", "N.G", "N.T")
              colnames(.Object@N.z) = 1:.Object@fpLen
            } else {
              .Object@N.z = N.z
              rownames(.Object@N.z) =  c("N.A", "N.C", "N.G", "N.T")
              colnames(.Object@N.z) = 1:.Object@fpLen
            } 
            
            # initialize N.sig slot
            if (missing(N.sig)) {
              .Object@N.sig = matrix(0, nrow = 4, ncol = .Object@fpLen)
              rownames(.Object@N.sig) =  c("N.A", "N.C", "N.G", "N.T")
              colnames(.Object@N.sig) = 1:.Object@fpLen
            } else  {
              .Object@N.sig = N.sig
              rownames(.Object@N.sig) =  c("N.A", "N.C", "N.G", "N.T")
              colnames(.Object@N.sig) = 1:.Object@fpLen
            } 
            
            # initialize N.oldValues slot
            if (missing(N.oldValues)) {
              .Object@N.oldValues = array(0,
                                          dim = c(4, .Object@fpLen, 0),
                                          dimnames = list(c("N.A", "N.C", "N.G", "N.T"),
                                                          c(1:.Object@fpLen),
                                                          NULL))
            } else  {
              if (dim(N.oldValues)[3] > 0) {
                dimnames(N.oldValues) = list(c("N.A", "N.C", "N.G", "N.T"),
                                             c(1:.Object@fpLen),
                                             c((dim(N.oldValues)[3]-1):0))
              } else {
                dimnames(N.oldValues) = list(c("N.A", "N.C", "N.G", "N.T"),
                                             c(1:.Object@fpLen),
                                             NULL)
              }

              .Object@N.oldValues = N.oldValues
            } 
            
            # initialize N.oldErrors slot
            if (missing(N.oldErrors)) {
              .Object@N.oldErrors = array(0,
                                          dim = c(4, .Object@fpLen, 0),
                                          dimnames = list(c("N.A", "N.C", "N.G", "N.T"),
                                                          c(1:.Object@fpLen),
                                                          NULL))
            } else {
              if (dim(N.oldErrors)[3] > 0) {
                dimnames(N.oldErrors) = list(c("N.A", "N.C", "N.G", "N.T"),
                                             c(1:.Object@fpLen),
                                             c((dim(N.oldErrors)[3]-1):0))
              } else {
                dimnames(N.oldErrors) = list(c("N.A", "N.C", "N.G", "N.T"),
                                             c(1:.Object@fpLen),
                                             NULL)
              }
              .Object@N.oldErrors = N.oldErrors
            } 
            
            # initialize N.oldZ slot
            if (missing(N.oldZ)) {
              .Object@N.oldZ = array(0,
                                     dim = c(4, .Object@fpLen, 0),
                                     dimnames = list(c("N.A", "N.C", "N.G", "N.T"),
                                                     c(1:.Object@fpLen),
                                                     NULL))
            } else {
              if (dim(N.oldZ)[3] > 0) {
                dimnames(N.oldZ) = list(c("N.A", "N.C", "N.G", "N.T"),
                                        c(1:.Object@fpLen),
                                        c((dim(N.oldZ)[3]-1):0))
              } else {
                dimnames(N.oldZ) = list(c("N.A", "N.C", "N.G", "N.T"),
                                             c(1:.Object@fpLen),
                                             NULL)
              }
              
              .Object@N.oldZ = N.oldZ
            } 
            
            # initialize N.oldSig slot
            if (missing(N.oldSig)) {
              .Object@N.oldSig = array(0,
                                       dim = c(4, .Object@fpLen, 0),
                                       dimnames = list(c("N.A", "N.C", "N.G", "N.T"),
                                                       c(1:.Object@fpLen),
                                                       NULL))
            } else {
              if (dim(N.oldSig)[3] > 0) {
                dimnames(N.oldSig) = list(c("N.A", "N.C", "N.G", "N.T"),
                                          c(1:.Object@fpLen),
                                          c((dim(N.oldSig)[3]-1):0))
              } else {
                dimnames(N.oldSig) = list(c("N.A", "N.C", "N.G", "N.T"),
                                        c(1:.Object@fpLen),
                                        NULL)
              }
              
              .Object@N.oldSig = N.oldSig
            } 
            validObject(.Object)
            .Object
          })


#' Summary method for \linkS4class{N}
#' 
#' Defines a summary method for class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "summary",
  signature = "N",
  definition = function(object) {
    cat("An object of class '", class(object), "'\n", sep = "")
    cat("Fits ",
        length(object@N.set),
        " nucleotides at positions ",
        paste(object@N.set, collapse = ", "),
        " for a feature model of length ",
        (object@fS.upFootprintExtend+object@fS.downFootprintExtend+object@seedLen),
        ".\n",
        sep = "")
    if (!all(object@N.equivMat == 0)) {
      if (all(object@N.equivMat[,1:(object@fpLen %/%2)] == 0)) {
        halfMatCheck = object@N.equivMat[1:(object@fpLen %/%2+object@fpLen %%2),
                                         (object@fpLen %/%2+1):object@fpLen]
        halfMatCheck = halfMatCheck[,ncol(halfMatCheck):1]
        if (all(diag(halfMatCheck) == -1)) {
          cat("Nucleotide features are reverse complement symmetric.\n")
        }
      }
    }
    cat("Nucleotide beta values:\n")
    print(getValues(object))
    cat("\n")
    cat("Nucleotide beta errors:\n")
    print(getErrors(object))
    cat("\n")
    cat("Nucleotide beta z-scores:\n")
    print(getZ(object))
    cat("\n")
    cat("Nucleotide beta p-values:\n")
    print(getSig(object))
    cat("\n")
    invisible(NULL)
  }
)

#' Show method for \linkS4class{N}
#' 
#' Defines a show method for class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "show",
  signature = "N",
  definition = function(object) {
    cat("An object of class '", class(object), "'\n", sep = "")
    cat("\n")
    cat('Slot "seedLen": ', object@seedLen, '\n')
    cat("\n")
    cat('Slot "N.upFootprintExtend": ', object@N.upFootprintExtend, '\n')
    cat("\n")
    cat('Slot "N.downFootprintExtend": ', object@N.downFootprintExtend, '\n')
    cat("\n")
    cat('Slot "fS.upFootprintExtend": ', object@fS.upFootprintExtend, '\n')
    cat("\n")
    cat('Slot "fS.downFootprintExtend": ', object@fS.downFootprintExtend, '\n')
    cat("\n")
    cat('Slot "fpLen": ', object@fpLen, '\n')
    cat("\n")
    cat('Slot "N.set": ', object@N.set, '\n')
    cat("\n")
    cat('Slot "N.equivMat":\n')
    if (all(object@N.equivMat == 0)) {
      cat(object@fpLen," x ", object@fpLen, " null equivalence matrix")
    } else {
      print(object@N.equivMat)
    }
    cat("\n")
    cat("\n")
    cat('Slot "N.values":\n')
    print(object@N.values)
    cat("\n")
    cat("\n")
    cat('Slot "N.errors":\n')
    print(object@N.errors)
    cat("\n")
    cat("\n")
    cat('Slot "N.z":\n')
    print(object@N.z)
    cat("\n")
    cat("\n")
    cat('Slot "N.sig":\n')
    print(object@N.sig)
    cat("\n")
    cat("\n")
    
    if (length(dim(object@N.oldValues)) == 3) {
      cat('Slot "N.oldValues":\n')
      if (dim(object@N.oldValues)[3] == 0) {
        cat('<',
            dim(object@N.oldValues)[1],
            ' x ',
            dim(object@N.oldValues)[2],
            ' x ',
            dim(object@N.oldValues)[3],
            ' array of double>','\n', sep = "")
      } else {
        print(object@N.oldValues)
        cat("\n")
      }
    } else {
      cat('Slot "N.oldValues": Not initialized')
      cat("\n")
    }
    
    cat("\n")
    
    if (length(dim(object@N.oldErrors)) == 3) {
      cat('Slot "N.oldErrors":\n')
      if (dim(object@N.oldErrors)[3] == 0) {
        cat('<',
            dim(object@N.oldErrors)[1],
            ' x ',
            dim(object@N.oldErrors)[2],
            ' x ',
            dim(object@N.oldErrors)[3],
            ' array of double>','\n', sep = "")
      } else {
        print(object@N.oldErrors)
        cat("\n")
      }
    } else {
      cat('Slot "N.oldErrors": Not initialized')
      cat("\n")
    }
    cat("\n")
    
    if (length(dim(object@N.oldZ)) == 3) {
      cat('Slot "N.oldZ":\n')
      if (dim(object@N.oldZ)[3] == 0) {
        cat('<',
            dim(object@N.oldZ)[1],
            ' x ',
            dim(object@N.oldZ)[2],
            ' x ',
            dim(object@N.oldZ)[3],
            ' array of double>','\n', sep = "")
      } else {
        print(object@N.oldZ)
        cat("\n")
      }
    } else {
      cat('Slot "N.oldZ": Not initialized')
      cat("\n")
    }
    cat("\n")
    
    if (length(dim(object@N.oldSig)) == 3) {
      cat('Slot "N.oldSig":\n')
      if (dim(object@N.oldSig)[3] == 0) {
        cat('<',
            dim(object@N.oldSig)[1],
            ' x ',
            dim(object@N.oldSig)[2],
            ' x ',
            dim(object@N.oldSig)[3],
            ' array of double>','\n', sep = "")
      } else {
        print(object@N.oldSig)
        cat("\n")
      }
    } else {
      cat('Slot "N.oldSig": Not initialized')
      cat("\n")
    }
    cat("\n")
    invisible(NULL)
  }
)

#' Print method for \linkS4class{N}
#' 
#' Defines a print method for class \linkS4class{N}.
#' 
#' @param x Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "print",
  signature = "N",
  definition = function(x) {
    show(x)
  }
)

#' Get Feature Design for \linkS4class{N}
#' 
#' Defines a print design method for class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getFeatureDesign",
  signature = "N",
  definition = function(object) {
    cat("An object of class '", class(object), "'\n", sep = "")
    cat("\n")
    cat('Slot "seedLen": ', object@seedLen, '\n')
    cat("\n")
    cat('Slot "N.upFootprintExtend": ', object@N.upFootprintExtend, '\n')
    cat("\n")
    cat('Slot "N.downFootprintExtend": ', object@N.downFootprintExtend, '\n')
    cat("\n")
    cat('Slot "fS.upFootprintExtend": ', object@fS.upFootprintExtend, '\n')
    cat("\n")
    cat('Slot "fS.downFootprintExtend": ', object@fS.downFootprintExtend, '\n')
    cat("\n")
    cat('Slot "N.set": ', object@N.set, '\n')
    cat("\n")
    cat("\n")
    invisible(NULL)
  }
)


#' Length method for \linkS4class{N}
#' 
#' Defines a length method for class \linkS4class{N}.
#' 
#' @param x Object of class \linkS4class{N}.
#' @return Footprint length for mono-nucleotides (i.e. the PSAM length).
#' @export
setMethod(
  f = "length",
  signature = "N",
  definition = function(x) {
    return (x@fpLen)
  }
)

#' Get beta values for \linkS4class{N}
#' 
#' Gets N.values slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getValues",
  signature = "N",
  definition = function(object) {
    return(object@N.values)
  }
)

#' Get beta errors for \linkS4class{N}
#' 
#' Gets N.errors slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getErrors",
  signature = "N",
  definition = function(object) {
    return(object@N.errors)
  }
)

#' Get beta z-scores for class \linkS4class{N}
#' 
#' Gets N.z slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getZ",
  signature = "N",
  definition = function(object) {
    return(object@N.z)
  }
)

#' Get beta p-values for \linkS4class{N}
#' 
#' Gets N.sig slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getSig",
  signature = "N",
  definition = function(object) {
    return(object@N.sig)
  }
)

#' Get beta values from previous iterations for \linkS4class{N}
#' 
#' Gets N.oldValues slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getOldValues",
  signature = "N",
  definition = function(object) {
    return(object@N.oldValues)
  }
)

#' Get beta errors from previous iterations for \linkS4class{N}
#' 
#' Gets N.oldErrors slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getOldErrors",
  signature = "N",
  definition = function(object) {
    return(object@N.oldErrors)
  }
)

#' Get beta z-scores from previous iterations for \linkS4class{N}
#' 
#' Gets N.oldZ slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getOldZ",
  signature = "N",
  definition = function(object) {
    return(object@N.oldZ)
  }
)

#' Get beta p-values from previous iterations for \linkS4class{N}
#' 
#' Gets N.oldSig slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getOldSig",
  signature = "N",
  definition = function(object) {
    return(object@N.oldSig)
  }
)

#' Get equivalence matrix for \linkS4class{N}
#' 
#' Gets N.equivMat slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getEquivMat",
  signature = "N",
  definition = function(object) {
    return(object@N.equivMat)
  }
)

#' Get Seed Length
#' 
#' Gets seedLen slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getSeedLen",
  signature = "N",
  definition = function(object) {
    return(object@seedLen)
  }
)

#' Get Length of Upstream Footprint Extension for Nucleotides
#' 
#' Gets N.upFootprintExtend slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getUpFootprintExtend",
  signature = "N",
  definition = function(object) {
    return(object@N.upFootprintExtend)
  }
)

#' Get Length of Upstream Footprint Extension for Full Feature Set
#' 
#' Gets fS.upFootprintExtend slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getFsUpFootprintExtend",
  signature = "N",
  definition = function(object) {
    return(object@fS.upFootprintExtend)
  }
)

#' Get Length of Downstream Footprint Extension for Nucleotides
#' 
#' Gets N.downFootprintExtend slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getDownFootprintExtend",
  signature = "N",
  definition = function(object) {
    return(object@N.downFootprintExtend)
  }
)


#' Get Length of Upstream Footprint Extension for Full Feature Set
#' 
#' Gets fS.downFootprintExtend slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getFsDownFootprintExtend",
  signature = "N",
  definition = function(object) {
    return(object@fS.downFootprintExtend)
  }
)

#' Get Footprint Length
#' 
#' Gets fpLen slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getFpLen",
  signature = "N",
  definition = function(object) {
    return(object@fpLen)
  }
)

#' Get Positions for Mononucleotides Fitting
#' 
#' Gets N.set slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export
setMethod(
  f = "getPositions",
  signature = "N",
  definition = function(object) {
    return(object@N.set)
  }
)

#' Set Values for \linkS4class{N}
#' 
#' Sets N.values slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Matrix of 4xfpLen containing beta values for mono-nucleotides in -ddG units.
#' @export
setReplaceMethod(
  f = "setValues",
  signature = "N",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@N.values)) {
      stop('Current and replacement versions of N.values must have the same dimensions.')
    }
    rownames(value) = c("N.A", "N.C", "N.G", "N.T")
    colnames(value) = 1:object@fpLen
    object@N.values <- value
    validObject(object)
    object
  }
)

#' Set Errors for \linkS4class{N}
#' 
#' Sets N.errors slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Matrix of 4xfpLen containing beta errors for mono-nucleotides in -ddG units
#' @export
setReplaceMethod(
  f = "setErrors",
  signature = "N",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@N.errors)) {
      stop('Current and replacement versions of N.errors must have the same dimensions.')
    }
    rownames(value) = c("N.A", "N.C", "N.G", "N.T")
    colnames(value) = 1:object@fpLen
    object@N.errors <- value
    validObject(object)
    object
  }
)

#' Set Z-scores for \linkS4class{N}
#' 
#' Sets N.z slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Matrix of 4xfpLen containing beta z-scores for mono-nucleotides in -ddG units
#' @export
setReplaceMethod(
  f = "setZ",
  signature = "N",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@N.z)) {
      stop('Current and replacement versions of N.z must have the same dimensions.')
    }
    rownames(value) = c("N.A", "N.C", "N.G", "N.T")
    colnames(value) = 1:object@fpLen
    object@N.z <- value
    validObject(object)
    object
  }
)


#' Set P-values for \linkS4class{N}
#' 
#' Sets N.sig slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Matrix of 4xfpLen containing beta p-values for mono-nucleotides in -ddG units
#' @export
setReplaceMethod(
  f = "setSig",
  signature = "N",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@N.sig)) {
      stop('Current and replacement versions of N.sig must have the same dimensions.')
    }
    rownames(value) = c("N.A", "N.C", "N.G", "N.T")
    colnames(value) = 1:object@fpLen
    object@N.sig <- value
    validObject(object)
    object
  }
)

#' Set Values for Previous Iterations for \linkS4class{N}
#' 
#' Sets N.oldValues slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Array of 4 x fpLen x (iterations-1) containing beta values for mono-nucleotides in -ddG units
#' @export
setReplaceMethod(
  f = "setOldValues",
  signature = "N",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@N.oldValues)) {
      stop('Current and replacement versions of N.oldValues must have the same dimensions.')
    }
    if (dim(value)[3] > 0) {
      dimnames(value) = list(c("N.A", "N.C", "N.G", "N.T"),
                             c(1:object@fpLen),
                             c((dim(value)[3]-1):0))
    } else {
      dimnames(value) = list(c("N.A", "N.C", "N.G", "N.T"),
                             c(1:object@fpLen),
                             NULL)
    }
    
    object@N.oldValues <- value
    validObject(object)
    object
  }
)

#' Set Errors for Previous Iterations for \linkS4class{N}
#' 
#' Sets N.oldErrors slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Array of 4 x fpLen x (iterations-1) containing beta errors for mono-nucleotides in -ddG units
#' @export
setReplaceMethod(
  f = "setOldErrors",
  signature = "N",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@N.oldErrors)) {
      stop('Current and replacement versions of N.oldErrors must have the same dimensions.')
    }
    if (dim(value)[3] > 0) {
      dimnames(value) = list(c("N.A", "N.C", "N.G", "N.T"),
                             c(1:object@fpLen),
                             c((dim(value)[3]-1):0))
    } else {
      dimnames(value) = list(c("N.A", "N.C", "N.G", "N.T"),
                             c(1:object@fpLen),
                             NULL)
    }
    object@N.oldErrors <- value
    validObject(object)
    object
  }
)

#' Set Z-scores for Previous Iterations for \linkS4class{N}
#' 
#' Sets N.oldZ slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Array of 4 x fpLen x (iterations-1) containing beta z-scores for mono-nucleotides in -ddG units
#' @export
setReplaceMethod(
  f = "setOldZ",
  signature = "N",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@N.oldZ)) {
      stop('Current and replacement versions of N.oldZ must have the same dimensions.')
    }
    if (dim(value)[3] > 0) {
      dimnames(value) = list(c("N.A", "N.C", "N.G", "N.T"),
                             c(1:object@fpLen),
                             c((dim(value)[3]-1):0))
    } else {
      dimnames(value) = list(c("N.A", "N.C", "N.G", "N.T"),
                             c(1:object@fpLen),
                             NULL)
    }
    object@N.oldZ <- value
    validObject(object)
    object
  }
)

#' Set P-values for Previous Iterations for \linkS4class{N}
#' 
#' Sets N.oldSig slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Array of 4 x fpLen x (iterations-1) containing beta p-values for mono-nucleotides in -ddG units
#' @export
setReplaceMethod(
  f = "setOldSig",
  signature = "N",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@N.oldSig)) {
      stop('Current and replacement versions of N.oldSig must have the same dimensions.')
    }
    if (dim(value)[3] > 0) {
      dimnames(value) = list(c("N.A", "N.C", "N.G", "N.T"),
                             c(1:object@fpLen),
                             c((dim(value)[3]-1):0))
    } else {
      dimnames(value) = list(c("N.A", "N.C", "N.G", "N.T"),
                             c(1:object@fpLen),
                             NULL)
    }
    object@N.oldSig <- value
    validObject(object)
    object
  }
)

#' Set Seed Length
#' 
#' Sets seedLen slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Seed length for 'model' object to which \linkS4class{N} object belongs.
#' @export
setReplaceMethod(
  f = "setSeedLen",
  signature = "N",
  definition = function(object, value) {
    stop('Changes to the seed length cannot be made at the level of "N" class objects. Note that any changes in the seed length must be made to "model" class objects before any beta information is added.')
    validObject(object)
    object
  }
)

#' Set Length of Upstream Footprint Extension for Nucleotides
#' 
#' Sets N.upFootprintExtend slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Number of bp upstream of seed model to fit mono-nucleotides.
#' @export
setReplaceMethod(
  f = "setUpFootprintExtend",
  signature = "N",
  definition = function(object, value) {
    if (value > object@fS.upFootprintExtend) {
      stop('N.upFootprintExtend must be <= fS.upFootprintExtend')
    }
    if (value < 0) {
      stop('N.upFootprintExtend must be >= 0.')
    }
    if (value == object@N.upFootprintExtend) {
      validObject(object)
      return(object)
    } else if (value > object@N.upFootprintExtend) {
      diff = value-object@N.upFootprintExtend
      min = object@fS.upFootprintExtend-value+1
      seqAdd = seq(min, min+diff-1)
      object@N.set = unique(c(seqAdd, object@N.set))[order(unique(c(seqAdd, object@N.set)))]
      if (value > 0) {
        fullSeqNew = min-1+seq(1, value)
        notIncluded = fullSeqNew[!fullSeqNew %in% object@N.set]
        if (length(notIncluded) > 0) {
          message(paste('Based on previous value for N.set, the following bases implied by the new value for N.upFootprintExtend are excluded from N.set: ',
                        paste(notIncluded, collapse = ", "),
                        sep = ""))
        }
      }
      
      
    } else {
      if (!all(object@N.values == 0)) {
        stop('Cannot deprecate value of N.upFootprintExtend after beta values have been added.')
      } else {
        newSetMin = object@fS.upFootprintExtend-value+1
        object@N.set = object@N.set[object@N.set >= newSetMin]
        if (!(newSetMin %in% object@N.set)) {
          message(paste('New N.upFootprintExtend adds the following nucleotide values not previously included in N.set: ',
                        newSetMin,
                        sep = ""))
        }
        object@N.set = unique(c(newSetMin, object@N.set))[order(unique(c(newSetMin, object@N.set)))]
        if (value > 0) {
          fullSeqNew = newSetMin-1+seq(1, value)
          notIncluded = fullSeqNew[!fullSeqNew %in% object@N.set]
          if (length(notIncluded) > 0) {
            message(paste('Based on previous value for N.set, the following bases implied by the new value for N.upFootprintExtend are excluded from N.set: ',
                          paste(notIncluded, collapse = ", "),
                          sep = ""))
          }
        }
        
      }
    }
    object@N.upFootprintExtend = value
    message('N.downFootprintExtend is not changed automatically to match N.upFootprintExtend')
    validObject(object)
    object
  }
)

#' Set Length of Upstream Footprint Extension for Full Feature Set
#' 
#' Sets fS.upFootprintExtend slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Number of bp upstream of seed model to fit any features.
#' @export
setReplaceMethod(
  f = "setFsUpFootprintExtend",
  signature = "N",
  definition = function(object, value) {
    stop('Changes to fS.upFootprintExtend cannot be made at the level of "N" class objects. Note that any changes in fS.upFootprintExtend must be made to "model" class objects before any beta information is added.')
    validObject(object)
    object
  }
)

#' Set Length of Downstream Footprint Extension for Nucleotides
#' 
#' Sets N.downFootprintExtend slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Number of bp downstream of seed model to fit mono-nucleotides.
#' @export
setReplaceMethod(
  f = "setDownFootprintExtend",
  signature = "N",
  definition = function(object, value) {
    if (value > object@fS.downFootprintExtend) {
      stop('N.downFootprintExtend must be <= fS.downFootprintExtend')
    }
    if (value < 0) {
      stop('N.downFootprintExtend must be >= 0.')
    }
    if (value == object@N.downFootprintExtend) {
      validObject(object)
      return(object)
    } else if (value > object@N.downFootprintExtend) {
      diff = value-object@N.downFootprintExtend
      min = object@fS.upFootprintExtend+object@seedLen+object@N.downFootprintExtend+1
      seqAdd = seq(min, min+diff-1)
      object@N.set = unique(c(seqAdd, object@N.set))[order(unique(c(seqAdd, object@N.set)))]
      
      fullSeqNew = seq(object@fS.upFootprintExtend+object@seedLen+1, min+diff-1)
      notIncluded = fullSeqNew[!fullSeqNew %in% object@N.set]
      if (length(notIncluded) > 0) {
        message(paste('Based on previous value for N.set, the following bases implied by the new value for N.downFootprintExtend are excluded from N.set: ',
                      paste(notIncluded, collapse = ", "),
                      sep = ""))
      }
      
    } else {
      if (!all(object@N.values == 0)) {
        stop('Cannot deprecate value of N.downFootprintExtend after beta values have been added.')
      } else {
        newSetMax = object@fS.upFootprintExtend+object@seedLen+value
        object@N.set = object@N.set[object@N.set <= newSetMax]
        if (!(newSetMax %in% object@N.set)) {
          message(paste('New N.downFootprintExtend adds the following nucleotide values not previously included in N.set: ',
                        newSetMax,
                        sep = ""))
        }
        object@N.set = unique(c(newSetMax, object@N.set))[order(unique(c(newSetMax, object@N.set)))]
        if (value > 0) {
          fullSeqNew = object@fS.upFootprintExtend+object@seedLen+seq(1, value)
          notIncluded = fullSeqNew[!fullSeqNew %in% object@N.set]
          if (length(notIncluded) > 0) {
            message(paste('Based on previous value for N.set, the following bases implied by the new value for N.downFootprintExtend are excluded from N.set: ',
                          paste(notIncluded, collapse = ", "),
                          sep = ""))
          } 
        }
        
      }
    } 
    object@N.downFootprintExtend = value
    message('N.upFootprintExtend is not changed automatically to match N.downFootprintExtend')
    validObject(object)
    object
  }
)

#' Set Length of Downstream Footprint Extension for Full Feature Set
#' 
#' Sets fS.downFootprintExtend slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Number of bp downstream of seed model to fit any features.
#' @export
setReplaceMethod(
  f = "setFsDownFootprintExtend",
  signature = "N",
  definition = function(object, value) {
    stop('Changes to fS.downFootprintExtend cannot be made at the level of "N" class objects. Note that any changes in fS.downFootprintExtend must be made to "model" class objects before any beta information is added.')
    validObject(object)
    object
  }
)

#' Set Footprint Length
#' 
#' Sets fpLen slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Footprint length for full \linkS4class{featureSet} object.
#' @export
setReplaceMethod(
  f = "setFpLen",
  signature = "N",
  definition = function(object, value) {
    stop('Changes to fpLen cannot be made at the level of "N" class objects. Note that any changes in fpLen must be made to "model" class objects before any beta information is added.')
    validObject(object)
    object
  }
)

#' Set Positions for Mononucleotides for \linkS4class{N}
#' 
#' Sets N.set slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Vector of positions within \linkS4class{featureSet} for which mononucleotide features should be fit. 
#' @export
setReplaceMethod(
  f = "setPositions",
  signature = "N",
  definition = function(object, value) {
    if (max(value) > object@fpLen) {
      stop('New N.set value is not consistent with fpLen. ')
    }
    if (!all(object@N.set %in% value)) {
      if (!all(object@N.values == 0)) {
        stop('Cannot deprecate value of N.set after beta values have been added.')
      }
    }
    minSet = min(value)
    if (minSet != min(object@N.set)) {
      object@N.upFootprintExtend = max(object@fS.upFootprintExtend-minSet+1,0)
      message('New N.set implies new value for N.upFootprintExtend: ', object@N.upFootprintExtend)
    }
    maxSet = max(value)
    if (maxSet != max(object@N.set)) {
      object@N.downFootprintExtend = max(maxSet-object@fS.upFootprintExtend-object@seedLen,0)
      message('New N.set implies new value for N.downFootprintExtend: ', object@N.downFootprintExtend)
    }
    
    
    object@N.set = value
    validObject(object)
    object
  }
)

#' Update Values for \linkS4class{N}
#' 
#' Updates N.values slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Matrix of 4xfpLen containing beta values for mono-nucleotides in -ddG units.
setReplaceMethod(
  f = "updateValues",
  signature = "N",
  definition = function(object, value){
    oldValues = getOldValues(object)
    # N object stores values in matrix and oldValues in an array
    updateDim = dim(oldValues)+c(0,0,1)
    oV.updated = array(0, dim = updateDim)
    if (dim(oldValues)[3] > 0) {
      oV.updated[,,2:(dim(oV.updated)[3])] = oldValues 
    }
    oV.updated[,,1] = getValues(object)
    # N specific dimensionality for labels
    dimnames(oV.updated) = list(rownames(getValues(object)), colnames(getValues(object)), c((updateDim[3]-1):0))
    colnames(value) = colnames(getValues(object))
    rownames(value) = rownames(getValues(object))
    object@N.oldValues <- oV.updated
    object@N.values <- value
    object
  }
)

#' Update Errors for \linkS4class{N}
#' 
#' Updates N.errors slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Matrix of 4xfpLen containing beta errors for mono-nucleotides in -ddG units.
setReplaceMethod(
  f = "updateErrors",
  signature = "N",
  definition = function(object, value){
    oldErrors = getOldErrors(object)
    # N object stores values in matrix and oldValues in an array
    updateDim = dim(oldErrors)+c(0,0,1)
    oV.updated = array(0, dim = updateDim)
    if (dim(oldErrors)[3] > 0) {
      oV.updated[,,2:(dim(oV.updated)[3])] = oldErrors
    }
    oV.updated[,,1] = getErrors(object)
    # N specific dimensionality for labels
    dimnames(oV.updated) = list(rownames(getErrors(object)), colnames(getErrors(object)), c((updateDim[3]-1):0))
    colnames(value) = colnames(getErrors(object))
    rownames(value) = rownames(getErrors(object))
    object@N.oldErrors <- oV.updated
    object@N.errors <- value
    object
  }
)

#' Update Z-scores for \linkS4class{N}
#' 
#' Updates N.z slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Matrix of 4xfpLen containing beta z-scores for mono-nucleotides in -ddG units.
setReplaceMethod(
  f = "updateZ",
  signature = "N",
  definition = function(object, value){
    oldValues = getOldZ(object)
    # N object stores values in matrix and oldValues in an array
    updateDim = dim(oldValues)+c(0,0,1)
    oV.updated = array(0, dim = updateDim)
    if (dim(oldValues)[3] > 0) {
      oV.updated[,,2:(dim(oV.updated)[3])] = oldValues 
    }
    oV.updated[,,1] = getZ(object)
    # N specific dimensionality for labels
    dimnames(oV.updated) = list(rownames(getZ(object)), colnames(getZ(object)), c((updateDim[3]-1):0))
    colnames(value) = colnames(getZ(object))
    rownames(value) = rownames(getZ(object))
    object@N.oldZ <- oV.updated
    object@N.z <- value
    object
  }
)

#' Update P-values for class \linkS4class{N}
#' 
#' Updates N.sig slot from object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Matrix of 4xfpLen containing beta p-values for mono-nucleotides in -ddG units. 
setReplaceMethod(
  f = "updateSig",
  signature = "N",
  definition = function(object, value){
    oldValues = getOldSig(object)
    # N object stores values in matrix and oldValues in an array
    updateDim = dim(oldValues)+c(0,0,1)
    oV.updated = array(0, dim = updateDim)
    if (dim(oldValues)[3] > 0) {
      oV.updated[,,2:(dim(oV.updated)[3])] = oldValues 
    }
    oV.updated[,,1] = getSig(object)
    # N specific dimensionality for labels
    dimnames(oV.updated) = list(rownames(getSig(object)), colnames(getSig(object)), c((updateDim[3]-1):0))
    colnames(value) = colnames(getSig(object))
    rownames(value) = rownames(getSig(object))
    object@N.oldSig <- oV.updated
    object@N.sig <- value
    object
  }
)

#' Add New Betas to \linkS4class{N}
#' 
#' Adds new betas to object of class \linkS4class{N} using the output of a glm fit. 
#' 
#' @param object Object of class \linkS4class{N}.
#' @param design Design table used to fit glm model. 
#' @param value Glm fit object. 
#' @param missingValueSuppression DdG suppression of features not occurring in the design dataframe relative to the lowest value betas that are fit. 
#' @return Object of class \linkS4class{N}, updated to include glm fit data. 
#' @export
setMethod(
  f = "addNewBetas",
  signature = "N",
  definition = function(object, design, value, missingValueSuppression, useFixedValuesOffset.N){
    betaSummary = summary(value)$coefficients
    N.set = object@N.set
    N.notSet = c(1:object@fpLen)[!(c(1:object@fpLen) %in% N.set)]
    N.ncol = object@fpLen
    N.values = matrix(NA, nrow = 4, ncol = N.ncol,
                      dimnames = list(c("N.A", "N.C", "N.G", "N.T"), 1:object@fpLen))
    N.errors = matrix(NA, nrow = 4, ncol = N.ncol,
                      dimnames = list(c("N.A", "N.C", "N.G", "N.T"), 1:object@fpLen))
    N.z = matrix(NA, nrow = 4, ncol = N.ncol,
                 dimnames = list(c("N.A", "N.C", "N.G", "N.T"), 1:object@fpLen))
    N.sig = matrix(NA, nrow = 4, ncol = N.ncol,
                   dimnames = list(c("N.A", "N.C", "N.G", "N.T"), 1:object@fpLen))
    if (length(N.notSet) > 0) {
      if (useFixedValuesOffset.N == TRUE) {
        N.values[,N.notSet] = getValues(object)[,N.notSet]
        N.errors[,N.notSet] = getErrors(object)[,N.notSet]
        N.z[,N.notSet] = getZ(object)[,N.notSet]
        N.sig[,N.notSet] = getSig(object)[,N.notSet]
      } else {
        N.values[,N.notSet] = 0
        N.errors[,N.notSet] = 0
        N.z[,N.notSet] = 0
        N.sig[,N.notSet] = 0
      }
    }
    if (!all(object@N.set == 0)) {
      sumN = getDesignMatrix(object, design)
      refN = apply(sumN, 1, which.max)
      zeroNs = matrix(c(refN, as.numeric(names(refN))), nrow = length(refN),
                      ncol = 2, byrow = FALSE)
      N.values[zeroNs] = 0
      N.errors[zeroNs] = 0
      N.z[zeroNs] = 0
      N.sig[zeroNs] = 0
      
      grStr =  "N.[ACGT][0-9]{1,6}"
      Ns = rownames(betaSummary)[grep(grStr, rownames(betaSummary))]
      cols = as.numeric(gsub("N.[ACGT]", "", Ns))
      rows = gsub("[0-9]{1,6}", "", Ns)
      names(cols) = Ns
      names(rows) = Ns
      for (n in Ns) {
        N.values[rows[n], cols[n]] = betaSummary[n,1]
        N.errors[rows[n], cols[n]] = betaSummary[n,2]
        N.z[rows[n], cols[n]] = betaSummary[n,3]
        N.sig[rows[n], cols[n]] = betaSummary[n,4]
      } 
      if (!(all(object@N.equivMat == 0))) {
        equivMat = object@N.equivMat
        for (i in rev(object@N.set)) {
          matchCol = equivMat[,i]
          matchIndices = which(matchCol != 0)
          if (length(matchIndices[matchIndices %in% object@N.set]) > 0) {
            matchIndex = min(matchIndices[matchIndices %in% object@N.set])
            if (matchCol[matchIndex] == 1) {
              N.values[,i] = N.values[,matchIndex]
              N.errors[,i] = N.errors[,matchIndex]
              N.z[,i] = N.z[,matchIndex]
              N.sig[,i] = N.sig[,matchIndex]
              equivMat[i,] = 0
            } else if (matchCol[matchIndex] == -1) {
              if (matchIndex != i) {
                N.values[,i] = rev(N.values[,matchIndex])
                N.errors[,i] = rev(N.errors[,matchIndex])
                N.z[,i] = rev(N.z[,matchIndex])
                N.sig[,i] = rev(N.sig[,matchIndex])
                equivMat[i,] = 0
              } else {
                N.values[3:4,i] = rev(N.values[1:2,i])
                N.errors[3:4,i] = rev(N.errors[1:2,i])
                N.z[3:4,i] = rev(N.z[1:2,i])
                N.sig[3:4,i] = rev(N.sig[1:2,i])
              }
            }
          }
        }
      }
      badN = which(is.na(N.values[,N.set]), arr.ind = TRUE)
      if (length(badN) > 0) {
        if (nrow(badN) < length(N.values[,N.set])) {
          minN = which(N.values[,N.set] == min(N.values[,N.set], na.rm = TRUE), arr.ind = TRUE)[1,]
          N.values[,N.set][badN] = c(N.values[,N.set][minN]-missingValueSuppression)[1]
          N.errors[,N.set][badN] = c(N.errors[,N.set][minN]+missingValueSuppression)[1]
          N.z[,N.set][badN] = 0
          N.sig[,N.set][badN] = 1
        }  else if (nrow(badN) == length(N.values[,N.set])) {
        minN = which(N.values == min(N.values, na.rm = TRUE), arr.ind = TRUE)[1,]
        N.values[,N.set][badN] = c(N.values[minN]-missingValueSuppression)[1]
        N.errors[,N.set][badN] = c(N.errors[minN]+missingValueSuppression)[1]
        N.z[,N.set][badN] = 0
        N.sig[,N.set][badN] = 1
        
        }
      }
      
    }
    
    
    updateValues(object) = N.values
    updateErrors(object) = N.errors
    updateZ(object)= N.z
    updateSig(object) = N.sig
    validObject(object)
    return (object)
  }
)

#' Build Reverse Complement Symmetric Equivalence Matrix for \linkS4class{N}
#' 
#' Builds equivalence matrix encoding reverse complement symmetry for an object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}
setMethod(
  f = "buildSymmetricEquivalenceMatrix",
  signature = "N",
  definition = function(object) {
    n = object@fpLen
    M = matrix(0, nrow = n, ncol = n)
    for (i in 1:(ceiling(n/2))) {
      M[i, (n+1-i)] = -1
    }
    colnames(M) = c(1:ncol(M))
    rownames(M) = c(1:nrow(M))
    return(M)
  }
) 


#' Build Null Equivalence Matrix for \linkS4class{N}
#' 
#' Builds equivalence matrix encoding no symmetries for an object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}
setMethod(
  f = "buildNullEquivalenceMatrix",
  signature = "N",
  definition = function(object) {
    n = object@fpLen
    M = matrix(0, nrow = n, ncol = n)
    colnames(M) = c(1:ncol(M))
    rownames(M) = c(1:nrow(M))
    return(M)
  }
) 


#' Formats Equivalent Matrix for \linkS4class{N}
#' 
#' Formats equivalence matrix for object of class \linkS4class{N}, checking for conflicting entries and getting rid of redundancies.
#' 
#' @param object Object of class \linkS4class{N}.
setMethod(
  f = "formatEquivalenceMatrix",
  signature = "N",
  definition = function(object) {
    equivMat = getEquivMat(object)
    if (all(abs(equivMat) <= 1) != TRUE) {
      stop('Only 1, 0, -1 can be used to indicate pairing in equivalence matrix.')
    }
    for (i in ncol(equivMat):1) {
      if (length(which(equivMat[,i] != 0)) > 1) {
        matchIndices = which(equivMat[,i] != 0)
        keepIndex = matchIndices[1]
        for (j in rev(matchIndices[2:length(matchIndices)])) {
          if (equivMat[keepIndex, j] == 0) {
            equivMat[keepIndex, j] = equivMat[keepIndex,i]*equivMat[j, i]
            equivMat[j,i] = 0
          } else if (equivMat[keepIndex, j] != equivMat[keepIndex,i]*equivMat[j, i]) {
            stop(sprintf('Equivalence Matrix contains conflicts between entry [', keepIndex, ",", j, "] and entry [", j, ",", i, "]."))
          }
          
        }
      } 
    } 
    colnames(equivMat) = c(1:ncol(equivMat))
    rownames(equivMat) = c(1:nrow(equivMat))
    return(equivMat)
  }
)


#' Sets Equivalence Matrix for \linkS4class{N}
#' 
#' Sets an equivalence matrix for object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param value Equivalence matrix for object of class \linkS4class{N}.
#' @export
setReplaceMethod(
  f = "setEquivMat",
  signature = "N",
  definition = function(object, value) {
    colnames(value) = c(1:ncol(value))
    rownames(value) = c(1:nrow(value))
    object@N.equivMat = value
    object@N.equivMat = formatEquivalenceMatrix(object)
    validObject(object)
    object
  }
)


#' Summary Design for Nucleotide Features
#' 
#' Gets nucleotide design from design matrix, i.e. gets counts for each mono-nucleotide feature.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param design Data frame with design matrix, as generated by addDesignMatrix function.
#' @export
setMethod(
  f = "getDesignMatrix",
  signature = "N",
  definition = function(object, design) {
    if (all(object@N.set == 0)) {
      mat = matrix(0, nrow = 0, ncol = 0)
      message('No nucleotides included in fit.')
      return (mat)
    }
    N.set = object@N.set
    if (!(all(object@N.equivMat == 0))) {
      N.equivMat = object@N.equivMat
      N.eMsub = N.equivMat[N.set,]
      N.eMsub = N.eMsub[,N.set]
      rownames(N.eMsub) = N.set
      colnames(N.eMsub) = N.set
      N.included = c(numeric(0))
      N.selfSym= c(numeric(0))
      for (i in N.set) {
        if (length(which(N.eMsub[,as.character(i)] == 0)) == nrow(N.eMsub)) {
          N.included = c(N.included, i)
        } else if (length(which(N.eMsub[,as.character(i)] == 0)) == (nrow(N.eMsub)-1)) {
          if (N.eMsub[as.character(i), as.character(i)] == -1) {
            N.selfSym = c(N.selfSym, i)
          }
        }
      }
    } else {
      N.included = N.set
      N.selfSym = c(numeric(0))
    }
    
    
    # Remove doubled nucleotides
    
    mat = matrix(0, ncol=4, nrow=(length(N.included)+length(N.selfSym)))
    colnames(mat) = c("N.A", "N.C", "N.G", "N.T")
    rownames(mat) = c(N.included, N.selfSym)
    for (i in N.included) {
      for (j in colnames(mat)) {
        col = paste(j, i, sep = "")
        mat[as.character(i),j] = sum(design[,col]/design$Round)
      }
    }
    
    for (i in N.selfSym) {
      for (j in colnames(mat)[1:2]) {
        col = paste(j, i, sep = "")
        mat[as.character(i),j] = sum(design[,col]/design$Round)
      }
    }
    mat = mat[order(as.numeric(rownames(mat))),]
    return(mat)
  }
)


#' Gets plot values for \linkS4class{N}
#' 
#' Reformats N.values and N.errors information for plotting.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param maxErr Maximum y-value to be allowed for error bars in nucleotide dot plot. Truncates error bars that are not contained within this limit.
#' @param iteration Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.  
#' @return Data frame formatted for plotting.
setMethod(
  f = "getPlotValues",
  signature = "N",
  definition = function(object, maxErr = 1, iteration = NULL) {
    minErr = -(maxErr-1)
    if (is.null(iteration)) {
      values = getValues(object)
      errors = getErrors(object)
    } else {
      if (iteration == dim(object@N.oldValues)[3]) {
        values = getValues(object)
        errors = getErrors(object)
      } else {
        values = getOldValues(object)[,,as.character(iteration)]
        errors = getOldErrors(object)[,,as.character(iteration)]
      }
      
    }
    
    values = t(values)
    errors = t(errors)
    maxInd = apply(values, 1, which.max)
    for (i in rownames(values)) {
      maxVal = values[i,maxInd[i]]
      maxSE = errors[i,maxInd[i]]
      for (j in colnames(values)) {
        values[i, j] = values[i, j]-maxVal
        errors[i,j] = errors[i,j]^2
        errors[i, j] = errors[i, j]+maxSE^2
        errors[i,j] = sqrt(errors[i,j])
      }
    }
    values = exp(values)
    values = as.data.frame(values)
    errors = as.data.frame(errors)
    colnames(values) = c("A", "C", "G", "T")
    colnames(errors) = c("A", "C", "G", "T")
    values$PosID = as.numeric(rownames(values))
    errors$PosID = as.numeric(rownames(errors))
    Values.Plot = reshape2::melt(values, id = c("PosID"))
    Values.Plot = Values.Plot[order(Values.Plot$PosID, -Values.Plot$value), ]
    rownames(Values.Plot) = NULL
    colnames(Values.Plot) = c("PosID", "N", "Affinity")
    Errors.Plot= reshape2::melt(errors, id = c("PosID"))
    rownames(Errors.Plot) = NULL
    colnames(Errors.Plot) = c("PosID", "N", "SE")
    Values.Plot = merge(Values.Plot, Errors.Plot, by = c("PosID", "N"), sort = FALSE)
    Values.Plot$PlotErrorMax = Values.Plot$Affinity*(1+Values.Plot$SE)
    Values.Plot$PlotErrorMin = Values.Plot$Affinity*(1-Values.Plot$SE)
    Values.Plot$PlotErrorMax[Values.Plot$PlotErrorMax> maxErr] = maxErr
    Values.Plot$PlotErrorMin[Values.Plot$PlotErrorMin<minErr] = minErr
    return (Values.Plot)
  }
)


#' Get Plot Labels for \linkS4class{N} Plot
#' 
#' Gets labels for plotting that correspond to the highest mono-nucleotide affinity sequence.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param iteration Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.  
#' @return List of nucleotides to be used for labeling position in model for nucleotide dot plot.
#' @export
setMethod(
  f = "getTopSeqLabels",
  signature = "N",
  definition = function(object, iteration = NULL) {
    PSAM = getPSAM(object, iteration)
    lab.plot = c()
    Ns = c("A", "C", "G", "T")
    for (i in 1:object@fpLen) {
      lab.plot = c(lab.plot, Ns[(PSAM[,i] == 1)][!is.na(Ns[(PSAM[,i] == 1)])][1])
    }
    return (lab.plot)
  }
)


#' Returns Standardized PSAM Plot for \linkS4class{N}
#' 
#' Gets PSAM dot plot from N.values slot.
#' 
#' @param x Object of class \linkS4class{N}.
#' @param Ntitle Title for plot.
#' @param maxErr Maximum height of error bars to use in plotting.
#' @param iter Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.  
#' @param regs Optional parameter giving labels to use for region shading. If used, lengths parameter must also be specified.
#' @param lengths Length of region parameters to be included for shading. 
#' @param ddG logical: If ddG == TRUE, ddG units are used on y-axis. If ddG == FALSE, relative affinities are used for y-axis.
#' @return Ggplot dot plot of mono-nucleotide values
#' @export
setMethod(
  f = "plot",
  signature = "N",
  definition = function(x, Ntitle="", maxErr=1, iter=NULL, regs=NULL, lengths = NULL, ddG = FALSE) {
    values.Plot = getPlotValues(x, maxErr, iter)
    values.Plot = values.Plot[!is.na(values.Plot$Affinity),]
    lab.Plot = getTopSeqLabels(x, iter)
    nc = ncol(getPSAM(x))
    if (ddG == FALSE) {
      if ((is.null(regs)) & (is.null(lengths))) {
        psamPlot = ggplot2::ggplot(values.Plot)+
          ggplot2::geom_point(data = values.Plot, ggplot2::aes(x=PosID, y=Affinity, colour=N), size = 1, position = ggplot2::position_dodge(width = .2)) +
          ggplot2::geom_errorbar(data = values.Plot, ggplot2::aes(x=PosID, ymax = PlotErrorMax, ymin=PlotErrorMin, colour = N), width = 1, position = ggplot2::position_dodge(width = .2))+
          ggplot2::scale_colour_manual(name = "", values=c("green","blue","#FFC000","red"), breaks=c("A", "C", "G", "T"),
                                       labels=c("A  ", "C   ", "G   ", "T   "))+
          ggplot2::scale_x_discrete(limits = 1:nc,
                                    labels=parse(text=lab.Plot))+
          ggplot2::coord_cartesian(ylim = c(-.04, 1.04), expand = FALSE)+
          ggplot2::theme_bw()+
          ggplot2::ggtitle(Ntitle)+
          ggplot2::ylab(expression(paste("Relative Affinity")))+
          ggplot2::theme(axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                         axis.title.y=ggplot2::element_text(size = 7, color = "black"),
                         plot.margin=ggplot2::unit(c(4,18,4,1),"pt"),
                         axis.ticks.length=ggplot2::unit(0.3, "mm"),
                         axis.ticks.x = ggplot2::element_line(size = .5),
                         axis.ticks.y = ggplot2::element_line(size = .5),
                         legend.text = ggplot2::element_text(size=6),
                         axis.text.x = ggplot2::element_text(size = 6,
                                                             vjust = -1),
                         axis.title.x = ggplot2::element_blank(),
                         plot.title = ggplot2::element_text(size = 8, 
                                                            face = "bold"),
                         legend.margin=ggplot2::unit(.1,"cm"), 
                         legend.position = "bottom",
                         legend.direction = "horizontal", 
                         legend.box= "horizontal",
                         legend.title = ggplot2::element_blank(),
                         legend.background = ggplot2::element_rect(
                           fill = "transparent", 
                           size = 1))+
          ggplot2::guides(position = "bottom",
                          colour=ggplot2::guide_legend(order = 2,
                                                       keywidth=0.6,
                                                       keyheight=0.4,
                                                       default.unit="cm"))
        
        
        
      } else {
        if (is.null(lengths)) {
          stop('If regs is specified, lengths for the regions must also be included.')
        } else if (is.null(regs)) {
          numReg = length(lengths[[1]])
          regs = list(as.character(1:numReg))
        }
        rects = getRegions(regs, lengths)
        
        psamPlot = ggplot2::ggplot(values.Plot)+
          ggplot2::geom_rect(data = rects, 
                             ggplot2::aes(xmin = xstart,
                                          xmax = xend,
                                          ymin = -Inf,
                                          ymax = Inf,
                                          fill = col), 
                             alpha = 0.4)+
          ggplot2::geom_point(data = values.Plot, ggplot2::aes(x=PosID, y=Affinity, colour=N), size = 1, position = ggplot2::position_dodge(width = .2)) +
          ggplot2::geom_errorbar(data = values.Plot, ggplot2::aes(x=PosID, ymax = PlotErrorMax, ymin=PlotErrorMin, colour = N), width = 1, position = ggplot2::position_dodge(width = .2))+
          ggplot2::scale_colour_manual(name = "", values=c("green","blue","#FFC000","red"), breaks=c("A", "C", "G", "T"),
                                       labels=c("A  ", "C   ", "G   ", "T   "))+
          
          ggplot2::scale_x_discrete(limits = 1:nc,
                                    labels=parse(text=lab.Plot))+
          ggplot2::coord_cartesian(ylim = c(-.04, 1.04), expand = FALSE)+
          ggplot2::theme_bw()+
          ggplot2::ggtitle(Ntitle)+
          ggplot2::ylab(expression(paste("Relative Affinity")))+
          ggplot2::theme(axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                         axis.title.y=ggplot2::element_text(size = 7, color = "black"),
                         plot.margin=ggplot2::unit(c(4,18,4,1),"pt"),
                         axis.ticks.length=ggplot2::unit(0.3, "mm"),
                         axis.ticks.x = ggplot2::element_line(size = .5),
                         axis.ticks.y = ggplot2::element_line(size = .5),
                         legend.text = ggplot2::element_text(size=6),
                         axis.text.x = ggplot2::element_text(size = 6,
                                                             vjust = -1),
                         axis.title.x = ggplot2::element_blank(),
                         plot.title = ggplot2::element_text(size = 8, 
                                                            face = "bold"),
                         legend.margin=ggplot2::unit(.1,"cm"), 
                         legend.position = "bottom",
                         legend.direction = "horizontal", 
                         legend.box= "vertical",
                         legend.title = ggplot2::element_blank(),
                         legend.background = ggplot2::element_rect(
                           fill = "transparent", 
                           size = 1))
        
        
        if (length(levels(rects$col)) <= 7) {
          fillVals = c(character(0))
          for (l in levels(rects$col)) {
            fillVals = c(fillVals, rects$hexCol[rects$col == l][1])
            names(fillVals)[length(fillVals)] = l
          }
          
          psamPlot <- psamPlot +ggplot2::scale_fill_manual(name = NULL,
                                                           values= fillVals,
                                                           breaks = levels(rects$col),
                                                           labels = levels(rects$col))+
            ggplot2::guides(position = "bottom",
                            colour=ggplot2::guide_legend(order = 1,
                                                         keywidth=0.6,
                                                         keyheight=0.4,
                                                         default.unit="cm"),
                            fill=ggplot2::guide_legend(order = 2,
                                                       keywidth=0.6,
                                                       keyheight=0.4,
                                                       default.unit="cm"))
        }
        
      }
    } else {
      values.Plot$ddG = log(values.Plot$Affinity)
      yrange = max(values.Plot$ddG+values.Plot$SE, na.rm = TRUE)-min(values.Plot$ddG-values.Plot$SE, na.rm = TRUE)
      if ((is.null(regs)) & (is.null(lengths))) {
        psamPlot = ggplot2::ggplot(values.Plot)+
          ggplot2::geom_point(data = values.Plot, ggplot2::aes(x=PosID, y=ddG, colour=N), size = 1, position = ggplot2::position_dodge(width = .2)) +
          ggplot2::geom_errorbar(data = values.Plot, ggplot2::aes(x=PosID, ymax = ddG+SE, ymin=ddG-SE, colour = N), width = 1, position = ggplot2::position_dodge(width = .2))+
          ggplot2::scale_colour_manual(name = "", values=c("green","blue","#FFC000","red"), breaks=c("A", "C", "G", "T"),
                                       labels=c("A  ", "C   ", "G   ", "T   "))+
          ggplot2::scale_x_discrete(limits = 1:nc,
                                    labels=parse(text=lab.Plot))+
          ggplot2::ylab(expression(paste('-', Delta, Delta, 'G')))+
          ggplot2::coord_cartesian(xlim = c(0, max(values.Plot$PosID)+1),
                                   ylim = c(min(values.Plot$ddG-values.Plot$SE)-yMarg((yrange/20), 2),max(values.Plot$ddG+values.Plot$SE)+yMarg((yrange/20), 2)),
                                   expand = FALSE)+
          ggplot2::theme_bw()+
          ggplot2::ggtitle(Ntitle)+
          ggplot2::theme(axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                         axis.title.y=ggplot2::element_text(size = 7, color = "black"),
                         plot.margin=ggplot2::unit(c(4,18,4,1),"pt"),
                         axis.ticks.length=ggplot2::unit(0.3, "mm"),
                         axis.ticks.x = ggplot2::element_line(size = .5),
                         axis.ticks.y = ggplot2::element_line(size = .5),
                         legend.text = ggplot2::element_text(size=6),
                         axis.text.x = ggplot2::element_text(size = 6,
                                                             vjust = -1),
                         axis.title.x = ggplot2::element_blank(),
                         plot.title = ggplot2::element_text(size = 8, 
                                                            face = "bold"),
                         legend.margin=ggplot2::unit(.1,"cm"), 
                         legend.position = "bottom",
                         legend.direction = "horizontal", 
                         legend.box= "horizontal",
                         legend.title = ggplot2::element_blank(),
                         legend.background = ggplot2::element_rect(
                           fill = "transparent", 
                           size = 1))+
          ggplot2::guides(position = "bottom",
                          colour=ggplot2::guide_legend(order = 2,
                                                       keywidth=0.6,
                                                       keyheight=0.4,
                                                       default.unit="cm"))
        
        
        
      } else {
        if (is.null(lengths)) {
          stop('If regs is specified, lengths for the regions must also be included.')
        } else if (is.null(regs)) {
          numReg = length(lengths[[1]])
          regs = list(as.character(1:numReg))
        }
        rects = getRegions(regs, lengths)
        
        psamPlot = ggplot2::ggplot(values.Plot)+
          ggplot2::geom_rect(data = rects, 
                             ggplot2::aes(xmin = xstart,
                                          xmax = xend,
                                          ymin = -Inf,
                                          ymax = Inf, 
                                          fill = col), 
                             alpha = 0.4)+
          ggplot2::geom_point(data = values.Plot, ggplot2::aes(x=PosID, y=ddG, colour=N), size = 1, position = ggplot2::position_dodge(width = .2)) +
          ggplot2::geom_errorbar(data = values.Plot, ggplot2::aes(x=PosID, ymax = ddG+SE, ymin=ddG-SE, colour = N), width = 1, position = ggplot2::position_dodge(width = .2))+
          ggplot2::scale_colour_manual(name = "", values=c("green","blue","#FFC000","red"), breaks=c("A", "C", "G", "T"),
                                       labels=c("A  ", "C   ", "G   ", "T   "))+
          ggplot2::scale_x_discrete(limits = 1:nc,
                                    labels=parse(text=lab.Plot))+
          ggplot2::ylab(expression(paste('-', Delta, Delta, 'G')))+
          ggplot2::coord_cartesian(xlim = c(0, max(values.Plot$PosID)+1),
                                   ylim = c(yMarg(min(values.Plot$ddG-values.Plot$SE, na.rm = TRUE))-yMarg((yrange/20), 2),yMarg(max(values.Plot$ddG+values.Plot$SE, na.rm = TRUE)+yMarg((yrange/20), 2))),
                                   expand = FALSE)+
          ggplot2::theme_bw()+
          ggplot2::ggtitle(Ntitle)+
          ggplot2::theme(axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                         axis.title.y=ggplot2::element_text(size = 7, color = "black"),
                         plot.margin=ggplot2::unit(c(4,18,4,1),"pt"),
                         axis.ticks.length=ggplot2::unit(0.3, "mm"),
                         axis.ticks.x = ggplot2::element_line(size = .5),
                         axis.ticks.y = ggplot2::element_line(size = .5),
                         legend.text = ggplot2::element_text(size=6),
                         axis.text.x = ggplot2::element_text(size = 6,
                                                             vjust = -1),
                         axis.title.x = ggplot2::element_blank(),
                         plot.title = ggplot2::element_text(size = 8, 
                                                            face = "bold"),
                         legend.margin=ggplot2::unit(.1,"cm"), 
                         legend.position = "bottom",
                         legend.direction = "horizontal", 
                         legend.box= "vertical",
                         legend.title = ggplot2::element_blank(),
                         legend.background = ggplot2::element_rect(
                           fill = "transparent", 
                           size = 1))
        
        
        if (length(levels(rects$col)) <= 7) {
          fillVals = c(character(0))
          for (l in levels(rects$col)) {
            fillVals = c(fillVals, rects$hexCol[rects$col == l][1])
            names(fillVals)[length(fillVals)] = l
          }
          
          psamPlot <- psamPlot +ggplot2::scale_fill_manual(name = NULL,
                                                           values= fillVals,
                                                           breaks = levels(rects$col),
                                                           labels = levels(rects$col))+
            ggplot2::guides(position = "bottom",
                            colour=ggplot2::guide_legend(order = 1,
                                                         keywidth=0.6,
                                                         keyheight=0.4,
                                                         default.unit="cm"),
                            fill=ggplot2::guide_legend(order = 2,
                                                       keywidth=0.6,
                                                       keyheight=0.4,
                                                       default.unit="cm"))
        }
        
      }
    }
    
    
    
    return(psamPlot)
  }
)


#' Get Mononucleotide PSAM. 
#' 
#' Reshapes N.values into standard PSAM format.
#' 
#' @param object Object of class \linkS4class{N}.
#' @param iteration Iteration of beta values to be used for PSAM. If iteration == NULL, the most recent round of beta values are used.  
#' @return Data.frame containing PSAM of relative affinity values (exp(-ddG))
#' @export
setMethod(
  f = "getPSAM",
  signature = "N",
  definition = function(object, iteration=NULL){
    if (is.null(iteration)) {
      PSAM = t(object@N.values)
    } else {
      oldIters = as.numeric(dimnames(object@N.oldValues)[[3]])
      if (iteration == max(oldIters)+1) {
        PSAM = t(object@N.values)
      } else if (iteration %in% oldIters) {
        PSAM = t(object@N.oldValues[,,as.character(iteration)])
      } else {
        stop('Invalid value for iteration')
      }
    }
    
    maxNs = apply(PSAM, 1, max, na.rm=TRUE)
    for (i in 1:4) {
      PSAM[,i] = PSAM[,i]-maxNs
    }
    PSAM = t(exp(PSAM))
    rownames(PSAM) = c("A", "C", "G", "T")
    return (PSAM)
  }
)







