#' @include Intercept-validity.R
#' @include generic-functions-definitions.R
NULL

#' 'Intercept' Class Definition
#' 
#' An S4 class to represent intercept features, including round-, view-, and strand-specific intercept contributions.
#' 
#' @slot I.values 3-D array of beta values with entries corresponding to (Strand)x(View)x(Round) intercept values.
#' @slot I.errors 3-D array of errors for I.values.
#' @slot I.z 3-D array of Z-scores for I.values.
#' @slot I.sig 3-D array of p-values for I.values.
#' @slot I.oldValues 4-D array for storing previous iteration values beta values for I.values where each entry corresponds to (Strand)x(View)x(Round)x(Iteration).
#' @slot I.oldErrors 4-D array of previous iteration I.errors.
#' @slot I.oldZ 4-D array of previous iteration I.z.
#' @slot I.oldSig 4-D array of previous iteration I.sig.
#' @export Intercept
Intercept <- setClass("Intercept", 
                      representation(I.values = "array",
                                     I.errors = "array", 
                                     I.z = "array", 
                                     I.sig = "array", 
                                     I.oldValues = "array",
                                     I.oldErrors = "array",
                                     I.oldZ = "array", 
                                     I.oldSig = "array"), 
                      validity = validIntercept)


#' Constructor for \linkS4class{Intercept} Class
#' 
#' Constructor for S4 class  \linkS4class{Intercept}.
#'
#' @param numViews number of views scored on each DNA strand.
#' @param rounds rounds of SELEX data used to evaluate 'model' to which 'Intercept' belongs.
#' @export
setMethod("initialize",
          "Intercept",
          function(.Object,
                   I.values,
                   I.errors,
                   I.z,
                   I.sig, 
                   I.oldValues,
                   I.oldErrors,
                   I.oldZ, 
                   I.oldSig,
                   numViews,
                   rounds, ...) {
            strandLabels = paste("Strand.", c("F", "R"), sep = "")
            if ((missing(numViews)) & (missing(I.values)) & (missing(I.oldValues))) {
              stop("Must give numViews if no number of windows is specified by I.values and/or I.oldValues")
            } else {
              if (!missing(numViews)) {
                viewLabels = paste("View.", 1:numViews, sep = "")
              } else if (!missing(I.values)) {
                numViews = ncol(I.values)
                viewLabels = paste("View.", 1:ncol(I.values), sep = "")
              } else {
                numViews = ncol(I.oldValues)
                viewLabels = paste("View..", 1:ncol(I.oldValues), sep = "")
              }
            }
            if (missing(rounds)) {
              stop("Must specify rounds explicitly")
            }else {
              if (length(rounds[[1]]) < 1) {
                stop("Must specify rounds as a list of length 1 with the 1st entry holding a vector of length >= 1.")
              }
              rounds[[1]] = rounds[[1]][order(rounds[[1]])]
              roundLabels = paste("Round.", rounds[[1]], sep = "")
              if (!missing(I.values)) {
                if (!all(c(2, numViews, length(rounds[[1]])) == dim(I.values))) {
                  stop(paste("dim(I.values) must be equal to c(2, ",numViews,", ",length(rounds[[1]]),").",sep = "" ))
                }
                if (length(dimnames(I.values)[[3]]) > 0) {
                  if (dim(I.values)[3] != length(roundLabels)) {
                    stop("Dimension of I.values does not match length of rounds[[1]].")
                  } else if (!all(unique(roundLabels)[order(unique(roundLabels))] == unique(dimnames(I.values)[[3]])[order(unique(dimnames(I.values)[[3]]))])) {
                    stop("Dimnames[[3]] of I.values does not match rounds[[1]].")
                  } else {
                    I.values = I.values[,,roundLabels,drop = FALSE]
                  }
                } else {
                  stop("If I.values is assigned, dimnames[[3]] must contain labels denoting rounds.")
                }
              }
              
              if (!missing(I.oldValues)) {
                if (!all(c(2, numViews, length(rounds[[1]])) == dim(I.oldValues)[1:3])) {
                  stop(paste("dim(I.oldValues)[1:3] must be equal to c(2, ",numViews,", ",length(rounds[[1]]),").",sep = "" ))
                }
                if (length(dimnames(I.oldValues)[[3]]) != 0) {
                  if (dim(I.oldValues)[3] != length(roundLabels)) {
                    stop("Dimension of I.oldValues does not match length of rounds[[1]].")
                  } else if (!all(unique(roundLabels)[order(unique(roundLabels))] == unique(dimnames(I.oldValues)[[3]])[order(unique(dimnames(I.oldValues)[[3]]))])) {
                    stop("Dimnames[[3]] of I.oldValues does not match rounds[[1]].")
                  } else {
                    I.oldValues = I.oldValues[,,roundLabels,,drop = FALSE]
                  }
                } else {
                  stop("If I.oldValues is assigned, dimnames[[3]] must contain labels denoting rounds.")
                }
              }
              
              if (!missing(I.errors)) {
                if (!all(c(2, numViews, length(rounds[[1]])) == dim(I.errors))) {
                  stop(paste("dim(I.error) must be equal to c(2, ",numViews,", ",length(rounds[[1]]),").",sep = "" ))
                }
                if (length(dimnames(I.errors)[[3]]) != 0) {
                  if (dim(I.errors)[3] != length(roundLabels)) {
                    stop("Dimension of I.errors does not match length of rounds[[1]].")
                  } else if (!all(unique(roundLabels)[order(unique(roundLabels))] == unique(dimnames(I.errors)[[3]])[order(unique(dimnames(I.errors)[[3]]))])) {
                    stop("Dimnames[[3]] of I.errors does not match rounds[[1]].")
                  } else {
                    I.errors = I.errors[,,roundLabels,drop = FALSE]
                  }
                } else {
                  stop("If I.errors is assigned, dimnames[[3]] must contain labels denoting rounds.")
                }
              }
              
              if (!missing(I.oldErrors)) {
                if (!all(c(2, numViews, length(rounds[[1]])) == dim(I.oldErrors)[1:3])) {
                  stop(paste("dim(I.oldErrors)[1:3] must be equal to c(2, ",numViews,", ",length(rounds[[1]]),").",sep = "" ))
                }
                if (length(dimnames(I.oldErrors)[[3]]) != 0) {
                  if (dim(I.oldErrors)[3] != length(roundLabels)) {
                    stop("Dimension of I.oldErrors does not match length of rounds[[1]].")
                  } else if (!all(unique(roundLabels)[order(unique(roundLabels))] == unique(dimnames(I.oldErrors)[[3]])[order(unique(dimnames(I.oldErrors)[[3]]))])) {
                    stop("Dimnames[[3]] of I.oldErrors does not match rounds[[1]].")
                  } else {
                    I.oldErrors = I.oldErrors[,,roundLabels,,drop = FALSE]
                  }
                } else {
                  stop("If I.oldErrors is assigned, dimnames[[3]] must contain labels denoting rounds.")
                }
              }
              
              if (!missing(I.z)) {
                if (!all(c(2, numViews, length(rounds[[1]])) == dim(I.z))) {
                  stop(paste("dim(I.z) must be equal to c(2, ",numViews,", ",length(rounds[[1]]),").",sep = "" ))
                }
                if (length(dimnames(I.z)[[3]]) != 0) {
                  if (dim(I.z)[3] != length(roundLabels)) {
                    stop("Dimension of I.z does not match length of rounds[[1]].")
                  } else if (!all(unique(roundLabels)[order(unique(roundLabels))] == unique(dimnames(I.z)[[3]])[order(unique(dimnames(I.z)[[3]]))])) {
                    stop("Dimnames[[3]] of I.z does not match rounds[[1]].")
                  } else {
                    I.z = I.z[,,roundLabels,drop = FALSE]
                  }
                } else {
                  stop("If I.z is assigned, dimnames[[3]] must contain labels denoting rounds.")
                }
              }
              
              if (!missing(I.oldZ)) {
                if (!all(c(2, numViews, length(rounds[[1]])) == dim(I.oldZ)[1:3])) {
                  stop(paste("dim(I.oldZ)[1:3] must be equal to c(2, ",numViews,", ",length(rounds[[1]]),").",sep = "" ))
                }
                if (length(dimnames(I.oldZ)[[3]]) != 0) {
                  if (dim(I.oldZ)[3] != length(roundLabels)) {
                    stop("Dimension of I.oldZ does not match length of rounds[[1]].")
                  } else if (!all(unique(roundLabels)[order(unique(roundLabels))] == unique(dimnames(I.oldZ)[[3]])[order(unique(dimnames(I.oldZ)[[3]]))])) {
                    stop("Dimnames[[3]] of I.oldZ does not match rounds[[1]].")
                  } else {
                    I.oldZ = I.oldZ[,,roundLabels,,drop = FALSE]
                  }
                } else {
                  stop("If I.oldZ is assigned, dimnames[[3]] must contain labels denoting rounds.")
                }
              }
              
              if (!missing(I.sig)) {
                if (!all(c(2, numViews, length(rounds[[1]])) == dim(I.sig))) {
                  stop(paste("dim(I.sig) must be equal to c(2, ",numViews,", ",length(rounds[[1]]),").",sep = "" ))
                }
                if (length(dimnames(I.sig)[[3]]) != 0) {
                  if (dim(I.sig)[3] != length(roundLabels)) {
                    stop("Dimension of I.sig does not match length of rounds[[1]].")
                  } else if (!all(unique(roundLabels)[order(unique(roundLabels))] == unique(dimnames(I.sig)[[3]])[order(unique(dimnames(I.sig)[[3]]))])) {
                    stop("Dimnames[[3]] of I.sig does not match rounds[[1]].")
                  } else {
                    I.sig = I.sig[,,roundLabels,drop = FALSE]
                  }
                } else {
                  stop("If I.sig is assigned, dimnames[[3]] must contain labels denoting rounds.")
                }
              }
              
              if (!missing(I.oldSig)) {
                if (!all(c(2, numViews, length(rounds[[1]])) == dim(I.oldSig)[1:3])) {
                  stop(paste("dim(I.oldSig)[1:3] must be equal to c(2, ",numViews,", ",length(rounds[[1]]),").",sep = "" ))
                }
                if (length(dimnames(I.oldSig)[[3]]) != 0) {
                  if (dim(I.oldSig)[3] != length(roundLabels)) {
                    stop("Dimension of I.oldSig does not match length of rounds[[1]].")
                  } else if (!all(unique(roundLabels)[order(unique(roundLabels))] == unique(dimnames(I.oldSig)[[3]])[order(unique(dimnames(I.oldSig)[[3]]))])) {
                    stop("Dimnames[[3]] of I.oldSig does not match rounds[[1]].")
                  } else {
                    I.oldSig = I.oldSig[,,roundLabels,,drop = FALSE]
                  }
                } else {
                  stop("If I.oldSig is assigned, dimnames[[3]] must contain labels denoting rounds.")
                }
              }
            }
            
            if (!(missing(I.values))) {
              if (ncol(I.values) != numViews) {
                stop(paste("ncol(I.values) must be equal to numViews = ",numViews,sep = "" ))
              }
              .Object@I.values = I.values
              dimnames(.Object@I.values) = list(strandLabels,
                                                paste("View.", 1:ncol(.Object@I.values), sep = ""),
                                                roundLabels)
            } else {
              .Object@I.values = array(0, dim = c(2, numViews, length(rounds[[1]])), dimnames = list(strandLabels,
                                                                                                     viewLabels,
                                                                                                     roundLabels))
            }
            if (!(missing(I.errors))) {
              if (ncol(I.errors) != numViews) {
                stop(paste("ncol(I.errors) must be equal to numViews = ",numViews,sep = "" ))
              }
              .Object@I.errors = I.errors
              dimnames(.Object@I.errors) = list(strandLabels,
                                                paste("View.", 1:ncol(.Object@I.errors), sep = ""),
                                                roundLabels)
            } else {
              .Object@I.errors = array(0, dim = c(2, numViews, length(rounds[[1]])), dimnames = list(strandLabels,
                                                                                                     viewLabels,
                                                                                                     roundLabels))
            }
            if (!(missing(I.z))) {
              if (ncol(I.z) != numViews) {
                stop(paste("ncol(I.z) must be equal to numViews = ",numViews,sep = "" ))
              }
              .Object@I.z= I.z
              dimnames(.Object@I.z) = list(strandLabels,
                                           paste("View.", 1:ncol(.Object@I.z), sep = ""),
                                           roundLabels)
            } else {
              .Object@I.z = array(0, dim = c(2, numViews, length(rounds[[1]])), dimnames = list(strandLabels,
                                                                                                viewLabels,
                                                                                                roundLabels))
            }
            if (!(missing(I.sig))) {
              if (ncol(I.sig) != numViews) {
                stop(paste("ncol(I.sig) must be equal to numViews = ",numViews,sep = "" ))
              }
              .Object@I.sig = I.sig
              dimnames(.Object@I.sig) = list(strandLabels,
                                             paste("View.", 1:ncol(.Object@I.sig), sep = ""),
                                             roundLabels)
            } else {
              .Object@I.sig = array(0, dim = c(2, numViews, length(rounds[[1]])), dimnames = list(strandLabels,
                                                                                                  viewLabels,
                                                                                                  roundLabels))
            }
            if (!(missing(I.oldValues))) {
              if (ncol(I.oldValues) != numViews) {
                stop(paste("ncol(I.oldValues) must be equal to numViews = ",numViews,sep = "" ))
              }
              .Object@I.oldValues = I.oldValues
              if (dim(.Object@I.oldValues)[4] == 0) {
                dimnames(.Object@I.oldValues) = list(strandLabels,
                                                     paste("View.", 1:ncol(.Object@I.oldValues), sep = ""),
                                                     roundLabels,
                                                     NULL)
              } else {
                dimnames(.Object@I.oldValues) = list(strandLabels,
                                                     paste("View.", 1:ncol(.Object@I.oldValues), sep = ""),
                                                     roundLabels,
                                                     c((dim(.Object@I.oldValues)[4]-1):0))
              }
              
            } else {
              .Object@I.oldValues = array(0, dim = c(2, numViews, length(rounds[[1]]), 0), dimnames = list(strandLabels,
                                                                                                           viewLabels,
                                                                                                           roundLabels,
                                                                                                           NULL))
            }
            
            if (!(missing(I.oldErrors))) {
              if (ncol(I.oldErrors) != numViews) {
                stop(paste("ncol(I.oldErrors) must be equal to numViews = ",numViews,sep = "" ))
              }
              .Object@I.oldErrors = I.oldErrors
              if (dim(.Object@I.oldErrors)[4] == 0) {
                dimnames(.Object@I.oldErrors) = list(strandLabels,
                                                     paste("View.", 1:ncol(.Object@I.oldErrors), sep = ""),
                                                     roundLabels,
                                                     NULL)
              } else {
                dimnames(.Object@I.oldErrors) = list(strandLabels,
                                                     paste("View.", 1:ncol(.Object@I.oldErrors), sep = ""),
                                                     roundLabels,
                                                     c((dim(.Object@I.oldErrors)[4]-1):0))
              }
              
            } else {
              .Object@I.oldErrors = array(0, dim = c(2, numViews, length(rounds[[1]]), 0), dimnames = list(strandLabels,
                                                                                                           viewLabels,
                                                                                                           roundLabels,
                                                                                                           NULL))
            }
            
            if (!(missing(I.oldZ))) {
              if (ncol(I.oldZ) != numViews) {
                stop(paste("ncol(I.oldZ) must be equal to numViews = ",numViews,sep = "" ))
              }
              .Object@I.oldZ = I.oldZ
              if (dim(.Object@I.oldZ)[4] == 0) {
                dimnames(.Object@I.oldZ) = list(strandLabels,
                                                paste("View.", 1:ncol(.Object@I.oldZ), sep = ""),
                                                roundLabels,
                                                NULL)
              } else {
                dimnames(.Object@I.oldZ) = list(strandLabels,
                                                paste("View.", 1:ncol(.Object@I.oldZ), sep = ""),
                                                roundLabels,
                                                c((dim(.Object@I.oldZ)[4]-1):0))
              }
              
            } else {
              .Object@I.oldZ = array(0, dim = c(2, numViews, length(rounds[[1]]), 0), dimnames = list(strandLabels,
                                                                                                      viewLabels,
                                                                                                      roundLabels,
                                                                                                      NULL))
            }
            
            if (!(missing(I.oldSig))) {
              if (ncol(I.oldSig) != numViews) {
                stop(paste("ncol(I.oldSig) must be equal to numViews = ",numViews,sep = "" ))
              }
              .Object@I.oldSig = I.oldSig
              if (dim(.Object@I.oldSig)[4] == 0) {
                dimnames(.Object@I.oldSig) = list(strandLabels,
                                                  paste("View.", 1:ncol(.Object@I.oldSig), sep = ""),
                                                  roundLabels,
                                                  NULL)
              } else {
                dimnames(.Object@I.oldSig) = list(strandLabels,
                                                  paste("View.", 1:ncol(.Object@I.oldSig), sep = ""),
                                                  roundLabels,
                                                  c((dim(.Object@I.oldSig)[4]-1):0))
              }
              
            } else {
              .Object@I.oldSig = array(0, dim = c(2, numViews, length(rounds[[1]]), 0), dimnames = list(strandLabels,
                                                                                                        viewLabels,
                                                                                                        roundLabels,
                                                                                                        NULL))
            }
            
            validObject(.Object)
            .Object
          })

#' Summary Method for  \linkS4class{Intercept}
#' 
#' Defines summary method for object of class  \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export 
setMethod(
  f = "summary",
  signature = "Intercept",
  definition = function(object) {
    cat("An object of class '", class(object), "'\n", sep = "")
    cat("Fits up to ",
        2*ncol(object@I.values),
        " views and ",
        dim(object@I.values)[3],
        " rounds.\n",
        sep = "")
    cat("\n")
    cat("Intercept beta values:\n")
    nR = dim(object@I.values)[3]
    for (i in 1:nR) {
      cat(paste(dimnames(object@I.values)[[3]][i], ":\n", sep = ""))
      print(object@I.values[,,i])
      cat("\n")
    }
    cat("\n")
    cat("Intercept beta errors:\n")
    for (i in 1:nR) {
      cat(paste(dimnames(object@I.errors)[[3]][i], ":\n", sep = ""))
      print(object@I.errors[,,i])
      cat("\n")
    }
    cat("\n")
    cat("Intercept beta z-scores:\n")
    for (i in 1:nR) {
      cat(paste(dimnames(object@I.z)[[3]][i], ":\n", sep = ""))
      print(object@I.z[,,i])
      cat("\n")
    }
    cat("\n")
    cat("Intercept beta p-values:\n")
    for (i in 1:nR) {
      cat(paste(dimnames(object@I.sig)[[3]][i], ":\n", sep = ""))
      print(object@I.sig[,,i])
      cat("\n")
    }
    invisible(NULL)
  }
)

#' Show Method for \linkS4class{Intercept}
#'
#' Defines show method for \linkS4class{Intercept} class.
#'
#' @param object Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "show",
  signature = "Intercept",
  definition = function(object) {
    cat("An object of class '", class(object), "'\n", sep = "")
    cat("\n")
    if (length(dim(object@I.values)) == 3) {
      nR = dim(object@I.values)[3]
      cat("\n")
      cat('Slot "I.values":\n')
      for (i in 1:nR) {
        cat(paste(dimnames(object@I.values)[[3]][i], ":\n", sep = ""))
        print(object@I.values[,,i])
        cat("\n")
      }
      cat("\n")
      cat('Slot "I.errors":\n')
      for (i in 1:nR) {
        cat(paste(dimnames(object@I.values)[[3]][i], ":\n", sep = ""))
        print(object@I.errors[,,i])
        cat("\n")
      }
      cat("\n")
      cat('Slot "I.z":\n')
      for (i in 1:nR) {
        cat(paste(dimnames(object@I.values)[[3]][i], ":\n", sep = ""))
        print(object@I.z[,,i])
        cat("\n")
      }
      cat("\n")
      cat('Slot "I.sig":\n')
      for (i in 1:nR) {
        cat(paste(dimnames(object@I.values)[[3]][i], ":\n", sep = ""))
        print(object@I.sig[,,i])
        cat("\n")
      }
      cat("\n")
      cat('Slot "I.oldValues":\n')
      if (dim(object@I.oldValues)[4] == 0) {
        cat('<',
            dim(object@I.oldValues)[1],
            ' x ',
            dim(object@I.oldValues)[2],
            ' x ',
            dim(object@I.oldValues)[3],
            ' x ',
            dim(object@I.oldValues)[4],
            ' array of double>\n', sep = "")
      } else {
        for (j in dimnames(object@I.oldValues)[[4]]) {
          for (i in 1:nR) {
            cat(paste('Iteration: ', j, ' ,', dimnames(object@I.oldValues)[[3]][i], ":\n", sep = ""))
            print(object@I.oldValues[,,i,j])
            cat("\n")
          }
        }
        
      }
      cat("\n")
      cat('Slot "I.oldErrors":\n')
      if (dim(object@I.oldErrors)[4] == 0) {
        cat('<',
            dim(object@I.oldErrors)[1],
            ' x ',
            dim(object@I.oldErrors)[2],
            ' x ',
            dim(object@I.oldErrors)[3],
            ' x ',
            dim(object@I.oldErrors)[4],
            ' array of double>\n', sep = "")
      } else {
        for (j in dimnames(object@I.oldErrors)[[4]]) {
          for (i in 1:nR) {
            cat(paste('Iteration: ', j, ' ,', dimnames(object@I.oldErrors)[[3]][i], ":\n", sep = ""))
            print(object@I.oldErrors[,,i,j])
            cat("\n")
          }
        }
      }
      cat("\n")
      cat('Slot "I.oldZ":\n')
      if (dim(object@I.oldZ)[4] == 0) {
        cat('<',
            dim(object@I.oldZ)[1],
            ' x ',
            dim(object@I.oldZ)[2],
            ' x ',
            dim(object@I.oldZ)[3],
            ' x ',
            dim(object@I.oldZ)[4],
            ' array of double>\n', sep = "")
      } else {
        for (j in dimnames(object@I.oldZ)[[4]]) {
          for (i in 1:nR) {
            cat(paste('Iteration: ', j, ' ,', dimnames(object@I.oldZ)[[3]][i], ":\n", sep = ""))
            print(object@I.oldZ[,,i,j])
            cat("\n")
          }
        }
      }
      cat("\n")
      cat('Slot "I.oldSig":\n')
      if (dim(object@I.oldSig)[4] == 0) {
        cat('<',
            dim(object@I.oldSig)[1],
            ' x ',
            dim(object@I.oldSig)[2],
            ' x ',
            dim(object@I.oldSig)[3],
            ' x ',
            dim(object@I.oldSig)[4],
            ' array of double>\n', sep = "")
      } else {
        for (j in dimnames(object@I.oldSig)[[4]]) {
          for (i in 1:nR) {
            cat(paste('Iteration: ', j, ' ,', dimnames(object@I.oldSig)[[3]][i], ":\n", sep = ""))
            print(object@I.oldSig[,,i,j])
            cat("\n")
          }
        }
      }
      cat("\n")
      invisible(NULL)
    } else {
      cat('Slot "I.values": Not initialized\n')
      cat('Slot "I.errors": Not initialized\n')
      cat('Slot "I.z": Not initialized\n')
      cat('Slot "I.sig": Not initialized\n')
      cat('Slot "I.oldValues": Not initialized\n')
      cat('Slot "I.oldErrors": Not initialized\n')
      cat('Slot "I.oldZ": Not initialized\n')
      cat('Slot "I.oldSig": Not initialized\n')
      cat("\n")
      invisible(NULL)
    }
  }
)

#' Print Method for \linkS4class{Intercept} 
#' 
#' Defines 'print' method for object of class \linkS4class{Intercept}
#' .
#' @param x Object of class \linkS4class{Intercept}.  
#' @export
setMethod(
  f = "print",
  signature = "Intercept",
  definition = function(x) {
    show(x)
  }
)

#' Length for \linkS4class{Intercept}
#'
#' Defines length method for object of class \linkS4class{Intercept}.
#' 
#' Length of \linkS4class{Intercept} object is defined as the number of views scored on each orientation of a DNA strand.
#'
#' @param x Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "length",
  signature = "Intercept",
  definition = function(x) {
    length = ncol(getValues(x))
    return (length)
  }
)

#' Get Feature Design for \linkS4class{Intercept} 
#' 
#' Defines a feature design print method for class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "getFeatureDesign",
  signature = "Intercept",
  definition = function(object) {
    cat("An object of class '", class(object), "'\n", sep = "")
    cat("Number of Views per Strand of DNA: ", ncol(object@I.values), '\n', sep = "")
    cat("Number of Rounds: ", dim(object@I.values)[3]," (", paste(gsub("Round.", "", dimnames(object@I.values)[[3]]), collapse = ", "),")\n", sep = "")
    cat('Number of iterations completed: ', dim(object@I.oldValues)[4], '\n', sep = "")
    cat("\n")
    invisible(NULL)
  }
)

#' Gets Plot Values for \linkS4class{Intercept}.
#' 
#' Reformats I.values and I.errors into a single data.frame for plotting.
#'  
#' @param object Object of class \linkS4class{Intercept}.
#' @param iteration Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.
setMethod(
  f = "getPlotValues",
  signature = "Intercept",
  definition = function(object, iteration = NULL) {
    if (is.null(iteration)) {
      values = getValues(object)
      errors = getErrors(object)
    } else {
      if (iteration == dim(object@I.oldValues)[4]) {
        values = getValues(object)
        errors = getErrors(object)
      } else {
        values = getOldValues(object)[,,as.character(iteration)]
        errors = getOldErrors(object)[,,as.character(iteration)]
      }
      
    }
    values = reshape2::melt(values, varnames = c("Strand", "View", "Round"),
                            na.rm = FALSE, as.is = FALSE, value.name = "ddG")
    errors = reshape2::melt(errors, varnames = c("Strand", "View", "Round"),
                            na.rm = FALSE, as.is = FALSE, value.name = "SE")
    
    values = merge(values, errors, sort = FALSE)
    values$Affinity = exp(values$ddG)
    values$SE.Aff = exp(values$ddG)*(values$SE)
    maxAff =  tapply(values$Affinity, values$Round, max, na.rm = TRUE)
    values$AffRef = 0
    values$SERef = 0
    for (l in levels(values$Round)) {
      maxAffInd = c(1:nrow(values))[values$Round == l][(values$Affinity[values$Round == l] == maxAff[l])][!is.na(c(1:nrow(values))[values$Round == l][(values$Affinity[values$Round == l] == maxAff[l])])]
      values$AffRef[values$Round == l] = values$Affinity[maxAffInd[1]]
      values$SERef[values$Round == l] = values$SE.Aff[maxAffInd[1]]
    } 
    values$Affinity = values$Affinity/values$AffRef
    values$SE.Aff = sqrt(values$SE.Aff^2+(values$Affinity/values$AffRef)^2*values$SERef^2)/values$AffRef
    values$AffRef = NULL
    values$SERef = NULL
    values$View = as.numeric(gsub("View.", "", as.character(values$View)))
    return (values)
  }
)



#' Plot \linkS4class{Intercept}
#'  
#' Plots object of class \linkS4class{Intercept} with specific methods for multi-round fitting and view- and strand-specific offsets. 
#' 
#' @param x Object of class \linkS4class{Intercept}.
#' @param sumI Intercept part of design matrix summary.
#' @param includeRound logical: indicating whether or not plots of ddG values for different rounds should be included. If no value is indicated, \linkS4class{Intercept} object will be used to infer whether multiple rounds were fit (includeRound = TRUE) or only one round was used (includeRound = FALSE). 
#' @param includeView logical: indicating whether or not affinities for different views should be plotted. If no value is indicated, \linkS4class{Intercept} object will be used to infer whether individual views were fit (includeView = TRUE) or only one round was used (includeView = FALSE).
#' @param includeDNAStrand logical: indicating whether or not affinities for different DNA strands should be plotted. If no value is indicated, \linkS4class{Intercept} object will be used to infer whether individual strands but not views were fit (includeDNAStrand = TRUE) or not (includeDNAStrand = FALSE).
#' @param rcSymmetric logical: indicating whether reverse complement symmetric version of plotting should be used. If no value is indicated, \linkS4class{Intercept} object will be used to infer whether reverse complement symmetrized values were fit (rcSymmetric = TRUE) or not (rcSymmetric = FALSE).
#' @param title Optional title parameter to be used for plot
#' @param iter Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.
#' @return ggplot Plots of \linkS4class{Intercept} values
#' @export
setMethod(
  f = "plot",
  signature = "Intercept",
  definition = function(x, sumI, includeRound, includeView, includeDNAstrand, rcSymmetric, title=NULL, iter=NULL) {
    values = getValues(x)
    if (missing(includeRound)) {
      if (dim(values)[3] == 1) {
        includeRound = FALSE
      } else {includeRound = TRUE}
    }
    
    
    if ((includeRound == TRUE) & (missing(sumI))) {
      stop('Design must be included to plot intercept affinities with round.')
    }
    
    if (missing(includeView)) {
      if (!all(values[1,,1] == values[1,1,1])) {
        includeView = TRUE
      } else {includeView = FALSE}
    }
    
    if (missing(rcSymmetric)) {
      if (all(values[2,,1] == 0)) {
        rcSymmetric = TRUE
      } else { rcSymmetric = FALSE }
    }
    
    if (missing(includeDNAstrand)) {
      if ((rcSymmetric == FALSE) & (includeView == FALSE)) {
        if (values[1,1,1] != values[2,1,1]) {
          includeDNAstrand = TRUE
        } else { includeDNAstrand = FALSE}
      } else {includeDNAstrand = FALSE}
    }
    
    
    if (includeView == TRUE) {
      if (includeRound == TRUE) {
        p1Title = "View-Specific Affinity Contributions"
        p2Title =  expression(paste('Round-Specific Intercept -', Delta, Delta, 'G values'))
        p1 = plotViews(x, rcSymmetric, title = p1Title)
        p2 = plotRounds(x, title = p2Title)
        g = gridExtra::grid.arrange(p1, p2, nrow = 2,  heights=c(3.0,2.0), 
                                    top=grid::textGrob(title,gp=grid::gpar(fontsize=12,fontface="bold")))
        
      } else {
        p1Title = "View-Specific Affinity Contributions"
        p1 = plotViews(x, rcSymmetric, title = p1Title)
        g = gridExtra::grid.arrange(p1, 
                                    top=grid::textGrob(title,gp=grid::gpar(fontsize=12,fontface="bold")))
      }
    } else {
      if (includeDNAstrand == TRUE) {
        if (includeRound == TRUE) {
          p1Title = "Strand-Specific Affinity Contributions"
          p2Title =  expression(paste('Round-Specific Intercept -', Delta, Delta, 'G values'))
          p1 = plotStrands(x, title = p1Title)
          p2 = plotRounds(x, title = p2Title)
          g = gridExtra::grid.arrange(p1, p2, nrow = 2,  heights=c(3.0,2.0), 
                                      top=grid::textGrob(title,gp=grid::gpar(fontsize=12,fontface="bold")))
        } else {
          p1Title = "Strand-Specific Affinity Contributions"
          p1 = plotStrands(x, title = p1Title)
          g = gridExtra::grid.arrange(p1, 
                                      top=grid::textGrob(title,gp=grid::gpar(fontsize=12,fontface="bold")))
        }
      } else {
        if (includeRound == TRUE) {
          p2Title =  expression(paste('Round-Specific Intercept -', Delta, Delta, 'G values'))
          p2 = plotRounds(x, title = p2Title)
          g = gridExtra::grid.arrange(p2, 
                                      top=grid::textGrob(title,gp=grid::gpar(fontsize=12,fontface="bold")))
        } else {
          p2Title =  expression(paste('Intercept -', Delta, Delta, 'G value'))
          p2 = plotBasicIntercept(x, title = p2Title)
          g = gridExtra::grid.arrange(p2, 
                                      top=grid::textGrob(title,gp=grid::gpar(fontsize=12,fontface="bold")))
        }
        
      } 
    }
    
    return(g)
    
  }
)



#' Plot Strand Offsets for \linkS4class{Intercept}
#' 
#' Plots strand-specific affinity contributions.
#' 
#' @param x Object of class \linkS4class{Intercept}.
#' @param title Title for plot.
#' @param iter Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.
plotStrands = function(x, title=NULL, iter=NULL) {
  values.Plot = getPlotValues(x, iter)
  values.Plot = values.Plot[!is.na(values.Plot$ddG),]
  rowInd = c(1:nrow(values.Plot))[values.Plot$ddG == max(values.Plot$ddG)]
  minErr = min(values.Plot$SE[rowInd])
  roundVar = values.Plot$Round[(values.Plot$ddG == max(values.Plot$ddG)) & (values.Plot$SE == minErr)][1]
  viewValues.Plot = values.Plot[values.Plot$Round == roundVar,]
  viewValues.Plot = viewValues.Plot[order(viewValues.Plot$Strand, viewValues.Plot$View),]
  viewValues.Plot$Strand = as.character(viewValues.Plot$Strand)
  viewValues.Plot$Strand[viewValues.Plot$Strand == "Strand.F"] = "Forward\nStrand"
  viewValues.Plot$Strand[viewValues.Plot$Strand == "Strand.R"] = "Rev. Comp.\nStrand"
  viewValues.Plot$Strand = as.factor(viewValues.Plot$Strand)
  viewValues.Plot = viewValues.Plot[viewValues.Plot$View == 1,]
  yrange = max(viewValues.Plot$ddG+viewValues.Plot$SE)-min(viewValues.Plot$ddG-viewValues.Plot$SE)
  interceptPlot = ggplot2::ggplot(viewValues.Plot)+
    ggplot2::geom_errorbar(data = viewValues.Plot, ggplot2::aes(x = Strand, ymax = ddG+SE, ymin = ddG-SE),col = "navy", width = .25)+
    ggplot2::geom_point(data = viewValues.Plot, ggplot2::aes(x = Strand, y = ddG, shape = Strand), col = "navy", fill = "navy", size = 2)+
    #ggplot2::scale_y_continuous(breaks = c(0, seq(max(viewValues.Plot$Affinity+viewValues.Plot$SE.Aff) %/% .25+1)*.25, .25))+
    ggplot2::ylab(expression(paste('-', Delta, Delta, 'G')))+
    ggplot2::ylim(min(viewValues.Plot$ddG-viewValues.Plot$SE)-yMarg((yrange/20), 2) ,max(viewValues.Plot$ddG+viewValues.Plot$SE)+yMarg((yrange/20), 2))+
    ggplot2::theme_bw()+
    ggplot2::scale_shape_manual(name = NULL,
                                values = c("Forward\nStrand" = 21,
                                           "Rev. Comp.\nStrand" = 23),
                                breaks = c("Forward\nStrand",
                                           "Rev. Comp.\nStrand"),
                                labels = c("Forward\nStrand",
                                           "Rev. Comp.\nStrand"))+
    ggplot2::theme(legend.position="none",
                   strip.text = ggplot2::element_text(face="bold", size = 7),
                   axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                   axis.title.y=ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                   axis.title.x=ggplot2::element_blank(),
                   plot.margin=ggplot2::unit(c(6,18,4,1),"pt"),
                   axis.ticks.length=ggplot2::unit(0.3, "mm"),
                   axis.ticks.x = ggplot2::element_line(size = .5),
                   axis.ticks.y = ggplot2::element_line(size = .5),
                   legend.text = ggplot2::element_text(size=6),
                   axis.text.x = ggplot2::element_text(size = 6,
                                                       vjust = -1),
                   plot.title = ggplot2::element_text(size = 8, 
                                                      face = "bold"),
                   legend.margin=ggplot2::unit(-.05,"cm"), 
                   #legend.position = "bottom",
                   #legend.direction = "vertical", 
                   #legend.box= "horizontal",
                   legend.title = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(
                     fill = "transparent", 
                     size = 1))+
    ggplot2::ggtitle(title)
  
  
  
  
  return(interceptPlot)
}


#' Plot View-specific ddG contributions for object of class \linkS4class{Intercept}
#' 
#' Plots view-specific ddG contributions for object of class \linkS4class{Intercept}. 
#' 
#' If rcSymmetric is TRUE, view-specific beta values are only modeled for one DNA strand orientation (1x numView - 1 beta values).
#' If rcSymmetric is FALSE, view-specific beta values are modeled for both DNA strand orientations (2x numView - 1 beta values).
#' 
#' @param x Object of class \linkS4class{Intercept}.
#' @param rcSymmetric logical: indicates that a reverse complement symmetric model was used for fitting beta values. If missing(rcSymmetric), a value is inferred from the I.values slot of x.
#' @param title Optional title parameter for plot.
#' @param iter Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.
plotViews = function(x, rcSymmetric, title=NULL, iter=NULL) {
  if (!missing(rcSymmetric)) {
    if (rcSymmetric == FALSE) {
      values.Plot = getPlotValues(x, iter)
      values.Plot = values.Plot[!is.na(values.Plot$ddG),]
      rowInd = c(1:nrow(values.Plot))[values.Plot$ddG == max(values.Plot$ddG)]
      minErr = min(values.Plot$SE[rowInd])[1]
      roundVar = values.Plot$Round[(values.Plot$ddG == max(values.Plot$ddG)) & (values.Plot$SE == minErr)]
      viewValues.Plot = values.Plot[values.Plot$Round == roundVar,]
      viewValues.Plot = viewValues.Plot[order(viewValues.Plot$Strand, viewValues.Plot$View),]
      viewValues.Plot$Strand = as.character(viewValues.Plot$Strand)
      viewValues.Plot$Strand[viewValues.Plot$Strand == "Strand.F"] = "Forward\nStrand"
      viewValues.Plot$Strand[viewValues.Plot$Strand == "Strand.R"] = "Rev. Comp.\nStrand"
      viewValues.Plot$Strand = as.factor(viewValues.Plot$Strand)
      yrange = max(viewValues.Plot$ddG+viewValues.Plot$SE)-min(viewValues.Plot$ddG-viewValues.Plot$SE)
      interceptPlot = ggplot2::ggplot(viewValues.Plot)+
        ggplot2::geom_errorbar(data = viewValues.Plot, ggplot2::aes(x = View, ymax = ddG+SE, ymin = ddG-SE),col = "black", width = .25)+
        ggplot2::geom_point(data = viewValues.Plot, ggplot2::aes(x = View, y = ddG, shape = Strand), col = "navy", fill = "navy", size = 2)+
        ggplot2::facet_grid(. ~ Strand)+
        ggplot2::scale_x_discrete(limits = seq(1, ncol(x@I.values), 1),
                                  labels= seq(1, ncol(x@I.values), 1))+
        ggplot2::ylab(expression(paste('-', Delta, Delta, 'G')))+
        ggplot2::coord_cartesian(xlim = c(0, ncol(x@I.values)+1),
                                 ylim = c(min(viewValues.Plot$ddG-viewValues.Plot$SE)-yMarg((yrange/20), 2) ,max(viewValues.Plot$ddG+viewValues.Plot$SE)+yMarg((yrange/20), 2)),
                                 expand = FALSE)+
        ggplot2::theme_bw()+
        ggplot2::scale_shape_manual(name = NULL,
                                    values = c("Forward\nStrand" = 21,
                                               "Rev. Comp.\nStrand" = 23),
                                    breaks = c("Forward\nStrand",
                                               "Rev. Comp.\nStrand"),
                                    labels = c("Forward\nStrand",
                                               "Rev. Comp.\nStrand"))+
        ggplot2::theme(legend.position="none",
                       strip.text = ggplot2::element_text(face="bold", size = 7),
                       axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                       axis.title.y=ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                       axis.title.x=ggplot2::element_text(size = 7, color = "black"),
                       plot.margin=ggplot2::unit(c(4,16,2,1),"pt"),
                       axis.ticks.length=ggplot2::unit(0.3, "mm"),
                       axis.ticks.x = ggplot2::element_line(size = .5),
                       axis.ticks.y = ggplot2::element_line(size = .5),
                       legend.text = ggplot2::element_text(size=6),
                       axis.text.x = ggplot2::element_text(size = 6,
                                                           vjust = -1),
                       plot.title = ggplot2::element_text(size = 8, 
                                                          face = "bold"),
                       legend.margin=ggplot2::unit(-.05,"cm"), 
                       #legend.position = "bottom",
                       #legend.direction = "vertical", 
                       #legend.box= "horizontal",
                       legend.title = ggplot2::element_blank(),
                       legend.background = ggplot2::element_rect(
                         fill = "transparent", 
                         size = 1))+
        
        ggplot2::ggtitle(title)
    } else {
      values.Plot = getPlotValues(x, iter)
      rowInd = c(1:nrow(values.Plot))[!is.na(values.Plot$ddG)][values.Plot$ddG[!is.na(values.Plot$ddG)] == max(values.Plot$ddG, na.rm = TRUE)]
      minErr = min(values.Plot$SE[rowInd])[1]
      roundVar = values.Plot$Round[(values.Plot$ddG[!is.na(values.Plot$ddG)] == max(values.Plot$ddG[!is.na(values.Plot$ddG)])) & (values.Plot$SE[!is.na(values.Plot$ddG)] == minErr)]
      viewValues.Plot = values.Plot[!is.na(values.Plot$ddG),][values.Plot$Round[!is.na(values.Plot$ddG)] == roundVar,]
      viewValues.Plot = viewValues.Plot[order(viewValues.Plot$Strand, viewValues.Plot$View),]
      viewValues.Plot$Strand = as.character(viewValues.Plot$Strand)
      viewValues.Plot = viewValues.Plot[viewValues.Plot$Strand == "Strand.F",]
      viewValues.Plot$Strand[viewValues.Plot$Strand == "Strand.F"] = "Forward\nStrand"
      yrange = max(viewValues.Plot$ddG+viewValues.Plot$SE)-min(viewValues.Plot$ddG-viewValues.Plot$SE)
      interceptPlot = ggplot2::ggplot(viewValues.Plot)+
        ggplot2::geom_errorbar(data = viewValues.Plot, ggplot2::aes(x = View, ymax = ddG+SE, ymin = ddG-SE),col = "black", width = .25)+
        ggplot2::geom_point(data = viewValues.Plot, ggplot2::aes(x = View, y = ddG), col = "navy", fill = "navy", size = 2)+
        ggplot2::scale_x_discrete(limits = seq(1, ncol(x@I.values), 1),
                                  labels= seq(1, ncol(x@I.values), 1))+
        ggplot2::ylab(expression(paste('-', Delta, Delta, 'G')))+
        ggplot2::coord_cartesian(xlim = c(0, ncol(x@I.values)+1),
                                 ylim = c(min(viewValues.Plot$ddG-viewValues.Plot$SE)-yMarg((yrange/20), 2) ,max(viewValues.Plot$ddG+viewValues.Plot$SE)+yMarg((yrange/20), 2)),
                                 expand = FALSE)+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position="none",
                       strip.text = ggplot2::element_text(face="bold", size = 7),
                       axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                       axis.title.y=ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                       axis.title.x=ggplot2::element_text(size = 7, color = "black"),
                       plot.margin=ggplot2::unit(c(6,18,4,1),"pt"),
                       axis.ticks.length=ggplot2::unit(0.3, "mm"),
                       axis.ticks.x = ggplot2::element_line(size = .5),
                       axis.ticks.y = ggplot2::element_line(size = .5),
                       legend.text = ggplot2::element_text(size=6),
                       axis.text.x = ggplot2::element_text(size = 6,
                                                           vjust = -1),
                       plot.title = ggplot2::element_text(size = 8, 
                                                          face = "bold"),
                       legend.margin=ggplot2::unit(-.05,"cm"), 
                       #legend.position = "bottom",
                       #legend.direction = "vertical", 
                       #legend.box= "horizontal",
                       legend.title = ggplot2::element_blank(),
                       legend.background = ggplot2::element_rect(
                         fill = "transparent", 
                         size = 1))+
        ggplot2::ggtitle(title)
    }
  } else {
    values.Plot = getPlotValues(x, iter)
    rowInd = c(1:nrow(values.Plot))[values.Plot$ddG == max(values.Plot$ddG)]
    minErr = min(values.Plot$SE[rowInd])[1]
    roundVar = values.Plot$Round[(values.Plot$ddG == max(values.Plot$ddG)) & (values.Plot$SE == minErr)]
    viewValues.Plot = values.Plot[values.Plot$Round == roundVar,]
    viewValues.Plot = viewValues.Plot[order(viewValues.Plot$Strand, viewValues.Plot$View),]
    viewValues.Plot$Strand = as.character(viewValues.Plot$Strand)
    viewValues.Plot$Strand[viewValues.Plot$Strand == "Strand.F"] = "Forward\nStrand"
    viewValues.Plot$Strand[viewValues.Plot$Strand == "Strand.R"] = "Rev. Comp.\nStrand"
    viewValues.Plot$Strand = as.factor(viewValues.Plot$Strand)
    yrange = max(viewValues.Plot$ddG+viewValues.Plot$SE)-min(viewValues.Plot$ddG-viewValues.Plot$SE)
    interceptPlot = ggplot2::ggplot(viewValues.Plot)+
      ggplot2::geom_errorbar(data = viewValues.Plot, ggplot2::aes(x = View, ymax = ddG+SE, ymin = ddG-SE),col = "black", width = .25)+
      ggplot2::geom_point(data = viewValues.Plot, ggplot2::aes(x = View, y = ddG, shape = Strand), col = "navy", fill = "navy", size = 2)+
      ggplot2::facet_grid(. ~ Strand)+
      ggplot2::scale_x_discrete(limits = seq(1, max(viewValues.Plot$View), 1),
                                labels= seq(1, max(viewValues.Plot$View), 1))+
      ggplot2::ylab(expression(paste('-', Delta, Delta, 'G')))+
      ggplot2::coord_cartesian(xlim = c(0, max(viewValues.Plot$View)+1),
                               ylim = c(min(viewValues.Plot$ddG-viewValues.Plot$SE)-yMarg((yrange/20), 2) ,max(viewValues.Plot$ddG+viewValues.Plot$SE)+yMarg((yrange/20), 2)),
                               expand = FALSE)+
      ggplot2::theme_bw()+
      ggplot2::scale_shape_manual(name = NULL,
                                  values = c("Forward\nStrand" = 21,
                                             "Rev. Comp.\nStrand" = 23),
                                  breaks = c("Forward\nStrand",
                                             "Rev. Comp.\nStrand"),
                                  labels = c("Forward\nStrand",
                                             "Rev. Comp.\nStrand"))+
      ggplot2::theme(legend.position="none",
                     strip.text = ggplot2::element_text(face="bold", size = 7),
                     axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                     axis.title.y=ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                     axis.title.x=ggplot2::element_text(size = 6, color = "black"),
                     plot.margin=ggplot2::unit(c(4,16,2,1),"pt"),
                     axis.ticks.length=ggplot2::unit(0.3, "mm"),
                     axis.ticks.x = ggplot2::element_line(size = .5),
                     axis.ticks.y = ggplot2::element_line(size = .5),
                     legend.text = ggplot2::element_text(size=6),
                     axis.text.x = ggplot2::element_text(size = 6,
                                                         vjust = -1),
                     plot.title = ggplot2::element_text(size = 8, 
                                                        face = "bold"),
                     legend.margin=ggplot2::unit(-.05,"cm"), 
                     
                     legend.title = ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(
                       fill = "transparent", 
                       size = 1))+
      
      ggplot2::ggtitle(title)
  }
  
  
  
  
  
  return(interceptPlot)
}


#' Plot round-specific -ddG values for object of class \linkS4class{Intercept}
#' 
#' Plots round-specific -ddG values for object of class \linkS4class{Intercept}. 
#' 
#' @param x Object of class \linkS4class{Intercept}.
#' @param title Optional title parameter for plot.
#' @param iter Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.
plotRounds = function(x, sumI, title=NULL, iter=NULL) {
  values.Plot = getPlotValues(x, iter)
  values.Plot = values.Plot[!is.na(values.Plot$ddG),]
  sumI = sumI[,paste("View.", 1:ncol(x@I.values), sep = "")]
  maxViewIndex = which(sumI == max(sumI), arr.ind = TRUE)
  refStrand = rownames(sumI)[maxViewIndex[[1]]]
  refView = maxViewIndex[[2]]
  values.Plot = values.Plot[(values.Plot$Strand == refStrand) & (values.Plot$View == refView),]
  maxDdg =  which.max(values.Plot$ddG)
  values.Plot$ddG = values.Plot$ddG-values.Plot$ddG[maxDdg]
  values.Plot$SE = sqrt(values.Plot$SE^2+values.Plot$SE[maxDdg]^2)
  values.Plot$Round = as.character(values.Plot$Round)
  values.Plot$Round = as.factor(gsub("Round.", "", values.Plot$Round))
  yrange = max(values.Plot$ddG+values.Plot$SE)-min(values.Plot$ddG-values.Plot$SE)
  roundPlot = ggplot2::ggplot(values.Plot)+
    ggplot2::geom_errorbar(data = values.Plot, ggplot2::aes(x = Round, ymax = ddG+SE, ymin = ddG-SE),col = "black", width = .05)+
    ggplot2::geom_point(data = values.Plot, ggplot2::aes(x = Round, y = ddG), col = "#330000", fill = "#330000", size = 2, shape = 21)+
    #ggplot2::scale_y_continuous(name = , seq((min(values.Plot$ddG-values.Plot$SE) %/% .25)*.25, (max(values.Plot$ddG+values.Plot$SE) %/% .25+1)*.25, .25))+
    ggplot2::theme_bw()+
    ggplot2::xlab("Round")+
    ggplot2::ylab(expression(paste('-', Delta, Delta, 'G')))+
    ggplot2::ylim(min(values.Plot$ddG-values.Plot$SE)-yMarg((yrange/20), 2) ,max(values.Plot$ddG+values.Plot$SE)+yMarg((yrange/20), 2))+
    ggplot2::theme(legend.position="none",
                   strip.text = ggplot2::element_text(face="bold", size = 7),
                   axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                   axis.title.y=ggplot2::element_text(size = 6, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                   axis.title.x=ggplot2::element_text(size = 6, color = "black"),
                   plot.margin=ggplot2::unit(c(6,18,4,1),"pt"),
                   axis.ticks.length=ggplot2::unit(0.3, "mm"),
                   axis.ticks.x = ggplot2::element_line(size = .5),
                   axis.ticks.y = ggplot2::element_line(size = .5),
                   legend.text = ggplot2::element_text(size=6),
                   axis.text.x = ggplot2::element_text(size = 6,
                                                       vjust = -1),
                   plot.title = ggplot2::element_text(size = 8, 
                                                      face = "bold"),
                   legend.margin=ggplot2::unit(-.05,"cm"), 
                   #legend.position = "bottom",
                   #legend.direction = "vertical", 
                   #legend.box= "horizontal",
                   legend.title = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(
                     fill = "transparent", 
                     size = 1))+
    ggplot2::ggtitle(title)
  
  
  return(roundPlot)
}


#' Plot simple intercept -ddG value for object of class \linkS4class{Intercept}
#' 
#' Plots simple intercept -ddG value for object of class \linkS4class{Intercept}. 
#' 
#' @param x Object of class \linkS4class{Intercept}.
#' @param title Optional title parameter for plot.
#' @param iter Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.
plotBasicIntercept = function(x, title=NULL, iter=NULL) {
  values.Plot = getPlotValues(x, iter)
  values.Plot = values.Plot[!is.na(values.Plot$ddG),]
  values.Plot = values.Plot[(values.Plot$Strand == "Strand.F"),]
  values.Plot = values.Plot[(values.Plot$View == min(values.Plot$View)),]
  values.Plot$Round = as.character(values.Plot$Round)
  values.Plot$Round = as.factor(gsub("Round.", "", values.Plot$Round))
  yrange = max(values.Plot$ddG+values.Plot$SE)-min(values.Plot$ddG-values.Plot$SE)
  roundPlot = ggplot2::ggplot(values.Plot)+
    ggplot2::geom_errorbar(data = values.Plot, ggplot2::aes(x = Round, ymax = ddG+SE, ymin = ddG-SE),col = "black", width = .05)+
    ggplot2::geom_point(data = values.Plot, ggplot2::aes(x = Round, y = ddG), col = "#330000", fill = "#330000", size = 2, shape = 21)+
    ggplot2::theme_bw()+
    ggplot2::ylab(expression(paste('-', Delta, Delta, 'G')))+
    ggplot2::ylim(min(values.Plot$ddG-values.Plot$SE)-yMarg((yrange/20), 2) ,max(values.Plot$ddG+values.Plot$SE)+yMarg((yrange/20), 2))+
    ggplot2::theme(legend.position="none",
                   strip.text = ggplot2::element_text(face="bold", size = 7),
                   axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                   axis.title.y=ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                   axis.title.x=ggplot2::element_blank(),
                   plot.margin=ggplot2::unit(c(6,18,18,1),"pt"),
                   axis.ticks.length=ggplot2::unit(0.3, "mm"),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_line(size = .5),
                   legend.text = ggplot2::element_text(size=6),
                   axis.text.x = ggplot2::element_text(size = 6,
                                                       vjust = -1),
                   plot.title = ggplot2::element_text(size = 8, 
                                                      face = "bold"),
                   legend.margin=ggplot2::unit(-.05,"cm"), 
                   #legend.position = "bottom",
                   #legend.direction = "vertical", 
                   #legend.box= "horizontal",
                   legend.title = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(
                     fill = "transparent", 
                     size = 1))+
    ggplot2::ggtitle(title)
  
  
  
  return(roundPlot)
}

#' Values for \linkS4class{Intercept}
#' 
#' Returns I.values slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "getValues",
  signature = "Intercept",
  definition = function(object) {
    return(object@I.values)
  }
)

#' Errors for class \linkS4class{Intercept}
#' 
#' Returns I.errors slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "getErrors",
  signature = "Intercept",
  definition = function(object) {
    return(object@I.errors)
  }
)

#' Z-scores for class \linkS4class{Intercept}
#' 
#' Returns I.z slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "getZ",
  signature = "Intercept",
  definition = function(object) {
    return(object@I.z)
  }
)

#' P-values for class \linkS4class{Intercept}
#' 
#' Returns I.sig slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "getSig",
  signature = "Intercept",
  definition = function(object) {
    return(object@I.sig)
  }
)

#' Values for previous iterations for \linkS4class{Intercept}
#' 
#' Returns I.oldValues slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "getOldValues",
  signature = "Intercept",
  definition = function(object) {
    return(object@I.oldValues)
  }
)

#' Errors for Previous Iterations for \linkS4class{Intercept}
#' 
#' Returns I.oldErrors slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "getOldErrors",
  signature = "Intercept",
  definition = function(object) {
    return(object@I.oldErrors)
  }
)


#' Z-scores for Previous Iterations for \linkS4class{Intercept}
#' 
#' Returns I.oldZ slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export
setMethod(
  f = "getOldZ",
  signature = "Intercept",
  definition = function(object) {
    return(object@I.oldZ)
  }
)

#' P-values for Previous Iterations for class \linkS4class{Intercept}
#' 
#' Returns I.oldSig slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class '\linkS4class{Intercept}.
#' @export
setMethod(
  f = "getOldSig",
  signature = "Intercept",
  definition = function(object) {
    return(object@I.oldSig)
  }
)

#' Set Values for class \linkS4class{Intercept}
#' 
#' Sets I.values slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value 3-D array with dimension of 2 x (numViews) x (number of rounds) containing -ddG values for 'Intercept' offsets.
#' @export
setReplaceMethod(
  f = "setValues",
  signature = "Intercept",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@I.values)) {
      stop('Current and replacement versions of I.values must have the same dimensions.')
    }
    if (dim(value)[3] != dim(object@I.values)[3]) {
      stop('Current and replacement versions of I.values must have the same dimensions.')
    }
    if (length(dimnames(value)[[3]]) == 0) {
      stop('New I.values must have rounds labeled.')
    } else if (length(dimnames(value)[[3]]) != length(dimnames(object@I.values)[[3]])) {
      stop('Dimension names for new I.values must have the same length as dimension names for the old I.values.')
    } else if (!all(unique(dimnames(value)[[3]])[order(unique(dimnames(value)[[3]]))] == unique(dimnames(object@I.values)[[3]])[order(unique(dimnames(object@I.values)[[3]]))])) {
      stop('Dimension names for new I.values must be the same as dimension names for the old I.values.')
    }
    value = value[,,dimnames(object@I.values)[[3]],drop=FALSE]
    dimnames(value) = list(c("Strand.F", "Strand.R"),
                           c(1:ncol(object@I.values)),
                           dimnames(object@I.values)[[3]])
    object@I.values <- value
    validObject(object)
    return(object)
  }
)

#' Set Errors for \linkS4class{Intercept}
#' 
#' Sets I.errors slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value 3-D array with dimension of 2 x (numViews) x (number of rounds) containing -ddG errors for 'Intercept' offsets.
#' @export
setReplaceMethod(
  f = "setErrors",
  signature = "Intercept",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@I.errors)) {
      stop('Current and replacement versions of I.errors must have the same dimensions.')
    }
    if (dim(value)[3] != dim(object@I.errors)[3]) {
      stop('Current and replacement versions of I.errors must have the same dimensions.')
    }
    if (length(dimnames(value)[[3]]) == 0) {
      stop('New I.errors must have rounds labeled.')
    } else if (length(dimnames(value)[[3]]) != length(dimnames(object@I.errors)[[3]])) {
      stop('Dimension names for new I.errors must have the same length as dimension names for the old I.errors.')
    } else if (!all(unique(dimnames(value)[[3]])[order(unique(dimnames(value)[[3]]))] == unique(dimnames(object@I.errors)[[3]])[order(unique(dimnames(object@I.errors)[[3]]))])) {
      stop('Dimension names for new I.errors must be the same as dimension names for the old I.errors.')
    }
    value = value[,,dimnames(object@I.errors)[[3]],drop=FALSE]
    dimnames(value) = list(c("Strand.F", "Strand.R"),
                           c(1:ncol(object@I.errors)),
                           dimnames(object@I.errors)[[3]])
    object@I.errors <- value
    validObject(object)
    return(object)
  }
)

#' Set Z-scores for \linkS4class{Intercept}
#' 
#' Sets I.z slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value 3-D array with dimension of 2 x (numViews) x (number of rounds) containing -ddG z-scores for 'Intercept' offsets.
#' @export
setReplaceMethod(
  f = "setZ",
  signature = "Intercept",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@I.z)) {
      stop('Current and replacement versions of I.z must have the same dimensions.')
    }
    if (dim(value)[3] != dim(object@I.z)[3]) {
      stop('Current and replacement versions of I.z must have the same dimensions.')
    }
    if (length(dimnames(value)[[3]]) == 0) {
      stop('New I.z must have rounds labeled.')
    } else if (length(dimnames(value)[[3]]) != length(dimnames(object@I.z)[[3]])) {
      stop('Dimension names for new I.z must have the same length as dimension names for the old I.z.')
    } else if (!all(unique(dimnames(value)[[3]])[order(unique(dimnames(value)[[3]]))] == unique(dimnames(object@I.z)[[3]])[order(unique(dimnames(object@I.z)[[3]]))])) {
      stop('Dimension names for new I.z must be the same as dimension names for the old I.z.')
    }
    value = value[,,dimnames(object@I.z)[[3]],drop=FALSE]
    dimnames(value) = list(c("Strand.F", "Strand.R"),
                           c(1:ncol(object@I.z)),
                           dimnames(object@I.z)[[3]])
    object@I.z <- value
    validObject(object)
    return(object)
  }
)

#' Set P-values for \linkS4class{Intercept}
#' 
#' Sets I.sig slot from object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value 3-D array with dimension of 2 x (numViews) x (number of rounds) containing -ddG p-values for 'Intercept' offsets.
#' @export
setReplaceMethod(
  f = "setSig",
  signature = "Intercept",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@I.sig)) {
      stop('Current and replacement versions of I.sig must have the same dimensions.')
    }
    if (dim(value)[3] != dim(object@I.sig)[3]) {
      stop('Current and replacement versions of I.sig must have the same dimensions.')
    }
    if (length(dimnames(value)[[3]]) == 0) {
      stop('New I.sig must have rounds labeled.')
    } else if (length(dimnames(value)[[3]]) != length(dimnames(object@I.sig)[[3]])) {
      stop('Dimension names for new I.sig must have the same length as dimension names for the old I.sig.')
    } else if (!all(unique(dimnames(value)[[3]])[order(unique(dimnames(value)[[3]]))] == unique(dimnames(object@I.sig)[[3]])[order(unique(dimnames(object@I.sig)[[3]]))])) {
      stop('Dimension names for new I.sig must be the same as dimension names for the old I.sig.')
    }
    value = value[,,dimnames(object@I.sig)[[3]],drop=FALSE]
    dimnames(value) = list(c("Strand.F", "Strand.R"),
                           c(1:ncol(object@I.sig)),
                           dimnames(object@I.sig)[[3]])
    object@I.sig <- value
    validObject(object)
    return(object)
  }
)

#' Set Values for Previous Iterations for \linkS4class{Intercept}
#' 
#' Sets I.oldValues slot for object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value 4-D array with dimension of 2 x (numViews) x (number of rounds) x (number of glm iterations) containing -ddG values for 'Intercept' offsets.
#' @export
setReplaceMethod(
  f = "setOldValues",
  signature = "Intercept",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@I.oldValues)) {
      stop('Current and replacement versions of I.oldValues must have the same dimensions.')
    }
    if (dim(value)[3] != dim(object@I.oldValues)[3]) {
      stop('Current and replacement versions of I.oldValues must have the same dimensions.')
    }
    if (dim(value)[4] != dim(object@I.oldValues)[4]) {
      stop('Current and replacement versions of I.oldValues must have the same dimensions.')
    }
    if (length(dimnames(value)[[3]]) == 0) {
      stop('New I.oldValues must have rounds labeled.')
    } else if (length(dimnames(value)[[3]]) != length(dimnames(object@I.oldValues)[[3]])) {
      stop('Dimension names for new I.oldValues must have the same length as dimension names for the old I.oldValues.')
    } else if (!all(unique(dimnames(value)[[3]])[order(unique(dimnames(value)[[3]]))] == unique(dimnames(object@I.oldValues)[[3]])[order(unique(dimnames(object@I.oldValues)[[3]]))])) {
      stop('Dimension names for new I.oldValues must be the same as dimension names for the old I.oldValues.')
    }
    value = value[,,dimnames(object@I.oldValues)[[3]],,drop=FALSE]
    
    if (dim(value)[4] > 0) {
      dimnames(value) = list(c("Strand.F", "Strand.R"),
                             c(1:ncol(object@I.oldValues)),
                             dimnames(object@I.oldValues)[[3]],
                             c((dim(value)[4]-1):0))
    } else {
      dimnames(value) = list(c("Strand.F", "Strand.R"),
                             c(1:ncol(object@I.oldValues)),
                             dimnames(object@I.oldValues)[[3]],
                             NULL)
    }
    
    object@I.oldValues <- value
    validObject(object)
    return(object)
  }
)

#' Set Errors for Previous Iterations for \linkS4class{Intercept}
#' 
#' Sets I.oldErrors slot for object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value 4-D array with dimension of 2 x (numViews) x (number of rounds) x (number of glm iterations) containing -ddG errors for 'Intercept' offsets.
#' @export
setReplaceMethod(
  f = "setOldErrors",
  signature = "Intercept",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@I.oldErrors)) {
      stop('Current and replacement versions of I.oldErrors must have the same dimensions.')
    }
    if (dim(value)[3] != dim(object@I.oldErrors)[3]) {
      stop('Current and replacement versions of I.oldErrors must have the same dimensions.')
    }
    if (dim(value)[4] != dim(object@I.oldErrors)[4]) {
      stop('Current and replacement versions of I.oldErrors must have the same dimensions.')
    }
    if (length(dimnames(value)[[3]]) == 0) {
      stop('New I.oldErrors must have rounds labeled.')
    } else if (length(dimnames(value)[[3]]) != length(dimnames(object@I.oldErrors)[[3]])) {
      stop('Dimension names for new I.oldErrors must have the same length as dimension names for the old I.oldErrors.')
    } else if (!all(unique(dimnames(value)[[3]])[order(unique(dimnames(value)[[3]]))] == unique(dimnames(object@I.oldErrors)[[3]])[order(unique(dimnames(object@I.oldErrors)[[3]]))])) {
      stop('Dimension names for new I.oldErrors must be the same as dimension names for the old I.oldErrors.')
    }
    value = value[,,dimnames(object@I.oldErrors)[[3]],,drop=FALSE]
    
    if (dim(value)[4] > 0) {
      dimnames(value) = list(c("Strand.F", "Strand.R"),
                             c(1:ncol(object@I.oldErrors)),
                             dimnames(object@I.oldErrors)[[3]],
                             c((dim(value)[4]-1):0))
    } else {
      dimnames(value) = list(c("Strand.F", "Strand.R"),
                             c(1:ncol(object@I.oldErrors)),
                             dimnames(object@I.oldErrors)[[3]],
                             NULL)
    }
    
    object@I.oldErrors <- value
    validObject(object)
    return(object)
  }
)

#' Set Z-scores for Previous Iterations for \linkS4class{Intercept}
#' 
#' Sets I.oldZ slot for object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value 4-D array with dimension of 2 x (numViews) x (number of rounds) x (number of glm iterations) containing -ddG z-scores for 'Intercept' offsets.
#' @export
setReplaceMethod(
  f = "setOldZ",
  signature = "Intercept",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@I.oldZ)) {
      stop('Current and replacement versions of I.oldZ must have the same dimensions.')
    }
    if (dim(value)[3] != dim(object@I.oldZ)[3]) {
      stop('Current and replacement versions of I.oldZ must have the same dimensions.')
    }
    if (dim(value)[4] != dim(object@I.oldZ)[4]) {
      stop('Current and replacement versions of I.oldZ must have the same dimensions.')
    }
    if (length(dimnames(value)[[3]]) == 0) {
      stop('New I.oldZ must have rounds labeled.')
    } else if (length(dimnames(value)[[3]]) != length(dimnames(object@I.oldZ)[[3]])) {
      stop('Dimension names for new I.oldZ must have the same length as dimension names for the old I.oldZ.')
    } else if (!all(unique(dimnames(value)[[3]])[order(unique(dimnames(value)[[3]]))] == unique(dimnames(object@I.oldZ)[[3]])[order(unique(dimnames(object@I.oldZ)[[3]]))])) {
      stop('Dimension names for new I.oldZ must be the same as dimension names for the old I.oldZ.')
    }
    value = value[,,dimnames(object@I.oldZ)[[3]],,drop=FALSE]
    
    if (dim(value)[4] > 0) {
      dimnames(value) = list(c("Strand.F", "Strand.R"),
                             c(1:ncol(object@I.oldZ)),
                             dimnames(object@I.oldZ)[[3]],
                             c((dim(value)[4]-1):0))
    } else {
      dimnames(value) = list(c("Strand.F", "Strand.R"),
                             c(1:ncol(object@I.oldZ)),
                             dimnames(object@I.oldZ)[[3]],
                             NULL)
    }
    object@I.oldZ <- value
    validObject(object)
    return(object)
  }
)

#' Set P-values for Previous Iterations for \linkS4class{Intercept}
#' 
#' Sets I.oldSig slot for object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value 4-D array with dimension of 2 x (numViews) x (number of rounds) x (number of glm iterations) containing -ddG p-values for 'Intercept' offsets.
#' @export
setReplaceMethod(
  f = "setOldSig",
  signature = "Intercept",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@I.oldSig)) {
      stop('Current and replacement versions of I.oldSig must have the same dimensions.')
    }
    if (dim(value)[3] != dim(object@I.oldSig)[3]) {
      stop('Current and replacement versions of I.oldSig must have the same dimensions.')
    }
    if (dim(value)[4] != dim(object@I.oldSig)[4]) {
      stop('Current and replacement versions of I.oldSig must have the same dimensions.')
    }
    if (length(dimnames(value)[[3]]) == 0) {
      stop('New I.oldSig must have rounds labeled.')
    } else if (length(dimnames(value)[[3]]) != length(dimnames(object@I.oldSig)[[3]])) {
      stop('Dimension names for new I.oldSig must have the same length as dimension names for the old I.oldSig.')
    } else if (!all(unique(dimnames(value)[[3]])[order(unique(dimnames(value)[[3]]))] == unique(dimnames(object@I.oldSig)[[3]])[order(unique(dimnames(object@I.oldSig)[[3]]))])) {
      stop('Dimension names for new I.oldSig must be the same as dimension names for the old I.oldSig.')
    }
    value = value[,,dimnames(object@I.oldSig)[[3]],,drop=FALSE]
    
    if (dim(value)[4] > 0) {
      dimnames(value) = list(c("Strand.F", "Strand.R"),
                             c(1:ncol(object@I.oldSig)),
                             dimnames(object@I.oldSig)[[3]],
                             c((dim(value)[4]-1):0))
    } else {
      dimnames(value) = list(c("Strand.F", "Strand.R"),
                             c(1:ncol(object@I.oldSig)),
                             dimnames(object@I.oldSig)[[3]],
                             NULL)
    }
    object@I.oldSig <- value
    validObject(object)
    return(object)
  }
)

#' Update Values for \linkS4class{Intercept}
#' 
#' Updates I.values and I.oldValues slots for object of class \linkS4class{Intercept} using output of glm model fit.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value Glm fit object.
setReplaceMethod(
  f = "updateValues",
  signature = "Intercept",
  definition = function(object, value){
    oldValues = getOldValues(object)
    new4thDim = dim(oldValues)[4]+1
    oV.updated = array(0, dim = c(dim(oldValues)[1:3], new4thDim))
    if (dim(oldValues)[4] > 0) {
      oV.updated[,,,2:new4thDim] = oldValues
    }
    oV.updated[,,,1] = getValues(object)
    dimnames(oV.updated) = list(dimnames(oldValues)[[1]],
                                dimnames(oldValues)[[2]],
                                dimnames(oldValues)[[3]],
                                c((new4thDim-1):0))
    object@I.oldValues <- oV.updated
    dimnames(value) = list(dimnames(oldValues)[[1]],
                           dimnames(oldValues)[[2]],
                           dimnames(oldValues)[[3]])
    object@I.values <- value
    return (object)
  }
)

#' Update Errors for \linkS4class{Intercept}
#' 
#' Updates I.errors and I.oldErrors slots for object of class \linkS4class{Intercept} using output of glm model fit.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value Glm fit object. 
setReplaceMethod(
  f = "updateErrors",
  signature = "Intercept",
  definition = function(object, value){
    oldErrors = getOldErrors(object)
    new4thDim = dim(oldErrors)[4]+1
    oV.updated = array(0, dim = c(dim(oldErrors)[1:3], new4thDim))
    if (dim(oldErrors)[4] > 0) {
      oV.updated[,,,2:new4thDim] = oldErrors
    }
    oV.updated[,,,1] = getErrors(object)
    dimnames(oV.updated) = list(dimnames(oldErrors)[[1]],
                                dimnames(oldErrors)[[2]],
                                dimnames(oldErrors)[[3]],
                                c((new4thDim-1):0))
    object@I.oldErrors <- oV.updated
    dimnames(value) = list(dimnames(oldErrors)[[1]],
                           dimnames(oldErrors)[[2]],
                           dimnames(oldErrors)[[3]])
    object@I.errors <- value
    return (object)
  }
)

#' Update Z-scores for \linkS4class{Intercept}
#' 
#' Updates I.z and I.oldZ slots for object of class \linkS4class{Intercept} using output of glm model fit.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value Glm fit object.
setReplaceMethod(
  f = "updateZ",
  signature = "Intercept",
  definition = function(object, value){
    oldZ = getOldZ(object)
    new4thDim = dim(oldZ)[4]+1
    oV.updated = array(0, dim = c(dim(oldZ)[1:3], new4thDim))
    if (dim(oldZ)[4] > 0) {
      oV.updated[,,,2:new4thDim] = oldZ
    }
    oV.updated[,,,1] = getZ(object)
    dimnames(oV.updated) = list(dimnames(oldZ)[[1]],
                                dimnames(oldZ)[[2]],
                                dimnames(oldZ)[[3]],
                                c((new4thDim-1):0))
    object@I.oldZ <- oV.updated
    dimnames(value) = list(dimnames(oldZ)[[1]],
                           dimnames(oldZ)[[2]],
                           dimnames(oldZ)[[3]])
    object@I.z <- value
    return (object)
  }
)

#' Update P-values for \linkS4class{Intercept}
#' 
#' Updates I.sig and I.oldSig slots for object of class \linkS4class{Intercept} using output of glm model fit.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param value Glm fit object.
setReplaceMethod(
  f = "updateSig",
  signature = "Intercept",
  definition = function(object, value){
    oldSig = getOldSig(object)
    new4thDim = dim(oldSig)[4]+1
    oV.updated = array(0, dim = c(dim(oldSig)[1:3], new4thDim))
    if (dim(oldSig)[4] > 0) {
      oV.updated[,,,2:new4thDim] = oldSig
    }
    oV.updated[,,,1] = getSig(object)
    dimnames(oV.updated) = list(dimnames(oldSig)[[1]],
                                dimnames(oldSig)[[2]],
                                dimnames(oldSig)[[3]],
                                c((new4thDim-1):0))
    object@I.oldSig <- oV.updated
    dimnames(value) = list(dimnames(oldSig)[[1]],
                           dimnames(oldSig)[[2]],
                           dimnames(oldSig)[[3]])
    object@I.sig <- value
    return (object)
  }
)


#' Add New Betas for \linkS4class{Intercept}
#' 
#' Updates all slots for object of class \linkS4class{Intercept} using glm model fit.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param design Data frame passed to glm that contains the feature design as output by the function addDesignMatrix.
#' @param value Glm fit object.
#' @export
setMethod(
  f = "addNewBetas",
  signature = "Intercept",
  definition = function(object, design, value, missingValueSuppression, includeDNAstrand, includeView, rcSymmetric){
    betaSummary = summary(value)$coefficients
    interceptData = betaSummary[grep("Intercept", rownames(betaSummary)),]
    newValues = array(interceptData[1],
                      dim = dim(getValues(object)),
                      dimnames = dimnames(getValues(object)))
    newErrors = array(interceptData[2],
                      dim = dim(getValues(object)),
                      dimnames = dimnames(getValues(object)))
    newZ = array(interceptData[3],
                 dim = dim(getValues(object)),
                 dimnames = dimnames(getValues(object)))
    newSig = array(interceptData[4],
                   dim = dim(getValues(object)),
                   dimnames = dimnames(getValues(object)))
    if (rcSymmetric == TRUE) {
      newValues["Strand.R",,] = 0
      newErrors["Strand.R",,] = 0
      newZ["Strand.R",,] = 0
      newSig["Strand.R",,] = 0
    }
    grStr = paste("Round.", "[0-9]{1,}", sep = "")
    round.inds = grep(grStr, rownames(betaSummary))
    betaSub = betaSummary[round.inds,,drop = FALSE]
    if (nrow(betaSub) > 0) {
      if (rcSymmetric == TRUE) {
        newValues["Strand.F",,rownames(betaSub)] = newValues["Strand.F",,rownames(betaSub)]+array(betaSub[rownames(betaSub), 1], 
                                                                                                  dim = c(1, dim(newValues)[2]))
        newErrors["Strand.F",,rownames(betaSub)] = sqrt(newErrors["Strand.F",,rownames(betaSub)]^2+array(betaSub[rownames(betaSub), 2], 
                                                                                                         dim = c(1, dim(newValues)[2]))^2)
        newZ["Strand.F",,rownames(betaSub)] = newValues["Strand.F",,rownames(betaSub)]/newErrors["Strand.F",,rownames(betaSub)]
        newSig["Strand.F",,rownames(betaSub)] = 2*pnorm(-abs(newZ["Strand.F",,rownames(betaSub)]))
      } else {
        newValues[,,rownames(betaSub)] = newValues[,,rownames(betaSub)]+array(betaSub[rownames(betaSub), 1], 
                                                                                                  dim = c(2, dim(newValues)[2]))
        newErrors[,,rownames(betaSub)] = sqrt(newErrors[,,rownames(betaSub)]^2+array(betaSub[rownames(betaSub), 2], 
                                                                                                         dim = c(2, dim(newValues)[2]))^2)
        newZ[,,rownames(betaSub)] = newValues[,,rownames(betaSub)]/newErrors[,,rownames(betaSub)]
        newSig[,,rownames(betaSub)] = 2*pnorm(-abs(newZ[,,rownames(betaSub)]))
      }
    }
    if (includeView == TRUE) {
      if (rcSymmetric == FALSE) {
        grView = paste("Strand.[FR][0-9]{1,}", sep = "")
      } else {
        grView = paste("Strand.F[0-9]{1,}", sep = "")
      }

      view.inds = grep(grView, rownames(betaSummary))
      V.values = array(NA, dim = dim(newValues),
                       dimnames = dimnames(newValues))
      V.errors = array(NA, dim = dim(newValues),
                       dimnames = dimnames(newValues))
      V.z = array(NA, dim = dim(newValues),
                  dimnames = dimnames(newValues))
      V.sig = array(NA, dim = dim(newValues),
                    dimnames = dimnames(newValues))
      if (rcSymmetric == TRUE) {
        V.values["Strand.R",,] = 0
        V.errors["Strand.R",,] = 0
        V.z["Strand.R",,] = 0
        V.sig["Strand.R",,] = 0
      }
      sumI = getDesignMatrix(object, design)
      maxViewIndex = which(sumI == max(sumI), arr.ind = TRUE)[1,]
      V.values[maxViewIndex[[1]], maxViewIndex[[2]],] = 0
      V.errors[maxViewIndex[[1]], maxViewIndex[[2]],] = 0
      V.z[maxViewIndex[[1]], maxViewIndex[[2]],] = 0
      V.sig[maxViewIndex[[1]], maxViewIndex[[2]],] = 0
      Vs = rownames(betaSummary)[grep(grView, rownames(betaSummary))]
      if (length(Vs) > 0) {
        for (vVar in Vs) {
          vStr = gsub("[0-9]{1,}", "", vVar)
          vV = paste("View.", gsub("Strand.[FR]", "", vVar), sep = "")
          V.values[vStr,vV,] = betaSummary[vVar,1]
          V.errors[vStr,vV,] = betaSummary[vVar,2]
          V.z[vStr,vV,] = betaSummary[vVar,3]
          V.sig[vStr,vV,] = betaSummary[vVar,4]
        } 
      }
      
      badV = which(is.na(V.values[,,1]), arr.ind = TRUE)
      if (length(badV) > 0) {
        if (rcSymmetric == FALSE) {
          minV = which(V.values[,,1] == min(V.values[,,1], na.rm = TRUE), arr.ind = TRUE)
          V.values[badV,,] = c(V.values[minV[[1]], minV[[2]], ]-missingValueSuppression)[1]
          V.errors[badV,,] = c(V.errors[minV[[1]], minV[[2]], ]+missingValueSuppression)[1]
          V.z[badV,,] = 0
          V.sig[badV,,] = 1
        } else {
          minV = which(V.values[,,1] == min(V.values[1,,1], na.rm = TRUE), arr.ind = TRUE)
          V.values[badV,,] = c(V.values[1, minV[[1]], ]-missingValueSuppression)[1]
          V.errors[badV,,] = c(V.errors[1, minV[[1]], ]+missingValueSuppression)[1]
          V.z[badV,,] = 0
          V.sig[badV,,] = 1
        }
        
      }
      
      
      newValues[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1] = newValues[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]+V.values[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]
      newErrors[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1] = sqrt(newErrors[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]^2+V.errors[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]^2)
      newZ[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1] = newValues[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]/newErrors[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]
      newSig[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1] = 2*pnorm(-abs(newZ[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]))
      
      
    } else if (includeDNAstrand == TRUE) {
      St.set = c("Strand.F", "Strand.R")
      St.values = array(NA, dim = dim(newValues),
                        dimnames = dimnames(newValues))
      St.errors = array(NA, dim = dim(newValues),
                        dimnames = dimnames(newValues))
      St.z = array(NA, dim = dim(newValues),
                   dimnames = dimnames(newValues))
      St.sig = array(NA, dim = dim(newValues),
                     dimnames = dimnames(newValues))
      sumI = getDesignMatrix(object, design)
      sumS = matrix(0, nrow = 2, ncol = 1, dimnames = list(rownames(sumI),
                                                           NULL))
      for (strand in rownames(sumI)) {
        sumS[strand,1] = sum(sumI[strand,])
      }
      
      maxSt = which.max(sumS[,1])
      
      St.values[maxSt,,] = 0
      St.errors[maxSt,,] = 0
      St.z[maxSt,,] = 0
      St.sig[maxSt,,] = 0
      
      
      
      
      if (length(St.set) > 0) {
        grStr =  "Strand.[FR]$"
        Sts = rownames(betaSummary)[grep(grStr, rownames(betaSummary))]
        if (length(Sts) > 0) {
          for (strVar in Sts) {
            St.values[strVar,,] = betaSummary[strVar,1]
            St.errors[strVar,,] = betaSummary[strVar,2]
            St.z[strVar,,]= betaSummary[strVar,3]
            St.sig[strVar,,] = betaSummary[strVar,4]
          } 
        }
        
      }
      badStr = which(is.na(St.values), arr.ind = TRUE)
      
      if (nrow(badStr) > 0) {
        minStr = which(St.values[,,1]== min(St.values[,,1], na.rm = TRUE), arr.ind = TRUE)
        St.values[badStr] = c(St.values[minStr[[1]], minStr[[2]], ]-missingValueSuppression)[1]
        St.errors[badStr] = c(St.errors[minStr[[1]], minStr[[2]], ]+missingValueSuppression)[1]
        St.z[badStr] = 0
        St.sig[badStr] = 0
      }
      
      newValues[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1] = newValues[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]+St.values[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]
      newErrors[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1] = sqrt(newErrors[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]^2+St.errors[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]^2)
      newZ[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1] = newValues[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]/newErrors[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]
      newSig[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1] = 2*pnorm(-abs(newZ[c("Strand.F", "Strand.R")[1:(1+1*(rcSymmetric == FALSE))],,1]))
      
      
      
    }
    
    
    
    updateValues(object) = newValues
    updateErrors(object) = newErrors
    updateZ(object) = newZ
    updateSig(object) = newSig
    validObject(object)
    return (object)
  }
)

#' Summary of Design Features for \linkS4class{Intercept}
#' 
#' Gets view design from design matrix, i.e. gets counts for each view feature.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @param design Data frame with design matrix, as generated by addDesignMatrix function.
#' @export
setMethod(
  f = "getDesignMatrix",
  signature = "Intercept",
  definition = function(object, design) {
    
    mat = matrix(0, ncol=dim(getValues(object))[2], nrow=2, dimnames = list(dimnames(getValues(object))[[1]],
                                                                            dimnames(getValues(object))[[2]]))
    for (st in c("F", "R")) {
      stLab = paste("Strand.", st, sep = "")
      for (bw in 1:ncol(getValues(object))) {
        bwLab = paste("View.", bw, sep = "")
        mat[stLab, bwLab] = nrow(design[((design$topMatchView == bw) & (design$topMatchStrand == st)),])
      }
    }
    
    
    return(mat)
  }
)


#' Summary of rounds design features
#' 
#' Gets round design from design matrix, i.e. gets counts for each view feature.
#' 
#' @param object object of class 'Intercept'.
#' @param design data frame with design matrix, as generated by addDesignMatrix function.
#' @export
getRoundDesign = function(object, design) {
  roundList = as.numeric(gsub("Round.", "", dimnames(object@I.values)[[3]]))
  roundCounts = hist(design$Round, breaks = seq(min(roundList)-.5, max(roundList)+.5), plot = FALSE)$counts
  rounds = hist(design$Round, breaks = seq(min(roundList)-.5, max(roundList)+.5), plot = FALSE)$mids
  roundsKeep = rounds[rounds %in% roundList]
  countsKeep = roundCounts[rounds %in% roundList]
  mat = matrix(0, ncol=length(roundList), nrow=1, dimnames = list(c("Round"),
                                                                  roundList))
  
  mat[,as.character(roundsKeep)] = countsKeep
  return(mat)
}







