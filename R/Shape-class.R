#' @include Shape-validity.R
#' @include generic-functions-definitions.R
#' @include plottingTools.R
NULL

#' Class 'Shape' for Shape Feature Coefficients 
#' 
#' Defines an S4 class to represent shape feature coefficients.
#' 
#' @slot seedLen number of base pairs in the seeding model.
#' @slot Shape.upFootprintExtend number of upstream positions to fit shape features beyond the seed.
#' @slot Shape.downFootprintExtend number of downstream positions to fit shape features beyond the seed.
#' @slot fS.upFootprintExtend maximum number of upstream positions to be fit beyond the footprint for any feature.
#' @slot fS.downFootprintExtend maximum number of downstream positions to be fit beyond the footprint for any feature. 
#' @slot fpLen footprint length for full set of features.
#' @slot Shape.set set of positions for which shape coefficients are fit where positions are given relative to the full feature set footprint.
#' @slot Shape.equivMat matrix encoding reverse complement symmetries between positions with -1 at Shape.equivMat[position1, position2] entry. All other entries are 0.  
#' @slot Shape.values beta values for positions in Shape.set.
#' @slot Shape.errors beta value errors for positions in Shape.set.
#' @slot Shape.z beta value z scores for positions in Shape.set.
#' @slot Shape.sig beta value significance for positions in Shape.set.
#' @slot Shape.oldValues previous iteration beta values for positions in Shape.set.
#' @slot Shape.oldErrors previous iteration beta values' errors for positions in Shape.set.
#' @slot Shape.oldZ previous iteration beta values' z scores for positions in Shape.set.
#' @slot Shape.oldSig previous iteration beta values' significance for positions in Shape.set.
#' @export Shape
Shape <- setClass("Shape", 
                  representation(seedLen = "numeric", 
                                 Shape.upFootprintExtend = "numeric",
                                 Shape.downFootprintExtend = "numeric",
                                 fS.upFootprintExtend = "numeric",
                                 fS.downFootprintExtend = "numeric",
                                 fpLen = "numeric",
                                 Shape.set = "numeric",
                                 shapeParamsUsed = "list",
                                 Shape.equivMat = "matrix",
                                 Shape.values = "matrix",
                                 Shape.errors = "matrix",
                                 Shape.z = "matrix",
                                 Shape.sig = "matrix",
                                 Shape.oldValues = "array",
                                 Shape.oldErrors = "array",
                                 Shape.oldZ = "array",
                                 Shape.oldSig = "array"), 
                  validity = validShape)

#' Constructor for S4 class 'Shape'
#' 
#' Initializes an S4 class to represent shape feature coefficients
#' 
#' @param shapeParamsUsed list of length 1 containing names of shape parameters to be included in model.
#' @param rcSymmetric logical: indicates the model described by the 'Shape' object will be reverse complement symmetric.
#' @export Shape
setMethod("initialize",
          "Shape", 
          function(.Object,
                   seedLen, 
                   Shape.upFootprintExtend,
                   Shape.downFootprintExtend,
                   fS.upFootprintExtend,
                   fS.downFootprintExtend,
                   fpLen,
                   Shape.set,
                   Shape.equivMat,
                   Shape.values,
                   Shape.errors,
                   Shape.z,
                   Shape.sig,
                   Shape.oldValues,
                   Shape.oldErrors,
                   Shape.oldZ,
                   Shape.oldSig,
                   shapeParamsUsed,
                   rcSymmetric,
                   ...) {
            if ((!missing(shapeParamsUsed)) & (!missing(Shape.set))) {
              if (length(shapeParamsUsed[[1]]) ==  0) {
                if ((length(Shape.set) != 1) | (Shape.set[1] != 0)) {
                  stop('If no shape parameters are specified in shapeParamsUsed, Shape.set must equal c(0).')
                }
              }
            }
            if ((missing(shapeParamsUsed)) & (missing(Shape.values)) & (missing(Shape.oldValues))) {
              stop("Must give 'shapeParamsUsed' explicitly if no shape parameters are specified by Shape.values and/or Shape.oldValues.")
            } else {
              if (!missing(shapeParamsUsed)) {
                if (is.list(shapeParamsUsed) == FALSE) {
                  stop("'shapeParamsUsed' must be a list of length 1.")
                } else if (length(shapeParamsUsed) != 1) {
                  stop("'shapeParamsUsed' must be a list of length 1.")
                }
                shapeParamsUsed[[1]] = unique(shapeParamsUsed[[1]][order(shapeParamsUsed[[1]])])
                .Object@shapeParamsUsed = shapeParamsUsed
              } else if (!missing(Shape.values)) {
                if (nrow(Shape.values) == 0) {
                  shapeParamsUsed = list(c(character(0)))
                  .Object@shapeParamsUsed = shapeParamsUsed
                } else if (is.null(rownames(Shape.values))) {
                  stop("'Shape.values' must have specified rownames.")
                } else {
                  shapeLabels = rownames(Shape.values)
                  shapeParams = c(numeric(0))
                  if ("Shape.HelTA" %in% shapeLabels) {
                    shapeLabels = shapeLabels[shapeLabels != "Shape.HelTA"]
                    if ("Shape.HelTB" %in% shapeLabels) {
                      shapeLabels = shapeLabels[shapeLabels != "Shape.HelTB"]
                      shapeParams = c(shapeParams, "HelT")
                    } else {
                      stop("'Shape.HelTA' is in rownames(Shape.values) but 'Shape.HelTB' is not in rownames(Shape.values).")
                    }
                  } else if ("Shape.HelTB" %in% shapeLabels) {
                    stop("'Shape.HelTB' is in rownames(Shape.values) but 'Shape.HelTA' is not in rownames(Shape.values).")
                  }
                  if ("Shape.MGW" %in% shapeLabels) {
                    shapeLabels = shapeLabels[shapeLabels != "Shape.MGW"]
                    shapeParams = c(shapeParams, "MGW")
                  }
                  if ("Shape.ProT" %in% shapeLabels) {
                    shapeLabels = shapeLabels[shapeLabels != "Shape.ProT"]
                    shapeParams = c(shapeParams, "ProT")
                  }
                  if ("Shape.RollA" %in% shapeLabels) {
                    shapeLabels = shapeLabels[shapeLabels != "Shape.RollA"]
                    if ("Shape.RollB" %in% shapeLabels) {
                      shapeLabels = shapeLabels[shapeLabels != "Shape.RollB"]
                      shapeParams = c(shapeParams, "Roll")
                    } else {
                      stop("'Shape.RollA' is in rownames(Shape.values) but 'Shape.RollB' is not in rownames(Shape.values).")
                    }
                  } else if ("Shape.RollB" %in% shapeLabels) {
                    stop("'Shape.RollB' is in rownames(Shape.values) but 'Shape.RollA' is not in rownames(Shape.values).")
                  }
                  if (length(shapeLabels) != 0) {
                    stop(paste("Some rownames(Shape.values) are not in c('Shape.MGW', 'Shape.ProT', 'Shape.HelTA', 'Shape.HelTB', 'Shape.RollA', 'Shape.RollB'): ",paste(shapeLabels, collapse = ", "), sep = ""))
                  }
                  shapeParamsUsed = list(shapeParams)
                  shapeParamsUsed[[1]] = shapeParamsUsed[[1]][order(shapeParamsUsed[[1]])]
                  .Object@shapeParamsUsed = shapeParamsUsed
                }
                
              } else if (!missing(Shape.oldValues)) {
                if (nrow(Shape.oldValues) == 0) {
                  shapeParamsUsed = list(c(character(0)))
                  .Object@shapeParamsUsed = shapeParamsUsed
                } else if (is.null(rownames(Shape.oldValues))) {
                  stop("'Shape.oldValues' must have specified rownames.")
                } else {
                  shapeLabels = rownames(Shape.oldValues)
                  shapeParams = c(numeric(0))
                  if ("Shape.HelTA" %in% shapeLabels) {
                    shapeLabels = shapeLabels[shapeLabels != "Shape.HelTA"]
                    if ("Shape.HelTB" %in% shapeLabels) {
                      shapeLabels = shapeLabels[shapeLabels != "Shape.HelTB"]
                      shapeParams = c(shapeParams, "HelT")
                    } else {
                      stop("'Shape.HelTA' is in rownames(Shape.values) but 'Shape.HelTB' is not in rownames(Shape.values).")
                    }
                  } else if ("Shape.HelTB" %in% shapeLabels) {
                    stop("'Shape.HelTB' is in rownames(Shape.values) but 'Shape.HelTA' is not in rownames(Shape.values).")
                  }
                  if ("Shape.MGW" %in% shapeLabels) {
                    shapeLabels = shapeLabels[shapeLabels != "Shape.MGW"]
                    shapeParams = c(shapeParams, "MGW")
                  }
                  if ("Shape.ProT" %in% shapeLabels) {
                    shapeLabels = shapeLabels[shapeLabels != "Shape.ProT"]
                    shapeParams = c(shapeParams, "ProT")
                  }
                  if ("Shape.RollA" %in% shapeLabels) {
                    shapeLabels = shapeLabels[shapeLabels != "Shape.RollA"]
                    if ("Shape.RollB" %in% shapeLabels) {
                      shapeLabels = shapeLabels[shapeLabels != "Shape.RollB"]
                      shapeParams = c(shapeParams, "Roll")
                    } else {
                      stop("'Shape.RollA' is in rownames(Shape.values) but 'Shape.RollB' is not in rownames(Shape.values).")
                    }
                  } else if ("Shape.RollB" %in% shapeLabels) {
                    stop("'Shape.RollB' is in rownames(Shape.values) but 'Shape.RollA' is not in rownames(Shape.values).")
                  }
                  if (length(shapeLabels) != 0) {
                    stop(paste("Some rownames(Shape.values) are not in c('Shape.MGW', 'Shape.ProT', 'Shape.HelTA', 'Shape.HelTB', 'Shape.RollA', 'Shape.RollB'): ",paste(shapeLabels, collapse = ", "), sep = ""))
                  }
                  shapeParamsUsed = list(shapeParams)
                  shapeParamsUsed[[1]] = shapeParamsUsed[[1]][order(shapeParamsUsed[[1]])]
                  .Object@shapeParamsUsed = shapeParamsUsed
                }
              }  else {
                stop("Must give 'shapeParamsUsed' explicitly if no shape parameters are specified by Shape.values and/or Shape.oldValues and/or Shape.set=c(0).")
              }
            }
            shapeParamsUsed[[1]] = unique(shapeParamsUsed[[1]][order(shapeParamsUsed[[1]])])
            .Object@shapeParamsUsed = shapeParamsUsed
            if (length(shapeParamsUsed[[1]]) == 0) {
              noShape = TRUE
            } else {noShape = FALSE}
            
            if (missing(rcSymmetric)) {rcSymmetric = FALSE}
            
            
            #initialize seedLen slot
            if (!missing(seedLen)) {
              .Object@seedLen <- seedLen
            } else {
              stop('seedLen must be specified')
            }
            
            # initialize fS.upFootprintExtend slot
            if (!missing(fS.upFootprintExtend)) {
              .Object@fS.upFootprintExtend<- fS.upFootprintExtend
            } else {
              stop('fS.upFootprintExtend must be specified')
            }
            
            # initialize fS.upFootprintExtend slot
            if (!missing(fS.downFootprintExtend)) {
              .Object@fS.downFootprintExtend<- fS.downFootprintExtend
            } else {
              .Object@fS.downFootprintExtend = .Object@fS.upFootprintExtend
            }
            
            # initialize fpLen
            if (!missing(fpLen)) {
              message("fpLen will be overwritten if inconsistent with fpLen = seedLen+fS.upFootprintExtend+fS.downFootprintExtend")
            }
            .Object@fpLen = .Object@seedLen+.Object@fS.upFootprintExtend+.Object@fS.downFootprintExtend
            
            # initialize Shape.upFootprintExtend, Shape.downFootprintExtend, and Shape.set
            if (noShape == FALSE) {
              if (missing(Shape.upFootprintExtend)) {
                if (missing(Shape.downFootprintExtend)) {
                  if (missing(Shape.set)) {
                    .Object@Shape.upFootprintExtend = .Object@fS.upFootprintExtend
                    .Object@Shape.downFootprintExtend = .Object@fS.downFootprintExtend
                    minSet = 1
                    maxSet = .Object@fpLen
                    .Object@Shape.set = c(minSet:maxSet)
                  } else {
                    .Object@Shape.set = Shape.set
                    minSet = min(.Object@Shape.set)
                    maxSet = max(.Object@Shape.set)
                    if (minSet == 0) {
                      .Object@Shape.upFootprintExtend = 0
                    } else if (minSet > .Object@fS.upFootprintExtend) {
                      .Object@Shape.upFootprintExtend = 0
                    } else  {
                      .Object@Shape.upFootprintExtend = .Object@fS.upFootprintExtend-(minSet-1)
                    }
                    
                    if (maxSet == 0) {
                      .Object@Shape.downFootprintExtend = 0
                    } else if (maxSet <= .Object@seedLen+.Object@fS.upFootprintExtend) {
                      .Object@Shape.downFootprintExtend = 0
                    } else {
                      .Object@Shape.downFootprintExtend = maxSet-(.Object@fS.upFootprintExtend+.Object@seedLen)
                    }
                  }
                } else {
                  .Object@Shape.downFootprintExtend = Shape.downFootprintExtend
                  if (missing(Shape.set)) {
                    .Object@Shape.upFootprintExtend = .Object@fS.upFootprintExtend
                    .Object@Shape.set = c(1:(.Object@fS.upFootprintExtend+.Object@seedLen+.Object@Shape.downFootprintExtend))
                  } else {
                    .Object@Shape.set = Shape.set
                    minSet = min(.Object@Shape.set)
                    if (minSet == 0) {
                      .Object@Shape.upFootprintExtend = 0
                    } else if (minSet <= .Object@fS.upFootprintExtend) {
                      .Object@Shape.upFootprintExtend = .Object@fS.upFootprintExtend-(minSet-1)
                    } else {
                      .Object@Shape.upFootprintExtend = 0
                    }
                  }
                }
              } else {
                .Object@Shape.upFootprintExtend = Shape.upFootprintExtend
                if (missing(Shape.downFootprintExtend)) {
                  if (missing(Shape.set)) {
                    .Object@Shape.downFootprintExtend = .Object@fS.downFootprintExtend
                    minSet = .Object@fS.upFootprintExtend-.Object@Shape.upFootprintExtend+1
                    maxSet = .Object@fS.upFootprintExtend+.Object@seedLen+.Object@Shape.downFootprintExtend
                    .Object@Shape.set = c(minSet:maxSet)
                  } else {
                    .Object@Shape.set = Shape.set
                    minSet = min(.Object@Shape.set)
                    maxSet = max(.Object@Shape.set)
                    if (maxSet > .Object@fS.upFootprintExtend+.Object@seedLen) {
                      .Object@Shape.downFootprintExtend = maxSet-(.Object@fS.upFootprintExtend+.Object@seedLen)
                    } else {
                      .Object@Shape.downFootprintExtend = 0
                    }
                  }
                } else {
                  .Object@Shape.downFootprintExtend = Shape.downFootprintExtend
                  if (missing(Shape.set)) {
                    minSet = .Object@fS.upFootprintExtend-.Object@Shape.upFootprintExtend+1
                    maxSet = .Object@fS.upFootprintExtend+.Object@seedLen+.Object@Shape.downFootprintExtend
                    .Object@Shape.set = c(minSet:maxSet)
                  } else {
                    .Object@Shape.set = Shape.set
                  }
                }
              }
            } else {
              if (!missing(Shape.set)) {
                if ((length(Shape.set) != 1) | (Shape.set[1] != 0)) {
                  stop("If no shape parameters are specified in shapeParamsUsed, Shape.set must equal c(0).")
                } else {
                  .Object@Shape.set = Shape.set
                }
              } else {
                .Object@Shape.set = c(0)
              }
              
              if (!missing(Shape.upFootprintExtend)) {
                if (Shape.upFootprintExtend != 0) {
                  stop("If no shape parameters are specified in shapeParamsUsed, Shape.set must equal c(0) and Shape.upFootprintExtend must equal 0.")
                } else {
                  .Object@Shape.upFootprintExtend = Shape.upFootprintExtend
                }
              } else {
                .Object@Shape.upFootprintExtend = 0
              }
              
              if (!missing(Shape.downFootprintExtend)) {
                if (Shape.downFootprintExtend != 0) {
                  stop("If no shape parameters are specified in shapeParamsUsed, Shape.set must equal c(0) and Shape.downFootprintExtend must equal 0.")
                } else {
                  .Object@Shape.downFootprintExtend = Shape.downFootprintExtend
                }
              } else {
                .Object@Shape.downFootprintExtend = 0
              }
            }
            
            .Object@Shape.set = unique(.Object@Shape.set)
            .Object@Shape.set = .Object@Shape.set[order(.Object@Shape.set)]
            
            
            if (rcSymmetric == TRUE) {
              if (!missing(Shape.equivMat)) {
                message('rcSymmetric == TRUE & Shape.equivMat specified. Shape.equivMat input will not be used.')
              }
              .Object@Shape.equivMat = buildSymmetricEquivalenceMatrix(.Object)
            } else {
              if (!missing(Shape.equivMat)) {
                colnames(Shape.equivMat) = c(1:ncol(Shape.equivMat))
                rownames(Shape.equivMat) = c(1:nrow(Shape.equivMat))
                .Object@Shape.equivMat = Shape.equivMat
              } else {
                .Object@Shape.equivMat = buildNullEquivalenceMatrix(.Object)
              }
            }
            
            rn = c(character(0))
            if (length(shapeParamsUsed[[1]]) > 0) {
              if ("HelT" %in% shapeParamsUsed[[1]]) {
                rn = c(rn, "Shape.HelTA", "Shape.HelTB")
              }
              if ("MGW" %in% shapeParamsUsed[[1]]) {
                rn = c(rn, "Shape.MGW")
              }
              if ("ProT" %in% shapeParamsUsed[[1]]) {
                rn = c(rn, "Shape.ProT")
              }
              if ("Roll" %in% shapeParamsUsed[[1]]) {
                rn = c(rn, "Shape.RollA", "Shape.RollB")
              }
            }
            cn = c(1:.Object@fpLen)
            
            
            if (!missing(Shape.values)) {
              if (nrow(Shape.values) != length(rn)) {
                stop("nrow('Shape.values') must be equivalent to the number of shape parameters indicated by 'shapeParamsUsed'.")
              }
              if (ncol(Shape.values) != length(cn)) {
                stop("ncol('Shape.values') must be equivalent to 'fpLen'='fS.upFootprintExtend'+'seedLen'+'fS.downFootprintExtend'.")
              }
              colnames(Shape.values) = cn
              if (!is.null(rownames(Shape.values))) {
                if (!all(rownames(Shape.values)[order(rownames(Shape.values))] == rn[order(rn)])) {
                  stop("Shape parameters indicated by 'Shape.values' do not match shape parameters indicated by 'shapeParamsUsed'.")
                }
              } else if (nrow(Shape.values) > 0) {
                stop("Inputs for 'Shape.values' must have rownames.")
              }
              if (nrow(Shape.values) > 0) {
                Shape.values = Shape.values[rn,]
              } else {
                rownames(Shape.values) = NULL
              }
              .Object@Shape.values <- Shape.values
            } else {
              .Object@Shape.values <- matrix(0, ncol = .Object@fpLen, nrow = length(rn), dimnames = list(rn, cn))
            }
            if (!missing(Shape.errors)) {
              if (nrow(Shape.errors) != length(rn)) {
                stop("nrow('Shape.errors') must be equivalent to the number of shape parameters indicated by 'shapeParamsUsed'.")
              }
              if (ncol(Shape.errors) != length(cn)) {
                stop("ncol('Shape.errors') must be equivalent to 'fpLen'='fS.upFootprintExtend'+'seedLen'+'fS.downFootprintExtend'.")
              }
              colnames(Shape.errors) = cn
              if (!is.null(rownames(Shape.errors))) {
                if (!all(rownames(Shape.errors)[order(rownames(Shape.errors))] == rn[order(rn)])) {
                  stop("Shape parameters indicated by 'Shape.errors' do not match shape parameters indicated by 'shapeParamsUsed'.")
                }
              } else if (nrow(Shape.errors) > 0) {
                stop("Inputs for 'Shape.errors' must have rownames.")
              }
              if (nrow(Shape.errors) > 0) {
                Shape.errors = Shape.errors[rn,]
              } else {
                rownames(Shape.errors) = NULL
              }
              .Object@Shape.errors <- Shape.errors
            } else {
              .Object@Shape.errors <- matrix(0, ncol = .Object@fpLen, nrow = length(rn),  dimnames = list(rn, cn))
            }
            if (!missing(Shape.z)) {
              if (nrow(Shape.z) != length(rn)) {
                stop("nrow('Shape.z') must be equivalent to the number of shape parameters indicated by 'shapeParamsUsed'.")
              }
              if (ncol(Shape.z) != length(cn)) {
                stop("ncol('Shape.z') must be equivalent to 'fpLen'='fS.upFootprintExtend'+'seedLen'+'fS.downFootprintExtend'.")
              }
              colnames(Shape.z) = cn
              if (!is.null(rownames(Shape.z))) {
                if (!all(rownames(Shape.z)[order(rownames(Shape.z))] == rn[order(rn)])) {
                  stop("Shape parameters indicated by 'Shape.z' do not match shape parameters indicated by 'shapeParamsUsed'.")
                }
              } else if (nrow(Shape.z) > 0) {
                stop("Inputs for 'Shape.z' must have rownames.")
              }
              if (nrow(Shape.z) > 0) {
                Shape.z = Shape.z[rn,]
              } else {
                rownames(Shape.z) = NULL
              }
              .Object@Shape.z <- Shape.z
            } else {
              .Object@Shape.z <- matrix(0, ncol = .Object@fpLen, nrow = length(rn),  dimnames = list(rn, cn))
            }
            if (!missing(Shape.sig)) {
              if (nrow(Shape.sig) != length(rn)) {
                stop("nrow('Shape.sig') must be equivalent to the shape parameters indicated by 'shapeParamsUsed'.")
              }
              if (ncol(Shape.sig) != length(cn)) {
                stop("ncol('Shape.sig') must be equivalent to 'fpLen'='fS.upFootprintExtend'+'seedLen'+'fS.downFootprintExtend'.")
              }
              colnames(Shape.sig) = cn
              if (!is.null(rownames(Shape.sig))) {
                if (!all(rownames(Shape.sig)[order(rownames(Shape.sig))] == rn[order(rn)])) {
                  stop("Shape parameters indicated by 'Shape.sig' do not match shape parameters indicated by 'shapeParamsUsed'.")
                }
              } else if (nrow(Shape.sig) > 0) {
                stop("Inputs for 'Shape.sig' must have rownames.")
              }
              if (nrow(Shape.sig) > 0) {
                Shape.sig = Shape.sig[rn,]
              } else {
                rownames(Shape.sig) = NULL
              }
              .Object@Shape.sig <- Shape.sig
            } else {
              .Object@Shape.sig <- matrix(0,ncol = .Object@fpLen, nrow = length(rn),  dimnames = list(rn, cn))
            }
            if (all((!missing(Shape.oldValues)), (!missing(Shape.oldZ)), (!missing(Shape.oldErrors)), (!missing(Shape.oldSig)))) {
              if (dim(Shape.oldValues)[3] == 0) {
                d3n = NULL
              } else {
                d3n = c((dim(Shape.oldValues)[3]-1):0)
              }
              
              dimnamesA = list(rn, cn, d3n)
              
              if (!missing(Shape.oldValues)) {
                if (ncol(Shape.oldValues) != length(cn)) {
                  stop("ncol('Shape.oldValues') must be equal to 'fpLen'='fS.upFootprintExtend'+'seedLen'+'fS.downFootprintExtend'.")
                }
                if (nrow(Shape.oldValues) != length(rn)) {
                  stop("nrow('Shape.oldValues') must be equivalent to the number of shape parameters indicated by 'shapeParamsUsed'.")
                }
                if (!is.null(rownames(Shape.oldValues))) {
                  if (!all(rownames(Shape.oldValues)[order(rownames(Shape.oldValues))] == rn[order(rn)])) {
                    stop("Shape parameters indicated by 'Shape.oldValues' do not match shape parameters indicated by 'shapeParamsUsed'.")
                  }
                } else if (nrow(Shape.oldValues) > 0) {
                  stop("Inputs for 'Shape.oldValues' must have rownames.")
                }
                
                if (nrow(Shape.oldValues) > 0) {
                  Shape.oldValues = Shape.oldValues[rn,,,drop=FALSE]
                } else {
                  rownames(Shape.oldValues) = NULL
                }
                .Object@Shape.oldValues <- Shape.oldValues
              }
              if (!missing(Shape.oldErrors)) {
                if (ncol(Shape.oldErrors) != length(cn)) {
                  stop("ncol('Shape.oldErrors') must be equal to 'fpLen'='fS.upFootprintExtend'+'seedLen'+'fS.downFootprintExtend'.")
                }
                if (nrow(Shape.oldErrors) != length(rn)) {
                  stop("nrow('Shape.oldErrors') must be equal to the number of shape parameters indicated by 'shapeParamsUsed'.")
                }
                if (!is.null(rownames(Shape.oldErrors))) {
                  if (!all(rownames(Shape.oldErrors)[order(rownames(Shape.oldErrors))] == rn[order(rn)])) {
                    stop("Shape parameters indicated by 'Shape.oldErrors' do not match shape parameters indicated by 'shapeParamsUsed'.")
                  }
                } else if (nrow(Shape.oldErrors) > 0) {
                  stop("Inputs for 'Shape.oldErrors' must have rownames.")
                }
                
                if (nrow(Shape.oldErrors) > 0) {
                  Shape.oldErrors = Shape.oldErrors[rn,,,drop=FALSE]
                } else {
                  rownames(Shape.oldErrors) = NULL
                }
                .Object@Shape.oldErrors <- Shape.oldErrors
              }
              if (!missing(Shape.oldZ)) {
                if (ncol(Shape.oldZ) != length(cn)) {
                  stop("ncol('Shape.oldZ') must be equivalent to 'fpLen'='fS.upFootprintExtend'+'seedLen'+'fS.downFootprintExtend'.")
                }
                if (nrow(Shape.oldZ) != length(rn)) {
                  stop("nrow('Shape.oldZ') must be equivalent to the shape parameters indicated by 'shapeParamsUsed'.")
                }
                if (!is.null(rownames(Shape.oldZ))) {
                  if (!all(rownames(Shape.oldZ)[order(rownames(Shape.oldZ))] == rn[order(rn)])) {
                    stop("Shape parameters indicated by 'Shape.oldZ' do not match shape parameters indicated by 'shapeParamsUsed'.")
                  }
                } else if (nrow(Shape.oldZ) > 0) {
                  stop("Inputs for 'Shape.oldZ' must have rownames.")
                }
                
                if (nrow(Shape.oldZ) > 0) {
                  Shape.oldZ = Shape.oldZ[rn,,,drop=FALSE]
                } else {
                  rownames(Shape.oldZ) = NULL
                }
                .Object@Shape.oldZ <- Shape.oldZ
              }
              if (!missing(Shape.oldSig)) {
                if (ncol(Shape.oldSig) != length(cn)) {
                  stop("ncol('Shape.oldSig') must be equal to 'fpLen'='fS.upFootprintExtend'+'seedLen'+'fS.downFootprintExtend'.")
                }
                if (nrow(Shape.oldSig) != length(rn)) {
                  stop("nrow('Shape.oldSig') must be equal to the number of shape parameters indicated by 'shapeParamsUsed'.")
                }
                if (!is.null(rownames(Shape.oldSig))) {
                  if (!all(rownames(Shape.oldSig)[order(rownames(Shape.oldSig))] == rn[order(rn)])) {
                    stop("Shape parameters indicated by 'Shape.oldSig' do not match shape parameters indicated by 'shapeParamsUsed'.")
                  }
                } else if (nrow(Shape.oldSig) > 0) {
                  stop("Inputs for 'Shape.oldSig' must have rownames.")
                }
                
                if (nrow(Shape.oldSig) > 0) {
                  Shape.oldSig = Shape.oldSig[rn,,,drop=FALSE]
                } else {
                  rownames(Shape.oldSig) = NULL
                }
                .Object@Shape.oldSig <- Shape.oldSig
              }
            } else if (all((missing(Shape.oldValues)), (missing(Shape.oldZ)), (missing(Shape.oldErrors)), (missing(Shape.oldSig)))) {
              d3n = NULL
              dimnamesA = list(rn, cn, d3n)
              .Object@Shape.oldValues = array(0,
                                              dim = c(length(rn),
                                                      .Object@fpLen,
                                                      0),
                                              dimnames = dimnamesA)
              .Object@Shape.oldErrors = array(0,
                                              dim = c(length(rn),
                                                      .Object@fpLen,
                                                      0),
                                              dimnames = dimnamesA)
              .Object@Shape.oldZ = array(0,
                                         dim = c(length(rn),
                                                 .Object@fpLen,
                                                 0),
                                         dimnames = dimnamesA)
              .Object@Shape.oldSig = array(0,
                                           dim = c(length(rn),
                                                   .Object@fpLen,
                                                   0),
                                           dimnames = dimnamesA)
            } else {
              stop('Must specify all or none of the following values: Shape.oldValues, Shape.oldErrors, Shape.oldZ, and Shape.oldSig.')
            }
            
            validObject(.Object)
            .Object
          })

#' Show method for \linkS4class{Shape} 
#' 
#' Defines a show method for class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "show",
  signature = "Shape",
  definition = function(object) {
    cat("An object of class '", class(object), "'\n", sep = "")
    cat("\n")
    cat('Slot "seedLen": ', object@seedLen, '\n') 
    cat("\n")
    if (nrow(object@Shape.values) > 0) {
      cat('Slot "Shape.upFootprintExtend": ', object@Shape.upFootprintExtend, '\n')
      cat("\n")
      cat('Slot "Shape.downFootprintExtend": ', object@Shape.downFootprintExtend, '\n')
      cat("\n")
    }
    cat('Slot "fS.upFootprintExtend": ', object@fS.upFootprintExtend, '\n')
    cat("\n")
    cat('Slot "fS.downFootprintExtend": ', object@fS.downFootprintExtend, '\n')
    cat("\n")
    cat('Slot "fpLen": ', object@fpLen, '\n')
    cat("\n")
    if (length(object@shapeParamsUsed) > 0) {
      cat('Slot "shapeParamsUsed[[1]]": ', object@shapeParamsUsed[[1]], '\n')
      cat("\n")
    }
    if (nrow(object@Shape.values) > 0) {
      cat('Slot "Shape.set": ', object@Shape.set, '\n')
      cat("\n")
      cat('Slot "Shape.equivMat":\n')
      cat("\n")
    }
    
    if (all(object@Shape.equivMat == 0)) {
      cat(object@fpLen," x ", object@fpLen, " null equivalence matrix")
      cat("\n")
    } else {
      print(object@Shape.equivMat)
      cat("\n")
    }
    cat("\n")
    if (nrow(object@Shape.values) > 0) {
      cat('Slot "Shape.values":\n')
      print(object@Shape.values)
      cat("\n")
      cat("\n")
      cat('Slot "Shape.errors":\n')
      print(object@Shape.errors)
      cat("\n")
      cat("\n")
      cat('Slot "Shape.z":\n')
      print(object@Shape.z)
      cat("\n")
      cat("\n")
      cat('Slot "Shape.sig":\n')
      print(object@Shape.sig)
      cat("\n")
      cat("\n")
      cat('Slot "Shape.oldValues":\n')
      if (dim(object@Shape.oldValues)[3] == 0) {
        cat('<',
            dim(object@Shape.oldValues)[1],
            ' x ',
            dim(object@Shape.oldValues)[2],
            ' x ',
            dim(object@Shape.oldValues)[3],
            ' array of double>\n', sep = "")
      } else {
        print(object@Shape.oldValues)
        cat("\n")
      }
      cat("\n")
      cat('Slot "Shape.oldErrors":\n')
      if (dim(object@Shape.oldErrors)[3] == 0) {
        cat('<',
            dim(object@Shape.oldErrors)[1],
            ' x ',
            dim(object@Shape.oldErrors)[2],
            ' x ',
            dim(object@Shape.oldErrors)[3],
            ' array of double>\n', sep = "")
      } else {
        print(object@Shape.oldErrors)
        cat("\n")
      }
      cat("\n")
      cat('Slot "Shape.oldZ":\n')
      if (dim(object@Shape.oldZ)[3] == 0) {
        cat('<',
            dim(object@Shape.oldZ)[1],
            ' x ',
            dim(object@Shape.oldZ)[2],
            ' x ',
            dim(object@Shape.oldZ)[3],
            ' array of double>\n', sep = "")
      } else {
        print(object@Shape.oldZ)
        cat("\n")
      }
      cat("\n")
      cat('Slot "Shape.oldSig":\n')
      if (dim(object@Shape.oldSig)[3] == 0) {
        cat('<',
            dim(object@Shape.oldSig)[1],
            ' x ',
            dim(object@Shape.oldSig)[2],
            ' x ',
            dim(object@Shape.oldSig)[3],
            ' array of double>\n', sep = "")
      } else {
        print(object@Shape.oldSig)
        cat("\n")
      }
      cat("\n")
      
    } else {
      cat("\nNo shape parameters included in 'model' fit. Shape.values, Shape.errors, Shape.z, Shape.sig, ",
          "Shape.oldValues, Shape.oldErrors, Shape.oldZ, and Shape.oldSig have zero rows.")
    }
    
    cat("\n")
    invisible(NULL)
  }
)


#' Summary method for \linkS4class{Shape}
#' 
#' Defines a summary method for class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "summary",
  signature = "Shape",
  definition = function(object) {
    cat("An object of class '", class(object), "'\n", sep = "")
    cat("Fits shape coefficients at ",
        length(object@Shape.set[object@Shape.set != 0]),
        " positions (",
        paste(object@Shape.set[object@Shape.set != 0], collapse = ", "),
        ") for shape parameter(s) (",
        paste(object@shapeParamsUsed[[1]], collapse = ", "),
        ") in a feature model of length ",
        (object@fpLen),
        ".\n",
        sep = "")
    
    if (nrow(object@Shape.values) > 0) {
      if (!all(object@Shape.equivMat == 0)) {
        if (all(object@Shape.equivMat[,1:(object@fpLen %/%2)] == 0)) {
          halfMatCheck = object@Shape.equivMat[1:(object@fpLen %/%2+object@fpLen %%2),
                                               (object@fpLen %/%2+1):object@fpLen]
          halfMatCheck = halfMatCheck[,ncol(halfMatCheck):1]
          if (all(diag(halfMatCheck) == -1)) {
            cat("Shape features are reverse complement symmetric.\n")
          }
        }
      }
      cat("Shape beta values:\n")
      print(getValues(object))
      cat("\n")
      cat("Shape beta errors:\n")
      print(getErrors(object))
      cat("\n")
      cat("Shape beta z scores:\n")
      print(getZ(object))
      cat("\n")
      cat("Shape beta p-values:\n")
      print(getSig(object))
      cat("\n")
    } 
    invisible(NULL)
  }
)

#' Print Method for \linkS4class{Shape}
#' 
#' Defines a print method for class \linkS4class{Shape}.
#' 
#' @param x Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "print",
  signature = "Shape",
  definition = function(x) {
    show(x)
  }
)

#' Get Feature Design for \linkS4class{Shape} 
#' 
#' Defines a print design method for class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getFeatureDesign",
  signature = "Shape",
  definition = function(object) {
    cat("An object of class '", class(object), "'\n", sep = "")
    cat('Slot "seedLen": ', object@seedLen, '\n')
    cat('Slot "Shape.upFootprintExtend": ', object@Shape.upFootprintExtend, '\n')
    cat('Slot "Shape.downFootprintExtend": ', object@Shape.downFootprintExtend, '\n')
    cat('Slot "fS.upFootprintExtend": ', object@fS.upFootprintExtend, '\n')
    cat('Slot "fS.downFootprintExtend": ', object@fS.downFootprintExtend, '\n')
    cat('Slot "Shape.set": ', object@Shape.set, '\n')
    cat('"ShapeParamsUsed": ', gsub("Shape.", "",rownames(object@Shape.values)), '\n')
    cat("\n")
    invisible(NULL)
  }
)

#' Length Method for \linkS4class{Shape}
#' 
#' Defines a length method for class \linkS4class{Shape}.
#' 
#' @param x Object of class \linkS4class{Shape}.
#' @return Footprint length for shape parameters.
#' @export
setMethod(
  f = "length",
  signature = "Shape",
  definition = function(x) {
    if (nrow(x@Shape.values) > 0) {
      return (x@fpLen)
    } else {
      return (0)
    }
    
  }
)

#' Seed Length for \linkS4class{Shape}
#' 
#' Gets seedLen slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getSeedLen",
  signature = "Shape",
  definition = function(object) {
    return(object@seedLen)
  }
)

#' Get Length (bp) of Upstream Footprint Extension for \linkS4class{Shape}
#' 
#' Gets Shape.upFootprintExtend slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getUpFootprintExtend",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.upFootprintExtend)
  }
)

#' Get Length (bp) of Downstream Footprint Extension for \linkS4class{Shape}
#' 
#' Gets Shape.downFootprintExtend slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getDownFootprintExtend",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.downFootprintExtend)
  }
)

#' Get Length (bp) of Upstream Footprint Extension for Full Feature Set
#' 
#' Gets fS.upFootprintExtend slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getFsUpFootprintExtend",
  signature = "Shape",
  definition = function(object) {
    return(object@fS.upFootprintExtend)
  }
)

#' Get Length (bp) of Downstream Footprint Extension for Full Feature Set
#' 
#' Gets fS.downFootprintExtend slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getFsDownFootprintExtend",
  signature = "Shape",
  definition = function(object) {
    return(object@fS.downFootprintExtend)
  }
)

#' Get Footprint Length for \linkS4class{Shape}
#' 
#' Gets fpLen slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getFpLen",
  signature = "Shape",
  definition = function(object) {
    return(object@fpLen)
  }
)

#' Get positions for Shape Parameter Fitting from \linkS4class{Shape}
#' 
#' Gets Shape.set slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getPositions",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.set)
  }
)


#' Get Shape Parameters Used from \linkS4class{Shape}
#' 
#' Gets 'shapeParamsUsed' slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getShapeParamsUsed",
  signature = "Shape",
  definition = function(object) {
    return(object@shapeParamsUsed)
  }
)



#' Get Equivalence Matrix for \linkS4class{Shape}
#' 
#' Gets Shape.equivMat slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getEquivMat",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.equivMat)
  }
)

#' Get Values for \linkS4class{Shape}
#' 
#' Gets Shape.values slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getValues",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.values)
  }
)



#' Get Errors for \linkS4class{Shape}
#' 
#' Gets Shape.errors slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getErrors",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.errors)
  }
)

#' Get Z-scores for \linkS4class{Shape}
#' 
#' Gets Shape.z slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getZ",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.z)
  }
)

#' Get P-values for \linkS4class{Shape}
#' 
#' Gets Shape.sig slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getSig",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.sig)
  }
)

#' Get Values from previous iterations for \linkS4class{Shape}
#' 
#' Gets Shape.oldValues slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getOldValues",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.oldValues)
  }
)

#' Get Errors from previous iterations for \linkS4class{Shape}
#' 
#' Gets Shape.oldErrors slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getOldErrors",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.oldErrors)
  }
)

#' Get Z-scores from previous iterations for \linkS4class{Shape}
#' 
#' Gets Shape.oldZ slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getOldZ",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.oldZ)
  }
)

#' Get P-values from previous iterations for \linkS4class{Shape}
#' 
#' Gets Shape.oldSig slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
setMethod(
  f = "getOldSig",
  signature = "Shape",
  definition = function(object) {
    return(object@Shape.oldSig)
  }
)

#' Set Seed Length for \linkS4class{Shape}
#' 
#' Sets seedLen slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Seed length for \linkS4class{model} object to which \linkS4class{Shape} object belongs.
#' @export
setReplaceMethod(
  f = "setSeedLen",
  signature = "Shape",
  definition = function(object, value) {
    stop('Changes to the seed length cannot be made at the level of "Shape" class objects. Note that any changes in the seed length must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)

#' Set Length (bp) of Upstream Footprint Extension for \linkS4class{Shape}
#' 
#' Sets Shape.upFootprintExtend slot for object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Number of bp upstream of seed model to fit shape coefficients for.
#' @export
setReplaceMethod(
  f = "setUpFootprintExtend",
  signature = "Shape",
  definition = function(object, value) {
    if (value > object@fS.upFootprintExtend) {
      stop('Shape.upFootprintExtend must be <= fS.upFootprintExtend')
    }
    if (value < 0) {
      stop('Shape.upFootprintExtend must be >= 0.')
    }
    if (nrow(object@Shape.values) == 0) {
      stop('No shape parameters specified, so shape.set cannot be expanded.')
    }
    if (value == object@Shape.upFootprintExtend) {
      validObject(object)
      return(object)
    } else if (value > object@Shape.upFootprintExtend) {
      diff = value-object@Shape.upFootprintExtend
      min = object@fS.upFootprintExtend-value+1
      seqAdd = seq(min, min+diff-1)
      object@Shape.set = unique(c(seqAdd, object@Shape.set))[order(unique(c(seqAdd, object@Shape.set)))]
      if (value > 0) {
        fullSeqNew = min-1+seq(1, value)
        notIncluded = fullSeqNew[!fullSeqNew %in% object@Shape.set]
        if (length(notIncluded) > 0) {
          message(paste('Based on previous value for Shape.set, the following bases implied by the new value for Shape.upFootprintExtend are excluded from N.set: ',
                        paste(notIncluded, collapse = ", "),
                        sep = ""))
        }
      }
      
      
    } else {
      if (!all(object@Shape.values == 0)) {
        stop('Cannot deprecate value of Shape.upFootprintExtend after beta values have been added.')
      } else {
        newSetMin = object@fS.upFootprintExtend-value+1
        object@Shape.set = object@Shape.set[object@Shape.set >= newSetMin]
        if (!(newSetMin %in% object@Shape.set)) {
          message(paste('New Shape.upFootprintExtend adds the following nucleotide values not previously included in Shape.set: ',
                        newSetMin,
                        sep = ""))
        }
        object@Shape.set = unique(c(newSetMin, object@Shape.set))[order(unique(c(newSetMin, object@Shape.set)))]
        if (value > 0) {
          fullSeqNew = newSetMin-1+seq(1, value)
          notIncluded = fullSeqNew[!fullSeqNew %in% object@Shape.set]
          if (length(notIncluded) > 0) {
            message(paste('Based on previous value for Shape.set, the following bases implied by the new value for Shape.upFootprintExtend are excluded from Shape.set: ',
                          paste(notIncluded, collapse = ", "),
                          sep = ""))
          }
        }
        
      }
    }
    object@Shape.upFootprintExtend = value
    message('Shape.downFootprintExtend is not changed automatically to match Shape.upFootprintExtend')
    validObject(object)
    return(object)
  }
)


#' Set Length (bp) of Upstream Footprint Extension for Full Feature Model
#' 
#' Sets fS.upFootprintExtend slot for object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value number of bp upstream of seed model over which the full 'featureSet' object is defined.
#' @export
setReplaceMethod(
  f = "setFsUpFootprintExtend",
  signature = "Shape",
  definition = function(object, value) {
    stop('Changes to fS.upFootprintExtend cannot be made at the level of "Shape" class objects. Note that any changes in fS.upFootprintExtend must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)

#' Set Length (bp) of Downstream Footprint Extension for \linkS4class{Shape}
#' 
#' Sets Shape.downFootprintExtend slot for object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value number of bp downstream of seed model over which to fit mono-nucleotides for.
#' @export
setReplaceMethod(
  f = "setDownFootprintExtend",
  signature = "Shape",
  definition = function(object, value) {
    if (value > object@fS.downFootprintExtend) {
      stop('Shape.downFootprintExtend must be <= fS.downFootprintExtend')
    }
    if (value < 0) {
      stop('Shape.downFootprintExtend must be >= 0.')
    }
    if (nrow(object@Shape.values) == 0) {
      stop('No shape parameters specified, so shape.set cannot be expanded.')
    }
    if (value == object@Shape.downFootprintExtend) {
      validObject(object)
      return(object)
    } else if (value > object@Shape.downFootprintExtend) {
      diff = value-object@Shape.downFootprintExtend
      min = object@fS.upFootprintExtend+object@seedLen+object@Shape.downFootprintExtend+1
      seqAdd = seq(min, min+diff-1)
      object@Shape.set = unique(c(seqAdd, object@Shape.set))[order(unique(c(seqAdd, object@Shape.set)))]
      
      fullSeqNew = seq(object@fS.upFootprintExtend+object@seedLen+1, min+diff-1)
      notIncluded = fullSeqNew[!fullSeqNew %in% object@Shape.set]
      if (length(notIncluded) > 0) {
        message(paste('Based on previous value for Shape.set, the following bases implied by the new value for Shape.downFootprintExtend are excluded from Shape.set: ',
                      paste(notIncluded, collapse = ", "),
                      sep = ""))
      }
      
    } else {
      if (!all(object@Shape.values == 0)) {
        stop('Cannot deprecate value of Shape.downFootprintExtend after beta values have been added.')
      } else {
        newSetMax = object@fS.upFootprintExtend+object@seedLen+value
        object@Shape.set = object@Shape.set[object@Shape.set <= newSetMax]
        if (!(newSetMax %in% object@Shape.set)) {
          message(paste('New Shape.downFootprintExtend adds the following nucleotide values not previously included in Shape.set: ',
                        newSetMax,
                        sep = ""))
        }
        object@Shape.set = unique(c(newSetMax, object@Shape.set))[order(unique(c(newSetMax, object@Shape.set)))]
        if (value > 0) {
          fullSeqNew = object@fS.upFootprintExtend+object@seedLen+seq(1, value)
          notIncluded = fullSeqNew[!fullSeqNew %in% object@Shape.set]
          if (length(notIncluded) > 0) {
            message(paste('Based on previous value for Shape.set, the following bases implied by the new value for Shape.downFootprintExtend are excluded from Shape.set: ',
                          paste(notIncluded, collapse = ", "),
                          sep = ""))
          } 
        }
        
      }
    } 
    object@Shape.downFootprintExtend = value
    message('Shape.upFootprintExtend is not changed automatically to match Shape.downFootprintExtend')
    validObject(object)
    return(object)
  }
)


#' Set Length (bp) of Downstream Footprint Extension for Full Feature Model
#' 
#' Sets fS.downFootprintExtend slot for object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value number of bp upstream of seed model over which the full 'featureSet' object is defined.
#' @export
setReplaceMethod(
  f = "setFsDownFootprintExtend",
  signature = "Shape",
  definition = function(object, value) {
    stop('Changes to fS.downFootprintExtend cannot be made at the level of "Shape" class objects. Note that any changes in fS.downFootprintExtend must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)

#' Set Footprint Length for \linkS4class{Shape}
#' 
#' Sets fpLen slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Footprint length for full \linkS4class{featureSet} object.
#' @export
setReplaceMethod(
  f = "setFpLen",
  signature = "Shape",
  definition = function(object, value) {
    stop('Changes to fpLen cannot be made at the level of "Shape" class objects. Note that any changes in fpLen must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)

#' Set positions for Shape Parameter Fitting from \linkS4class{Shape}
#' 
#' Sets Shape.set slot for object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Vector of positions within \linkS4class{featureSet} for which shape parameter coefficients should be fit. 
#' @export
setReplaceMethod(
  f = "setPositions",
  signature = "Shape",
  definition = function(object, value) {
    if (max(value) > object@fpLen) {
      stop('New Shape.set value is not consistent with fpLen.')
    }
    if (nrow(object@Shape.values) == 0) {
      stop('No shape parameters specified, so shape.set cannot be expanded.')
    }
    if (!all(object@Shape.set %in% value)) {
      if (!all(object@Shape.values == 0)) {
        stop('Cannot deprecate value of Shape.set after beta values have been added.')
      }
    }
    minSet = min(value)
    if (minSet != min(object@Shape.set)) {
      object@Shape.upFootprintExtend = max(object@fS.upFootprintExtend-minSet+1, 0)
      message('New Shape.set implies new value for Shape.upFootprintExtend: ', object@Shape.upFootprintExtend)
    }
    maxSet = max(value)
    if (maxSet != max(object@Shape.set)) {
      object@Shape.downFootprintExtend = max(maxSet-object@fS.upFootprintExtend-object@seedLen, 0)
      message('New Shape.set implies new value for Shape.downFootprintExtend: ', object@Shape.downFootprintExtend)
    }
    
    
    object@Shape.set = value
    validObject(object)
    return(object)
  }
)

#' Sets Equivalence Matrix for \linkS4class{Shape}
#' 
#' Sets an equivalence matrix for object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Equivalence matrix for object of class \linkS4class{Shape}.
#' @export
setReplaceMethod(
  f = "setEquivMat",
  signature = "Shape",
  definition = function(object, value) {
    colnames(value) = c(1:ncol(value))
    rownames(value) = c(1:nrow(value))
    object@Shape.equivMat = value
    object@Shape.equivMat = formatEquivalenceMatrix(object)
    validObject(object)
    return(object)
  }
)

#' Set Shape Parameters Used for \linkS4class{Shape}
#' 
#' Sets shapeParamsUsed slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value ShapeParamsUsed for \linkS4class{model} object to which \linkS4class{Shape} object belongs.
#' @export
setReplaceMethod(
  f = "setShapeParamsUsed",
  signature = "Shape",
  definition = function(object, value) {
    stop('Changes to shapeParamsUsed slot cannot be made at the level of "Shape" class objects. "shapeParamsUsed" slot should be changed in the full "model" object.  Note that any changes in the shapeParamsUsed slot must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)


#' Set Values for \linkS4class{Shape}
#' 
#' Sets Shape.values slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Matrix of length(shapeParamsUsed[[1]]) x fpLen containing beta values for shape coefficients in -ddG/angstrom units.
#' @export
setReplaceMethod(
  f = "setValues",
  signature = "Shape",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@Shape.values)) {
      stop('Current and replacement versions of Shape.values must have the same dimensions.')
    }
    if (nrow(value) != nrow(object@Shape.values)) {
      stop('Current and replacement versions of Shape.values must have the same dimensions.')
    }
    if (!is.null(rownames(value))) {
      if (!all(unique(rownames(value))[unique(order(rownames(value)))] == unique(rownames(object@Shape.values))[order(unique(rownames(object@Shape.values)))])) {
        stop('Replacement shape parameters do not match shape parameters in current Shape.values.')
      } 
    } else if (nrow(value) > 0) {
      stop('Replacement shape parameters must have rownames that match the old rownames.')
    }
    if (nrow(value) > 0) {
      value = value[rownames(object@Shape.values),]
    } else {
      rownames(value) = NULL
    }
    colnames(value) = 1:object@fpLen
    object@Shape.values <- value
    validObject(object)
    return(object)
  }
)

#' Set Errors for \linkS4class{Shape}
#' 
#' Sets Shape.errors slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Matrix of length(shapeParamsUsed[[1]]) x fpLen containing beta errors for shape coefficients in -ddG/angstrom units.
#' @export
setReplaceMethod(
  f = "setErrors",
  signature = "Shape",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@Shape.errors)) {
      stop('Current and replacement versions of Shape.errors must have the same dimensions.')
    }
    if (nrow(value) != nrow(object@Shape.errors)) {
      stop('Current and replacement versions of Shape.errors must have the same dimensions.')
    }
    if (!is.null(rownames(value))) {
      if (!all(unique(rownames(value))[unique(order(rownames(value)))] == unique(rownames(object@Shape.errors))[order(unique(rownames(object@Shape.errors)))])) {
        stop('Replacement shape parameters do not match shape parameters in current Shape.errors.')
      } 
    } else if (nrow(value) > 0) {
      stop('Replacement shape parameters must have rownames that match the old rownames.')
    }
    if (nrow(value) > 0) {
      value = value[rownames(object@Shape.errors),]
    } else {
      rownames(value) = NULL
    }
    colnames(value) = 1:object@fpLen
    object@Shape.errors <- value
    validObject(object)
    return(object)
  }
)


#' Set Z-scores for \linkS4class{Shape}
#' 
#' Sets Shape.z slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Matrix of length(shapeParamsUsed[[1]]) x fpLen containing beta z for shape coefficients in -ddG/angstrom units.
#' @export
setReplaceMethod(
  f = "setZ",
  signature = "Shape",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@Shape.z)) {
      stop('Current and replacement versions of Shape.z must have the same dimensions.')
    }
    if (nrow(value) != nrow(object@Shape.z)) {
      stop('Current and replacement versions of Shape.z must have the same dimensions.')
    }
    if (!is.null(rownames(value))) {
      if (!all(unique(rownames(value))[unique(order(rownames(value)))] == unique(rownames(object@Shape.z))[order(unique(rownames(object@Shape.z)))])) {
        stop('Replacement shape parameters do not match shape parameters in current Shape.z.')
      } 
    } else if (nrow(value) > 0){
      stop('Replacement shape parameters must have rownames that match the old rownames.')
    }
    if (nrow(value) > 0) {
      value = value[rownames(object@Shape.z),]
    } else {
      rownames(value) = NULL
    }
    colnames(value) = 1:object@fpLen
    object@Shape.z <- value
    validObject(object)
    return(object)
  }
)

#' Set P-values for \linkS4class{Shape}
#' 
#' Sets Shape.sig slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Matrix of length(shapeParamsUsed[[1]]) x fpLen containing beta sig for shape coefficients in -ddG/angstrom units.
#' @export
setReplaceMethod(
  f = "setSig",
  signature = "Shape",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@Shape.sig)) {
      stop('Current and replacement versions of Shape.sig must have the same dimensions.')
    }
    if (nrow(value) != nrow(object@Shape.sig)) {
      stop('Current and replacement versions of Shape.sig must have the same dimensions.')
    }
    if (!is.null(rownames(value))) {
      if (!all(unique(rownames(value))[unique(order(rownames(value)))] == unique(rownames(object@Shape.sig))[order(unique(rownames(object@Shape.sig)))])) {
        stop('Replacement shape parameters do not match shape parameters in current Shape.sig.')
      } 
    } else if (nrow(value) > 0) {
      stop('Replacement shape parameters must have rownames that match the old rownames.')
    }
    if (nrow(value) > 0) {
      value = value[rownames(object@Shape.sig),]
    } else {
      rownames(value) = NULL
    }
    colnames(value) = 1:object@fpLen
    object@Shape.sig <- value
    validObject(object)
    return(object)
  }
)

#' Set Values for Previous Iterations for \linkS4class{Shape}
#' 
#' Sets Shape.oldValues slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Array of length(shapeParamsUsed[[1]]) x fpLen x (iterations-1) containing beta values for shape coefficients in -ddG/angstrom units.
#' @export
setReplaceMethod(
  f = "setOldValues",
  signature = "Shape",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@Shape.oldValues)) {
      stop('Current and replacement versions of Shape.oldValues must have the same dimensions.')
    }
    if (nrow(value) != nrow(object@Shape.oldValues)) {
      stop('Current and replacement versions of Shape.oldValues must have the same dimensions.')
    }
    if (!is.null(rownames(value))) {
      if (!all(unique(rownames(value))[unique(order(rownames(value)))] == unique(rownames(object@Shape.oldValues))[order(unique(rownames(object@Shape.oldValues)))])) {
        stop('Replacement shape parameters do not match shape parameters in current Shape.oldValues.')
      }
    } else if (nrow(value) > 0) {
      stop('Replacement shape parameters must have rownames that match the old rownames.')
    }
    if (dim(value)[1] > 0) {
      if (dim(value)[3] > 0) {
        value = value[dimnames(object@Shape.oldValues)[[1]],,,drop=FALSE]
        dimnames(value) = list(c(dimnames(object@Shape.oldValues)[[1]]),
                               c(1:object@fpLen),
                               c((dim(value)[3]-1):0))
      } else {
        value = value[dimnames(object@Shape.oldValues)[[1]],,,drop=FALSE]
        dimnames(value) = list(c(dimnames(object@Shape.oldValues)[[1]]),
                               c(1:object@fpLen),
                               NULL)
      }
    } else {
      if (dim(value)[3] > 0) {
        dimnames(value) = list(NULL,
                               c(1:object@fpLen),
                               c((dim(value)[3]-1):0))
      } else {
        dimnames(value) = list(NULL,
                               c(1:object@fpLen),
                               NULL)
      }
    }
    
    object@Shape.oldValues <- value
    validObject(object)
    return(object)
  }
)

#' Set Errors for Previous Iterations for \linkS4class{Shape}
#' 
#' Sets Shape.oldErrors slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Array of length(shapeParamsUsed[[1]]) x fpLen x (iterations-1) containing beta errors for shape coefficients in -ddG/angstrom units.
#' @export
setReplaceMethod(
  f = "setOldErrors",
  signature = "Shape",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@Shape.oldErrors)) {
      stop('Current and replacement versions of Shape.oldErrors must have the same dimensions.')
    }
    if (nrow(value) != nrow(object@Shape.oldErrors)) {
      stop('Current and replacement versions of Shape.oldErrors must have the same dimensions.')
    }
    if (!is.null(rownames(value))) {
      if (!all(unique(rownames(value))[unique(order(rownames(value)))] == unique(rownames(object@Shape.oldErrors))[order(unique(rownames(object@Shape.oldErrors)))])) {
        stop('Replacement shape parameters do not match shape parameters in current Shape.oldErrors.')
      }
    } else if (nrow(value) > 0) {
      stop('Replacement shape parameters must have rownames that match the old rownames.')
    }
    
    if (dim(value)[1] > 0) {
      if (dim(value)[3] > 0) {
        value = value[dimnames(object@Shape.oldErrors)[[1]],,,drop=FALSE]
        dimnames(value) = list(c(dimnames(object@Shape.oldErrors)[[1]]),
                               c(1:object@fpLen),
                               c((dim(value)[3]-1):0))
      } else {
        value = value[dimnames(object@Shape.oldErrors)[[1]],,,drop=FALSE]
        dimnames(value) = list(c(dimnames(object@Shape.oldErrors)[[1]]),
                               c(1:object@fpLen),
                               NULL)
      }
    } else {
      if (dim(value)[3] > 0) {
        dimnames(value) = list(NULL,
                               c(1:object@fpLen),
                               c((dim(value)[3]-1):0))
      } else {
        dimnames(value) = list(NULL,
                               c(1:object@fpLen),
                               NULL)
      }
    }
    object@Shape.oldErrors <- value
    validObject(object)
    return(object)
  }
)

#' Set Z-scores for Previous Iterations for \linkS4class{Shape}
#' 
#' Sets Shape.oldZ slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Array of length(shapeParamsUsed[[1]]) x fpLen x (iterations-1) containing beta z-scores for shape coefficients in -ddG/angstrom units.
#' @export
setReplaceMethod(
  f = "setOldZ",
  signature = "Shape",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@Shape.oldZ)) {
      stop('Current and replacement versions of Shape.oldZ must have the same dimensions.')
    }
    if (nrow(value) != nrow(object@Shape.oldZ)) {
      stop('Current and replacement versions of Shape.oldZ must have the same dimensions.')
    }
    if (!is.null(rownames(value))) {
      if (!all(unique(rownames(value))[unique(order(rownames(value)))] == unique(rownames(object@Shape.oldZ))[order(unique(rownames(object@Shape.oldZ)))])) {
        stop('Replacement shape parameters do not match shape parameters in current Shape.oldZ.')
      }
    } else if (nrow(value) > 0) {
      stop('Replacement shape parameters must have rownames that match the old rownames.')
    }
    
    if (dim(value)[1] > 0) {
      if (dim(value)[3] > 0) {
        value = value[dimnames(object@Shape.oldZ)[[1]],,,drop=FALSE]
        dimnames(value) = list(c(dimnames(object@Shape.oldZ)[[1]]),
                               c(1:object@fpLen),
                               c((dim(value)[3]-1):0))
      } else {
        value = value[dimnames(object@Shape.oldZ)[[1]],,,drop=FALSE]
        dimnames(value) = list(c(dimnames(object@Shape.oldZ)[[1]]),
                               c(1:object@fpLen),
                               NULL)
      }
    } else {
      if (dim(value)[3] > 0) {
        dimnames(value) = list(NULL,
                               c(1:object@fpLen),
                               c((dim(value)[3]-1):0))
      } else {
        dimnames(value) = list(NULL,
                               c(1:object@fpLen),
                               NULL)
      }
    }
    object@Shape.oldZ <- value
    validObject(object)
    return(object)
  }
)

#' Set P-values from Previous Iterations for \linkS4class{Shape}
#' 
#' Sets Shape.oldSig slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Array of length(shapeParamsUsed[[1]]) x fpLen x (iterations-1) containing beta p-values for shape coefficients in -ddG/angstrom units.
#' @export
setReplaceMethod(
  f = "setOldSig",
  signature = "Shape",
  definition = function(object, value) {
    if (ncol(value) != ncol(object@Shape.oldSig)) {
      stop('Current and replacement versions of Shape.oldSig must have the same dimensions.')
    }
    if (nrow(value) != nrow(object@Shape.oldSig)) {
      stop('Current and replacement versions of Shape.oldSig must have the same dimensions.')
    }
    if (!is.null(rownames(value))) {
      if (!all(unique(rownames(value))[unique(order(rownames(value)))] == unique(rownames(object@Shape.oldSig))[order(unique(rownames(object@Shape.oldSig)))])) {
        stop('Replacement shape parameters do not match shape parameters in current Shape.oldSig.')
      }
    }else if (nrow(value) > 0) {
      stop('Replacement shape parameters must have rownames that match the old rownames.')
    }
    
    if (dim(value)[1] > 0) {
      if (dim(value)[3] > 0) {
        value = value[dimnames(object@Shape.oldSig)[[1]],,,drop=FALSE]
        dimnames(value) = list(c(dimnames(object@Shape.oldSig)[[1]]),
                               c(1:object@fpLen),
                               c((dim(value)[3]-1):0))
      } else {
        value = value[dimnames(object@Shape.oldSig)[[1]],,,drop=FALSE]
        dimnames(value) = list(c(dimnames(object@Shape.oldSig)[[1]]),
                               c(1:object@fpLen),
                               NULL)
      }
    } else {
      if (dim(value)[3] > 0) {
        dimnames(value) = list(NULL,
                               c(1:object@fpLen),
                               c((dim(value)[3]-1):0))
      } else {
        dimnames(value) = list(NULL,
                               c(1:object@fpLen),
                               NULL)
      }
    }
    
    
    object@Shape.oldSig <- value
    validObject(object)
    return(object)
  }
)


#' Update Values for \linkS4class{Shape}
#' 
#' Updates Shape.values slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Matrix of 4xfpLen containing beta values for shape coefficients in -ddG/angstrom units.
setReplaceMethod(
  f = "updateValues",
  signature = "Shape",
  definition = function(object, value){
    if (nrow(object@Shape.values) != nrow(value)) {
      stop("Dimensions for updated values are not consistent with 'Shape' object.")
    }
    oldValues = getOldValues(object)
    # N object stores values in matrix and oldValues in an array
    updateDim = dim(oldValues)+c(0,0,1)
    oV.updated = array(0, dim = updateDim)
    if (dim(oldValues)[3] > 0) {
      oV.updated[,,2:(dim(oV.updated)[3])] = oldValues 
    }
    oV.updated[,,1] = getValues(object)
    #Shape specific dimensionality for labels
    if (nrow(value) > 0) {
      rn = rownames(getValues(object))
    }  else {
      rn = NULL
    }
    cn = colnames(getValues(object))
    d3n = c((updateDim[3]-1):0)
    rownames(value) = rn 
    colnames(value) = cn
    dimnames = list(rn, cn, d3n)
    dimnames(oV.updated) = dimnames
    object@Shape.oldValues <- oV.updated
    object@Shape.values <- value
    return (object)
  }
)

#' Update Errors for \linkS4class{Shape}
#' 
#' Updates Shape.errors slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Matrix of 4xfpLen containing beta errors for shape coefficients in -ddG/angstrom units.
setReplaceMethod(
  f = "updateErrors",
  signature = "Shape",
  definition = function(object, value){
    if (nrow(object@Shape.errors) != nrow(value)) {
      stop("Dimensions for updated errors are not consistent with 'Shape' object.")
    }
    oldValues = getOldErrors(object)
    # N object stores values in matrix and oldValues in an array
    updateDim = dim(oldValues)+c(0,0,1)
    oV.updated = array(0, dim = updateDim)
    if (dim(oldValues)[3] > 0) {
      oV.updated[,,2:(dim(oV.updated)[3])] = oldValues 
    }
    oV.updated[,,1] = getErrors(object)
    #Shape specific dimensionality for labels
    if (nrow(value) > 0) {
      rn = rownames(getErrors(object))
    }  else {
      rn = NULL
    }
    cn = colnames(getErrors(object))
    d3n = c((updateDim[3]-1):0)
    rownames(value) = rn 
    colnames(value) = cn
    dimnames = list(rn, cn, d3n)
    dimnames(oV.updated) = dimnames
    object@Shape.oldErrors <- oV.updated
    object@Shape.errors <- value
    return (object)
  }
)

#' Update Z-scores for \linkS4class{Shape}
#' 
#' Updates Shape.z slot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Matrix of 4xfpLen containing beta z-scores for shape coefficients in -ddG/angstrom units.
setReplaceMethod(
  f = "updateZ",
  signature = "Shape",
  definition = function(object, value){
    if (nrow(object@Shape.z) != nrow(value)) {
      stop("Dimensions for updated z-scores are not consistent with 'Shape' object.")
    }
    oldValues = getOldZ(object)
    # N object stores values in matrix and oldValues in an array
    updateDim = dim(oldValues)+c(0,0,1)
    oV.updated = array(0, dim = updateDim)
    if (dim(oldValues)[3] > 0) {
      oV.updated[,,2:(dim(oV.updated)[3])] = oldValues 
    }
    oV.updated[,,1] = getZ(object)
    #Shape specific dimensionality for labels
    if (nrow(value) > 0) {
      rn = rownames(getZ(object))
    }  else {
      rn = NULL
    }
    cn = colnames(getZ(object))
    d3n = c((updateDim[3]-1):0)
    rownames(value) = rn 
    colnames(value) = cn
    dimnames = list(rn, cn, d3n)
    dimnames(oV.updated) = dimnames
    object@Shape.oldZ <- oV.updated
    object@Shape.z <- value
    return (object)
  }
)


#' Update P-values for class \linkS4class{Shape}
#' 
#' Updates Shape.sig slot from object of class \linkS4class{Shape}
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param value Matrix of 4xfpLen containing beta p-values for shape coefficients in -ddG/angstrom units.
setReplaceMethod(
  f = "updateSig",
  signature = "Shape",
  definition = function(object, value){
    if (nrow(object@Shape.sig) != nrow(value)) {
      stop("Dimensions for updated p-values are not consistent with 'Shape' object.")
    }
    oldValues = getOldSig(object)
    # N object stores values in matrix and oldValues in an array
    updateDim = dim(oldValues)+c(0,0,1)
    oV.updated = array(0, dim = updateDim)
    if (dim(oldValues)[3] > 0) {
      oV.updated[,,2:(dim(oV.updated)[3])] = oldValues 
    }
    oV.updated[,,1] = getSig(object)
    #Shape specific dimensionality for labels
    if (nrow(value) > 0) {
      rn = rownames(getSig(object))
    }  else {
      rn = NULL
    }
    cn = colnames(getSig(object))
    d3n = c((updateDim[3]-1):0)
    rownames(value) = rn 
    colnames(value) = cn
    dimnames = list(rn, cn, d3n)
    dimnames(oV.updated) = dimnames
    object@Shape.oldSig <- oV.updated
    object@Shape.sig <- value
    return (object)
  }
)






#' Add New Betas to \linkS4class{Shape} 
#' 
#' Adds new betas to \linkS4class{Shape} object using the output of a glm fit. 
#' 
#' @param object Object of class \linkS4class{Shape}
#' @param design Design table used to fit glm model. 
#' @param value Glm fit object. 
setMethod(
  f = "addNewBetas",
  signature = "Shape",
  definition = function(object, design, value, useFixedValuesOffset.Shape){
    if (length(object@shapeParamsUsed[[1]]) == 0) {
      return (object)
    }
    betaSummary = summary(value)$coefficients
    Shape.set = object@Shape.set
    Shape.notSet = c(1:object@fpLen)[!(c(1:object@fpLen) %in% object@Shape.set)]
    Shape.ncol = object@fpLen
    rowNames = row.names(getValues(object))
    numRows = nrow(getValues(object))
    Shape.values = matrix(NA, nrow = numRows, ncol = Shape.ncol,
                          dimnames = list(rowNames, 1:object@fpLen))
    Shape.errors = matrix(NA, nrow = numRows, ncol = Shape.ncol,
                          dimnames = list(rowNames, 1:object@fpLen))
    Shape.z = matrix(NA, nrow = numRows, ncol = Shape.ncol,
                     dimnames = list(rowNames, 1:object@fpLen))
    Shape.sig = matrix(NA, nrow = numRows, ncol = Shape.ncol,
                       dimnames = list(rowNames, 1:object@fpLen))
    if (length(Shape.notSet) > 0) {
      if (useFixedValuesOffset.Shape == TRUE) {
        Shape.values[,Shape.notSet] = getValues(object)[,Shape.notSet]
        Shape.errors[,Shape.notSet] = getErrors(object)[,Shape.notSet]
        Shape.z[,Shape.notSet] = getZ(object)[,Shape.notSet]
        Shape.sig[,Shape.notSet] = getSig(object)[,Shape.notSet]
      } else {
        Shape.values[,Shape.notSet] = 0
        Shape.errors[,Shape.notSet] = 0
        Shape.z[,Shape.notSet] = 0
        Shape.sig[,Shape.notSet] = 0
      }
    }
    
    if (!all(Shape.set == 0)) {
      grStr =  "^Shape.[A-Za-z]{1,}[0-9]{1,}"
      shapeParamBetas = rownames(betaSummary)[grep(grStr, rownames(betaSummary))]
      cols = as.numeric(gsub("^Shape.[A-Za-z]{1,}", "", shapeParamBetas))
      rows = gsub("[0-9]{1,}", "", shapeParamBetas)
      names(cols) = shapeParamBetas
      names(rows) = shapeParamBetas
      for (s in shapeParamBetas) {
        Shape.values[rows[s], cols[s]] = betaSummary[s,1]
        Shape.errors[rows[s], cols[s]] = betaSummary[s,2]
        Shape.z[rows[s], cols[s]] = betaSummary[s,3]
        Shape.sig[rows[s], cols[s]] = betaSummary[s,4]
      }  
      if (!(all(object@Shape.equivMat == 0))) {
        equivMat = object@Shape.equivMat
        for (i in rev(object@Shape.set)) {
          matchCol = equivMat[,i]
          matchIndices = which(matchCol != 0)
          if (length(matchIndices[matchIndices %in% object@Shape.set]) > 0) {
            matchIndex = min(matchIndices[matchIndices %in% object@Shape.set])
            if (matchCol[matchIndex] == 1) {
              S = matchIndex
              S.equiv = i
              if (S != S.equiv) {
                Shape.values[,S.equiv] = Shape.values[,S]
                Shape.errors[,S.equiv] = Shape.errors[,S]
                Shape.z[,S.equiv] = Shape.z[,S]
                Shape.sig[,S.equiv] = Shape.sig[,S]
              }
              equivMat[i,] = 0
              equivMat[i,] = 0
            } else if (matchCol[matchIndex] == -1) {
              if (matchIndex != i) {
                S = matchIndex
                S.RC= i
                shapeParams = object@shapeParamsUsed[[1]]
                if ("MGW" %in% shapeParams) {
                  Shape.values["Shape.MGW", S.RC] = Shape.values["Shape.MGW", S]
                  Shape.errors["Shape.MGW", S.RC] = Shape.errors["Shape.MGW", S]
                  Shape.z["Shape.MGW", S.RC] = Shape.z["Shape.MGW", S]
                  Shape.sig["Shape.MGW", S.RC] = Shape.sig["Shape.MGW", S]
                }
                if ("ProT" %in% shapeParams) {
                  Shape.values["Shape.ProT", S.RC] = Shape.values["Shape.ProT", S]
                  Shape.errors["Shape.ProT", S.RC] = Shape.errors["Shape.ProT", S]
                  Shape.z["Shape.ProT", S.RC] = Shape.z["Shape.ProT", S]
                  Shape.sig["Shape.ProT", S.RC] = Shape.sig["Shape.ProT", S]
                }
                if ("HelT" %in% shapeParams) {
                  Shape.values["Shape.HelTB", S.RC] = Shape.values["Shape.HelTA", S]
                  Shape.errors["Shape.HelTB", S.RC] = Shape.errors["Shape.HelTA", S]
                  Shape.z["Shape.HelTB", S.RC] = Shape.z["Shape.HelTA", S]
                  Shape.sig["Shape.HelTB", S.RC] = Shape.sig["Shape.HelTA", S]
                  Shape.values["Shape.HelTA", S.RC] = Shape.values["Shape.HelTB", S]
                  Shape.errors["Shape.HelTA", S.RC] = Shape.errors["Shape.HelTB", S]
                  Shape.z["Shape.HelTA", S.RC] = Shape.z["Shape.HelTB", S]
                  Shape.sig["Shape.HelTA", S.RC] = Shape.sig["Shape.HelTB", S]
                }
                if ("Roll" %in% shapeParams) {
                  Shape.values["Shape.RollB", S.RC] = Shape.values["Shape.RollA", S]
                  Shape.errors["Shape.RollB", S.RC] = Shape.errors["Shape.RollA", S]
                  Shape.z["Shape.RollB", S.RC] = Shape.z["Shape.RollA", S]
                  Shape.sig["Shape.RollB", S.RC] = Shape.sig["Shape.RollA", S]
                  Shape.values["Shape.RollA", S.RC] = Shape.values["Shape.RollB", S]
                  Shape.errors["Shape.RollA", S.RC] = Shape.errors["Shape.RollB", S]
                  Shape.z["Shape.RollA", S.RC] = Shape.z["Shape.RollB", S]
                  Shape.sig["Shape.RollA", S.RC] = Shape.sig["Shape.RollB", S]
                }
                equivMat[i,] = 0
              } else {
                S = i
                shapeParams = object@shapeParamsUsed[[1]]
                if ("HelT" %in% shapeParams) {
                  Shape.values["Shape.HelTB", S] = Shape.values["Shape.HelTA", S]
                  Shape.errors["Shape.HelTB", S] = Shape.errors["Shape.HelTA", S]
                  Shape.z["Shape.HelTB", S] = Shape.z["Shape.HelTA", S]
                  Shape.sig["Shape.HelTB", S] = Shape.sig["Shape.HelTA", S]
                }
                if ("Roll" %in% shapeParams) {
                  Shape.values["Shape.RollB", S] = Shape.values["Shape.RollA", S]
                  Shape.errors["Shape.RollB", S] = Shape.errors["Shape.RollA", S]
                  Shape.z["Shape.RollB", S] = Shape.z["Shape.RollA", S]
                  Shape.sig["Shape.RollB", S] = Shape.sig["Shape.RollA", S]
                }
              }
            }
          }
        }
      }
      badS = which(is.na(Shape.values[,Shape.set]), arr.ind = TRUE)
      if (length(badS) > 0) {
        if (nrow(badS) < length(Shape.values[,Shape.set])) {
          maxError = which(Shape.errors[,Shape.set] == max(Shape.errors[,Shape.set], na.rm = TRUE), arr.ind = TRUE)[1,]
          Shape.values[,Shape.set][badS] = 0
          Shape.errors[,Shape.set][badS] = 2*Shape.errors[,Shape.set][maxError["row"], maxError["col"]]
          Shape.z[badS] = 0
          Shape.sig[badS] = 1
        } else if (nrow(badS) == length(Shape.values[,Shape.set])) {
          maxError = which(Shape.errors == max(Shape.errors, na.rm = TRUE), arr.ind = TRUE)[1,]
          Shape.values[,Shape.set][badS] = 0
          Shape.values[,Shape.set][badS] = 2*Shape.errors[maxError["row"], maxError["col"]]
          Shape.z[,Shape.set][badS] = 0
          Shape.sig[,Shape.set][badS] = 1
        }
      }
    }
    
    updateValues(object) = Shape.values
    updateErrors(object) = Shape.errors
    updateZ(object)= Shape.z
    updateSig(object) = Shape.sig
    validObject(object)
    return (object)
  }
  
  
)


#' Build Null Equivalence Matrix for \linkS4class{Shape}
#' 
#' Builds equivalence matrix encoding no symmetries for an object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
setMethod(
  f = "buildNullEquivalenceMatrix",
  signature = "Shape",
  definition = function(object) {
    n = object@fpLen
    M = matrix(0, nrow = n, ncol = n, dimnames = list(c(1:n),
                                                      c(1:n)))
    return(M)
  }
) 

#' Build Reverse Complement Symmetric Equivalence Matrix for \linkS4class{Shape}
#' 
#' Builds equivalence matrix encoding reverse complement symmetry for an object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
setMethod(
  f = "buildSymmetricEquivalenceMatrix",
  signature = "Shape",
  definition = function(object) {
    n = object@fpLen
    M = matrix(0, nrow = n, ncol = n, dimnames = list(c(1:n),
                                                      c(1:n)))
    for (i in 1:(ceiling(n/2))) {
      M[i, (n+1-i)] = -1
    }
    return(M)
  }
) 

#' Formats Equivalence matrix for \linkS4class{Shape}
#' 
#' Formats an equivalence matrix for an object of class \linkS4class{Shape}, checking for conflicting instructions and getting rid of redundancies.
#' 
#' @param object Object of class \linkS4class{Shape}.
setMethod(
  f = "formatEquivalenceMatrix",
  signature = "Shape",
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

#' Get Plot Values for \linkS4class{Shape}
#' 
#' Reformats Shape.values and Shape.errors into a single data.frame for plotting.
#'  
#' @param object Object of class \linkS4class{Shape}.
#' @param iteration Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.
setMethod(
  f = "getPlotValues",
  signature = "Shape",
  definition = function(object, iteration = NULL) {
    if (is.null(iteration)) {
      values = getValues(object)
      errors = getErrors(object)
    } else {
      if (iteration == dim(object@Shape.oldValues)[3]) {
        values = getValues(object)
        errors = getErrors(object)
      } else {
        values = getOldValues(object)[,,as.character(iteration)]
        errors = getOldErrors(object)[,,as.character(iteration)]
      }
      
    }
    if (nrow(values) == 0) {
      message('No shape parameters to be plotted.')
      return ()
    }
    
    shapeLabels = rownames(values)
    shapesIncluded = gsub("Shape.", "", shapeLabels)
    
    values = t(values)
    errors = t(errors)
    values = as.data.frame(values)
    errors = as.data.frame(errors)
    colnames(values) = shapesIncluded 
    colnames(errors) = shapesIncluded 
    values$PosID = as.numeric(rownames(values))
    errors$PosID = as.numeric(rownames(errors))
    
    Values.Plot = reshape2::melt(values, id = c("PosID"))
    Values.Plot = Values.Plot[order(Values.Plot$PosID, -Values.Plot$value), ]
    rownames(Values.Plot) = NULL
    colnames(Values.Plot) = c("PosID", "Shape", "ddG")
    Errors.Plot= reshape2::melt(errors, id = c("PosID"))
    rownames(Errors.Plot) = NULL
    colnames(Errors.Plot) = c("PosID", "Shape", "SE")
    Values.Plot = merge(Values.Plot, Errors.Plot, by = c("PosID", "Shape"), sort = FALSE)
    Values.Plot$PlotErrorMax = Values.Plot$SE+Values.Plot$ddG
    Values.Plot$PlotErrorMin = Values.Plot$ddG-Values.Plot$SE
    return (Values.Plot)
    
  }
)


#' Get Number of Rows for \linkS4class{Shape} Plot
#' 
#' Compute number of rows to be used with facet_wrap for plotting beta values for objects of class 'Shape'.
#' Number of rows vary depending on the number of shape parameters fit for the object of class 'Shape'.
#' 
#' @param x Object of class \linkS4class{Shape}.
#' @param verticalPlots logical: indicates that only one shape parameter should be plotted per row. This format is ideal for markdown, where the height of the output plot can be very flexibly accommodated.
shapePlotRows = function(x, verticalPlots = FALSE) {
  sPs = nrow(x@Shape.values)
  if (verticalPlots == TRUE) {
    nR = sPs
  } else  if (sPs <= 3) {
    nR = 1
  } else if (sPs <= 5) {
    nR = 2
  } else {nR = 3}
  
  return(nR)
}

#' Get Spacing for X-labels for \linkS4class{Shape} Plot
#' 
#' Gets spacing for x-labels used for shape plotting based on the number of shape parameters per row.
#' 
#' @param x Object of class \linkS4class{Shape}.
#' @param sPs Number of rows to be used for plot(\linkS4class{Shape}) (output of \code{shapePlotRows}).
shapeXlabSpacing = function(x, sPs) {
  if (sPs >= 4) {
    step = 1
  } else if (sPs == 3) {
    if (nrow(x@Shape.values) == 3) {
      step = 1
    } else {
      step = 2
    }
  } else if (sPs == 2) {
    if (nrow(x@Shape.values) > 4) {
      step = 5
    } else if (nrow(x@Shape.values) > 2) {
      step = 2
    } else {
      step = 1
    }
  } else {
    if (nrow(x@Shape.values) > 2) {
      step = 5
    } else if (nrow(x@Shape.values) == 2) {
      step = 2
    } else {
      step = 1
    }
  }
  return (step)
}


#' Plot \linkS4class{Shape}
#' 
#' Gets shape plot from object of class \linkS4class{Shape}.
#' 
#' @param x Object of class \linkS4class{Shape}.
#' @param title Title for plot.
#' @param iter Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.  
#' @param regs Optional parameter giving labels to use for region shading. If used, lengths parameter must also be specified.
#' @param lengths Length of region parameters to be included for shading.
#' @param verticalPlots logical: If TRUE, all plots are stacked vertically. If false, plots are fit into a roughly square block. 
#' @export
setMethod(
  f = "plot",
  signature = "Shape",
  definition = function(x, title="",  iter=NULL, regs=NULL, lengths = NULL, verticalPlots = FALSE) {
    plotNumRow = shapePlotRows(x, verticalPlots)
    xlabStepSize = shapeXlabSpacing(x, plotNumRow)
    values.Plot = getPlotValues(x, iteration = iter)
    values.Plot = values.Plot[!is.na(values.Plot$ddG),]
    ddGerrs= values.Plot$SE
    ddGerrs = -1*sort(unique(-ddGerrs))
    if (length(ddGerrs) > 1) {
      if (max(ddGerrs) - 2*ddGerrs[2] < 10^(-10)) {
        values.Plot = values.Plot[!((values.Plot$ddG == 0) & (abs(values.Plot$SE- max(ddGerrs)) < 10^(-10))), ]
      }
    }
    
    nc = x@fpLen
    shapeVarList = rep(c(21, 23, 22, 24, 25), times = 5)
    shapeList = gsub("Shape.", "", rownames(x@Shape.values))
    shapeVarList = shapeVarList[1:length(shapeList)]
    names(shapeVarList) = shapeList
    yrange = max(values.Plot$PlotErrorMax)-min(values.Plot$PlotErrorMin)
    
    if ((is.null(regs)) & (is.null(lengths))) {
      shapePlot = ggplot2::ggplot(values.Plot)+
        ggplot2::geom_errorbar(data = values.Plot, ggplot2::aes(x = PosID, ymax = PlotErrorMax, ymin = PlotErrorMin, group =Shape),col = "black", width = .25)+
        #ggplot2::geom_line(ggplot2::aes(x=PosID, y=ddG, z = Shape), col = "navy")+
        ggplot2::geom_point(data = values.Plot, ggplot2::aes(x = PosID, y = ddG, shape = Shape, group = Shape), col = "navy", fill = "navy", size = 2)+
        ggplot2::facet_wrap(~ Shape, nrow = plotNumRow)+
        ggplot2::scale_x_discrete(limits = seq(0, nc, xlabStepSize)[2:((max(values.Plot$PosID) %/% xlabStepSize)+1)],
                                  labels= seq(0, nc, xlabStepSize)[2:((max(values.Plot$PosID) %/% xlabStepSize)+1)])+
        #ggplot2::scale_y_continuous(breaks = seq((min(values.Plot$PlotErrorMin) %/% .1), (max(values.Plot$PlotErrorMax) %/% .1+1)*.1, .1))+
        ggplot2::ylab(expression(paste('Shape parameters -', Delta, Delta, 'G', ' per ', ring(A))))+
        ggplot2::coord_cartesian(xlim = c(0, nc+1),
                                 ylim = c(-max(-1*min(values.Plot$PlotErrorMin), max(values.Plot$PlotErrorMax))-yMarg((yrange/20), 2),max(max(values.Plot$PlotErrorMax), -1*min(values.Plot$PlotErrorMin))+yMarg((yrange/20), 2)),
                                 expand = FALSE)+
        ggplot2::theme_bw()+
        ggplot2::scale_shape_manual(name = NULL,
                                    values = shapeVarList,
                                    breaks = shapeList,
                                    labels = shapeList)+
        ggplot2::theme(legend.position="none",
                       strip.text = ggplot2::element_text(face="bold", size = 7),
                       axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                       axis.title.y=ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                       axis.title.x=ggplot2::element_blank(),
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
      if (is.null(lengths)) {
        stop('If regs is specified, lengths for the regions must also be included.')
      } else if (is.null(regs)) {
        numReg = length(lengths[[1]])
        regs = list(as.character(1:numReg))
      }
      rects = getRegions(regs, lengths)
      shapePlot = ggplot2::ggplot(values.Plot)+
        ggplot2::geom_rect(data = rects, 
                           ggplot2::aes(xmin = xstart,
                                        xmax = xend,
                                        ymin = -Inf,
                                        ymax = Inf, 
                                        fill = col), 
                           alpha = 0.4)+
        ggplot2::geom_errorbar(data = values.Plot, ggplot2::aes(x = PosID, ymax = PlotErrorMax, ymin = PlotErrorMin, group = Shape),col = "black", width = .25)+
        #ggplot2::geom_line(aes(x=PosID, y=ddG, z = Shape), col = "navy")+
        ggplot2::geom_point(data = values.Plot, ggplot2::aes(x = PosID, y = ddG, shape = Shape, group = Shape), col = "navy", fill = "navy", size = 2)+
        ggplot2::facet_wrap(~ Shape, nrow = plotNumRow)+
        ggplot2::scale_x_discrete(limits = seq(0, max(values.Plot$PosID), xlabStepSize)[2:((max(values.Plot$PosID) %/% xlabStepSize)+1)],
                                  labels= seq(0, max(values.Plot$PosID), xlabStepSize)[2:((max(values.Plot$PosID) %/% xlabStepSize)+1)])+
        #ggplot2::scale_y_continuous(breaks = seq((min(values.Plot$PlotErrorMin) %/% .1), (max(values.Plot$PlotErrorMax) %/% .1+1)*.1, .1))+
        ggplot2::ylab(expression(paste('Shape parameters -', Delta, Delta, 'G', ' per ', ring(A))))+
        ggplot2::coord_cartesian(xlim = c(0, max(values.Plot$PosID)+1),
                                 ylim = c(-max(-1*min(values.Plot$PlotErrorMin), max(values.Plot$PlotErrorMax))-yMarg((yrange/20), 2),max(max(values.Plot$PlotErrorMax), -1*min(values.Plot$PlotErrorMin))+yMarg((yrange/20), 2)),
                                 expand = FALSE)+
        ggplot2::theme_bw()+
        ggplot2::scale_shape_manual(name = NULL,
                                    values = shapeVarList,
                                    breaks = shapeList,
                                    labels = shapeList)+
        ggplot2::theme(strip.text = ggplot2::element_text(face="bold", size = 7),
                       axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                       axis.title.y=ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                       axis.title.x=ggplot2::element_blank(),
                       plot.margin=ggplot2::unit(c(4,16,4,1),"pt"),
                       axis.ticks.length=ggplot2::unit(0.3, "mm"),
                       axis.ticks.x = ggplot2::element_line(size = .5),
                       axis.ticks.y = ggplot2::element_line(size = .5),
                       legend.text = ggplot2::element_text(size=6),
                       axis.text.x = ggplot2::element_text(size = 6,
                                                           vjust = -1),
                       plot.title = ggplot2::element_text(size = 8, 
                                                          face = "bold"),
                       legend.margin=ggplot2::unit(.01,"cm"), 
                       legend.position = "bottom",
                       legend.direction = "horizontal", 
                       legend.box= "horizontal",
                       legend.title = ggplot2::element_blank(),
                       legend.background = ggplot2::element_rect(
                         fill = "transparent", 
                         size = 1))+
        ggplot2::ggtitle(title)
      if (length(levels(rects$col)) <= 7) {
        fillVals = rep(character(0), times = length(levels(rects$col)))
        for (i in 1:length(levels(rects$col))) {
          fillVals[i] = rects$hexCol[rects$col == levels(rects$col)[i]][1]
        }
        names(fillVals) = levels(rects$col)
        
        shapePlot <- shapePlot +ggplot2::scale_fill_manual(name = NULL,
                                                           values= fillVals,
                                                           breaks = levels(rects$col),
                                                           labels = levels(rects$col))+
          ggplot2::guides(position = "bottom",
                          colour = FALSE,
                          shape = FALSE,
                          fill=ggplot2::guide_legend(order = 2,
                                                     keywidth=0.6,
                                                     keyheight=0.4,
                                                     default.unit="cm"))
      }
    }
    
    return(shapePlot)
  }
)



#' Plot Z-scores for \linkS4class{Shape}
#' 
#' Gets shape plot from object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param title Title for plot.
#' @param maxErr Maximum height of error bars to use in plotting.
#' @param iter Iteration of beta values to be used in plot. If iteration == NULL, the most recent round of beta values are used.  
#' @param regs Optional parameter giving labels to use for region shading. If used, lengths parameter must also be specified.
#' @param lengths Length of region parameters to be included.
#' @param verticalPlots logical: If TRUE, all plots are stacked vertically. If false, plots are fit into a roughly square block. 
#' @export
setMethod(
  f = "plotZs",
  signature = "Shape",
  definition = function(object, title="", iter=NULL, regs=NULL, lengths = NULL, verticalPlots=FALSE) {
    plotNumRow = shapePlotRows(x, verticalPlots)
    xlabStepSize = shapeXlabSpacing(x, plotNumRow)
    values.Plot = getPlotValues(object, iter)
    ddGerrs= values.Plot$SE
    ddGerrs = -1*sort(unique(-ddGerrs))
    if (length(ddGerrs) > 1) {
      if (max(ddGerrs) - 2*ddGerrs[2] < 10^(-10)) {
        values.Plot = values.Plot[!((values.Plot$ddG == 0) & (abs(values.Plot$SE- max(ddGerrs)) < 10^(-10))), ]
      }
    }
    if (all(values.Plot$SE == 0)) {
      stop('Cannot plot z-scores equal to infinity.')
    }
    values.Plot$Z = values.Plot$ddG/values.Plot$SE
    nc = object@fpLen
    shapeVarList = rep(c(21, 23, 22, 24, 25), times = 5)
    shapeList = gsub("Shape.", "", rownames(object@Shape.values))
    shapeVarList = shapeVarList[1:length(shapeList)]
    names(shapeVarList) = shapeList
    yrange = max(values.Plot$Z)-min(values.Plot$Z)
    if ((is.null(regs)) & (is.null(lengths))) {
      zPlot = ggplot2::ggplot(values.Plot)+
        ggplot2::geom_point(data = values.Plot, ggplot2::aes(x = PosID, y = Z, shape = Shape, group = Shape), col = "navy", fill = "navy", size = 2)+
        ggplot2::facet_wrap(~ Shape, nrow = plotNumRow)+
        ggplot2::scale_x_discrete(limits = seq(0, max(values.Plot$PosID), xlabStepSize)[2:((max(values.Plot$PosID) %/% xlabStepSize)+1)],
                                  labels= seq(0, max(values.Plot$PosID), xlabStepSize)[2:((max(values.Plot$PosID) %/% xlabStepSize)+1)])+
        #ggplot2::scale_y_continuous(breaks = seq((min(values.Plot$Z) %/% .1), (max(values.Plot$Z) %/% .1+1)*.1, .1))+
        ggplot2::ylab(expression('Z'[paste('-', Delta, Delta, 'G', ' per ', ring(A))]))+
        ggplot2::coord_cartesian(xlim = c(0, max(values.Plot$PosID)+1),
                                 ylim = c(-max(-1*min(values.Plot$Z), max(values.Plot$Z))-yMarg((yrange/20), 2),max(max(values.Plot$Z), -1*min(values.Plot$Z))+yMarg((yrange/20), 2)),
                                 expand = FALSE)+
        ggplot2::theme_bw()+
        ggplot2::scale_shape_manual(name = NULL,
                                    values = shapeVarList,
                                    breaks = shapeList,
                                    labels = shapeList)+
        
        ggplot2::theme(legend.position="none",
                       strip.text = ggplot2::element_text(face="bold", size = 7),
                       axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                       axis.title.y=ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                       axis.title.x=ggplot2::element_blank(),
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
      if (is.null(lengths)) {
        stop('If regs is specified, lengths for the regions must also be included.')
      } else if (is.null(regs)) {
        numReg = length(lengths[[1]])
        regs = list(as.character(1:numReg))
      }
      
      rects = getRegions(regs, lengths)
      zPlot = ggplot2::ggplot(values.Plot)+
        ggplot2::geom_rect(data = rects, 
                           ggplot2::aes(xmin = xstart,
                                        xmax = xend,
                                        ymin = -Inf,
                                        ymax = Inf, 
                                        fill = col), 
                           alpha = 0.4)+
        ggplot2::geom_point(data = values.Plot, ggplot2::aes(x = PosID, y = Z, shape = Shape, group = Shape), col = "navy", fill = "navy", size = 2)+
        ggplot2::facet_wrap(~ Shape, nrow = plotNumRow)+
        ggplot2::scale_x_discrete(limits = seq(0, max(values.Plot$PosID), xlabStepSize)[2:((max(values.Plot$PosID) %/% xlabStepSize)+1)],
                                  labels= seq(0, max(values.Plot$PosID), xlabStepSize)[2:((max(values.Plot$PosID) %/% xlabStepSize)+1)])+
        #ggplot2::scale_y_continuous(breaks = seq((min(values.Plot$Z) %/% .1), (max(values.Plot$Px) %/% .1+1)*.1, .1))+
        ggplot2::ylab(expression('Z'[paste('-', Delta, Delta, 'G', ' per ', ring(A))]))+
        ggplot2::coord_cartesian(xlim = c(0, max(values.Plot$PosID)+1),
                                 ylim = c(-max(-1*min(values.Plot$Z), max(values.Plot$Z))-yMarg((yrange/20), 2),max(max(values.Plot$Z), -1*min(values.Plot$Z))+yMarg((yrange/20), 2)),
                                 expand = FALSE)+
        ggplot2::theme_bw()+
        ggplot2::scale_shape_manual(name = NULL,
                                    values = shapeVarList,
                                    breaks = shapeList,
                                    labels = shapeList)+
        ggplot2::theme(strip.text = ggplot2::element_text(face="bold", size = 7),
                       plot.margin=ggplot2::unit(c(4,16,2,1),"pt"),
                       axis.text.y=ggplot2::element_text(size = 6, color = "black", angle = 45, hjust = 1, vjust = 1),
                       axis.title.y=ggplot2::element_text(size = 7, color = "black", margin = ggplot2::margin(b = 2, t = 3)),
                       axis.title.x=ggplot2::element_blank(),
                       axis.ticks.length=ggplot2::unit(0.3, "mm"),
                       axis.ticks.x = ggplot2::element_line(size = .5),
                       axis.ticks.y = ggplot2::element_line(size = .5),
                       legend.text = ggplot2::element_text(size=6),
                       axis.text.x = ggplot2::element_text(size = 6,
                                                           vjust = -1),
                       plot.title = ggplot2::element_text(size = 8, 
                                                          face = "bold"),
                       legend.margin=ggplot2::unit(-.05,"cm"), 
                       legend.position = "bottom",
                       legend.direction = "horizontal", 
                       legend.box= "horizontal",
                       legend.title = ggplot2::element_blank(),
                       legend.background = ggplot2::element_rect(
                         fill = "transparent", 
                         size = 1))+
        ggplot2::ggtitle(title)
      
      if (length(levels(rects$col)) <= 7) {
        fillVals = rep(character(0), times = length(levels(rects$col)))
        for (i in 1:length(levels(rects$col))) {
          fillVals[i] = rects$hexCol[rects$col == levels(rects$col)[i]][1]
        }
        names(fillVals) = levels(rects$col)
        
        zPlot <- zPlot +ggplot2::scale_fill_manual(name = NULL,
                                                   values= fillVals,
                                                   breaks = levels(rects$col),
                                                   labels = levels(rects$col))+
          ggplot2::guides(position = "bottom",
                          colour = FALSE,
                          shape = FALSE,
                          fill=ggplot2::guide_legend(order = 2,
                                                     keywidth=0.6,
                                                     keyheight=0.4,
                                                     default.unit="cm"))
      }
    }
    
    return(zPlot)
  }
)

#' Design Matrix for \linkS4class{Shape} Features
#' 
#' Gets shape design from design matrix, i.e. gets counts for each mono-nucleotide feature.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @param design Data frame with design matrix, as generated by \code{addDesignMatrix} function.
#' @export
setMethod(
  f = "getDesignMatrix",
  signature = "Shape",
  definition = function(object, design) {
    if (all(object@Shape.set == 0)) {
      mat = matrix(0, nrow = 0, ncol = 0)
      message('No shape parameters included in fit.')
      return (mat)
    }
    S.set = object@Shape.set
    S.newSet = c()
    shapes = rownames(object@Shape.values)
    if (!(all(object@Shape.equivMat == 0))) {
      equivMat = object@Shape.equivMat
      for (i in rev(S.set)) {
        matchCol = equivMat[,i]
        matchIndices = which(matchCol != 0)
        if (length(matchIndices[matchIndices %in% object@Shape.set]) > 0) {
          matchIndex = min(matchIndices[matchIndices %in% object@Shape.set])
          equivMat[i,] = 0
        } else {
          S.newSet = c(S.newSet, i)
        }
      }
      S.newSet = unique(S.newSet[order(S.newSet)])
    } else {
      S.newSet = S.set
    }
    
    mat = matrix(0, nrow=length(S.newSet), ncol=length(shapes), dimnames = list(S.newSet, shapes))
    
    for (i in S.newSet) {
      for (s in shapes) {
        r = as.character(i)
        mat[r,s] = mean(design[,paste(s, r, sep = "")]/design$Round)
      }
    }
    
    
    return(mat)
  }
)


