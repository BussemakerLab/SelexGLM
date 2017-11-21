#' @include featureSet-validity.R
#' @include N-class.R
#' @include Intercept-class.R
#' @include Shape-class.R
#' @include generic-functions-definitions.R
NULL

#' Class 'featureSet' Definition
#' 
#' Defines an S4 class to represent features,
#' 
#' @slot seedLen number of base pairs in the seeding model.
#' @slot upFootprintExtend maximum number of upstream positions to be fit beyond the footprint for any feature.
#' @slot downFootprintExtend maximum number of downstream positions to be fit beyond the footprint for any feature. 
#' @slot numViews number of views scored on each DNA strand.
#' @slot rounds rounds of SELEX data used to evaluate 'model' to which 'Intercept' belongs.
#' @slot rcSymmetric logical: indicating whether reverse complement symmetric model fit should be used. 
#' @slot shapeParamsUsed list of shape features to be included in 'model'. 
#' @slot N object of class 'N'
#' @slot Intercept object of class 'Intercept'
#' @slot Shape object of class 'Shape'
#' @slot glmFits List of glm fit summaries.
#' @slot designMatrixSummary List of design matrices from getDesignMatrix corresponding to the glmFits.
#' @export featureSet
featureSet <- setClass("featureSet", 
                       representation(seedLen = "numeric",
                                      upFootprintExtend = "numeric", 
                                      downFootprintExtend = "numeric",
                                      numViews = "numeric",
                                      rounds = "list",
                                      shapeParamsUsed = "list",
                                      rcSymmetric = "logical",
                                      N = "N",
                                      Intercept = "Intercept",
                                      Shape = "Shape",
                                      glmFits = "list",
                                      designMatrixSummary = "list"), 
                       validity = validFeatureSet)

#' Initialize object of class \linkS4class{featureSet}.
#' 
#' Initializes an S4 class to represent full model features.
#' 
#' @param featureSet Object of class \linkS4class{featureSet}.
#' @export 
setMethod("initialize", 
          "featureSet", 
          function(.Object,
                   seedLen,
                   upFootprintExtend,
                   downFootprintExtend,
                   rcSymmetric,
                   glmFits,
                   designMatrixSummary,
                   N,
                   N.upFootprintExtend, 
                   N.downFootprintExtend,
                   N.set, 
                   N.values,
                   N.errors, 
                   N.z, 
                   N.sig,
                   N.equivMat,
                   N.oldValues,
                   N.oldErrors,
                   N.oldZ,
                   N.oldSig,
                   Intercept,
                   numViews,
                   rounds,
                   I.values, 
                   I.errors, 
                   I.z, 
                   I.sig,
                   I.oldValues,
                   I.oldErrors,
                   I.oldZ,
                   I.oldSig,
                   Shape,
                   shapeParamsUsed, 
                   Shape.values, 
                   Shape.errors,
                   Shape.z, 
                   Shape.sig,
                   Shape.oldValues,
                   Shape.oldErrors,
                   Shape.oldZ,
                   Shape.oldSig,
                   Shape.set, 
                   Shape.upFootprintExtend,
                   Shape.downFootprintExtend, 
                   Shape.equivMat,
                   ...) {
            
            # initialize seedLen slot
            if (missing(seedLen)) {
              stop("Must specify a value for 'seedLen'.")
            } else {
              .Object@seedLen = seedLen
            }
            
            # initialize upFootprintExtend
            if (missing(upFootprintExtend)) {
              stop("Must specify a value for 'upFootprintExtend'.")
            } else {
              .Object@upFootprintExtend = upFootprintExtend
            }
            
            # initialize downFootprintExtend
            if (missing(downFootprintExtend)) {
              .Object@downFootprintExtend = .Object@upFootprintExtend
            } else {
              .Object@downFootprintExtend = downFootprintExtend
            }
            # initialize glmFits slot
            if (missing(glmFits)) {
              .Object@glmFits = list()
            } else {
              .Object@glmFits = glmFits
            }
            
            # initialize designMatrixSummary slots
            if (missing(designMatrixSummary)) {
              .Object@designMatrixSummary = list()
            } else {
              .Object@designMatrixSummary = designMatrixSummary
            }
            # initialize rcSymmetric
            if (missing(rcSymmetric)) {
              .Object@rcSymmetric = FALSE
            } else {
              .Object@rcSymmetric = rcSymmetric
            }
            
            # initialize numViews 
            if (missing(numViews)) {
              if (missing(Intercept)) {
                if (missing(I.values)) {
                  if (missing(I.oldValues)) {
                    stop("Must either specify a value for 'numViews' or specify 'numViews' implicitly through the dimensions of 'I' or its slots 'I.values' or 'I.oldValues'.")
                  } else {
                    .Object@numViews = ncol(I.oldValues)
                  }
                } else {
                  .Object@numViews = ncol(I.values)
                }
              } else {
                .Object@numViews = ncol(getValues(Intercept))
              }
            } else {
              .Object@numViews = numViews
            }
            
            # initialize rounds
            if (missing(rounds)) {
              stop("Must specify a value for 'rounds'.")
            } else {
              rounds[[1]] = rounds[[1]][order(rounds[[1]])]
              .Object@rounds = rounds
            }
            
            # initialize shapeParamsUsed 
            if (missing(shapeParamsUsed)) {
              if (missing(Shape)) {
                if (missing(Shape.values)) {
                  if (missing(Shape.oldValues)) {
                    if (missing(Shape.set)) {
                      stop("Must either specify a value for 'shapeParamsUsed' or specify 'shapeParamsUsed' implicitly through the dimension names of 'Shape' or its slots 'Shape.values', 'Shape.oldValues', or 'Shape.set' (only in the case of 'Shape.set' = c(0)).")
                    } else {
                      if (all(Shape.set == 0)) {
                        .Object@shapeParamsUsed = list(c(character(0)))
                      } else {
                        stop("Must either specify a value for 'shapeParamsUsed' or specify 'shapeParamsUsed' implicitly through the dimension names of 'Shape' or its slots 'Shape.values', 'Shape.oldValues', or 'Shape.set' (only in the case of 'Shape.set' = c(0)).")
                      }
                    }
                  } else {
                    if (nrow(Shape.oldValues) == 0) {
                      .Object@shapeParamsUsed = list(c(character(0)))
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
                      .Object@shapeParamsUsed = list(shapeParams)
                      .Object@shapeParamsUsed[[1]] = unique(.Object@shapeParamsUsed[[1]][order(.Object@shapeParamsUsed[[1]])])
                    }
                  }
                } else {
                  if (nrow(Shape.values) == 0) {
                    .Object@shapeParamsUsed = list(c(character(0)))
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
                    .Object@shapeParamsUsed = list(shapeParams)
                    .Object@shapeParamsUsed[[1]] = unique(.Object@shapeParamsUsed[[1]][order(.Object@shapeParamsUsed[[1]])])
                  }
                }
              } else {
                .Object@shapeParamsUsed = Shape@shapeParamsUsed
                .Object@shapeParamsUsed[[1]] = unique(.Object@shapeParamsUsed[[1]][order(.Object@shapeParamsUsed[[1]])])
              }
            } else {
              .Object@shapeParamsUsed = shapeParamsUsed
              .Object@shapeParamsUsed[[1]] = unique(.Object@shapeParamsUsed[[1]][order(.Object@shapeParamsUsed[[1]])])
              
              
            }
            .Object@shapeParamsUsed[[1]] = unique(.Object@shapeParamsUsed[[1]][order(.Object@shapeParamsUsed[[1]])])
            # initialize N
            if (!missing(N)) {
              if (!(missing(N.upFootprintExtend))) {
                stop ("If 'N' is specified, the 'N.upFootprintExtend' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.downFootprintExtend))) {
                stop ("If 'N' is specified, the 'N.downFootprintExtend' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.set))) {
                stop ("If 'N' is specified, the 'N.set' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.equivMat))) {
                stop ("If 'N' is specified, the 'N.equivMat' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.values))) {
                stop ("If 'N' is specified, the 'N.values' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.errors))) {
                stop ("If 'N' is specified, the 'N.errors' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.z))) {
                stop ("If 'N' is specified, the 'N.z' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.sig))) {
                stop ("If 'N' is specified, the 'N.sig' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.oldValues))) {
                stop ("If 'N' is specified, the 'N.oldValues' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.oldErrors))) {
                stop ("If 'N' is specified, the 'N.oldErrors' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.oldZ))) {
                stop ("If 'N' is specified, the 'N.oldZ' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              if (!(missing(N.oldSig))) {
                stop ("If 'N' is specified, the 'N.oldSig' slot of 'N' cannot be specified outside of the class 'N' object.")
              }
              .Object@N <- N(seedLen = N@seedLen,
                             N.upFootprintExtend = N@N.upFootprintExtend,
                             N.downFootprintExtend = N@N.downFootprintExtend, 
                             fS.upFootprintExtend = N@fS.upFootprintExtend,
                             fS.downFootprintExtend = N@fS.downFootprintExtend,
                             N.set = N@N.set,
                             N.equivMat = N@N.equivMat,
                             N.values = N@N.values,
                             N.errors= N@N.errors,
                             N.z= N@N.z,
                             N.sig = N@N.sig,
                             N.oldValues = N@N.oldValues,
                             N.oldErrors = N@N.oldErrors,
                             N.oldZ = N@N.oldZ,
                             N.oldSig = N@N.oldSig)
            } else {
              
              N.args = list("N", seedLen = .Object@seedLen, fS.upFootprintExtend = .Object@upFootprintExtend,
                            fS.downFootprintExtend = .Object@downFootprintExtend, rcSymmetric = .Object@rcSymmetric)
              
              # initialize N.upFootprintExtend
              if (!(missing(N.upFootprintExtend))) {
                N.args = c(N.args, N.upFootprintExtend = N.upFootprintExtend)
              }
              # initialize N.downFootprintExtend
              if (!(missing(N.downFootprintExtend))) {
                N.args = c(N.args, N.downFootprintExtend = N.downFootprintExtend)
              }
              
              # initialize N.set
              if (!(missing(N.set))) {
                N.args = append(N.args, list(N.set = N.set))
              }
              
              # initialize N.equivMat
              if (!(missing(N.equivMat))) {
                N.args = append(N.args, list(N.equivMat=N.equivMat))
              }
              
              # initialize N.values
              if (!(missing(N.values))) {
                N.args = append(N.args, list(N.values = N.values))
              }
              
              # initialize N.errors
              if (!(missing(N.errors))) {
                N.args = append(N.args, list(N.errors = N.errors))
              }
              
              # initialize N.z
              if (!(missing(N.z))) {
                N.args = append(N.args, list(N.z = N.z))
              }
              
              # initialize N.sig
              if (!(missing(N.sig))) {
                N.args = append(N.args, list(N.sig = N.sig))
              }
              
              # initialize N.oldValues
              if (!(missing(N.oldValues))) {
                N.args = append(N.args, list(N.oldValues = N.oldValues))
              }
              
              # initialize N.oldErrors
              if (!(missing(N.oldErrors))) {
                N.args = append(N.args, list(N.oldErrors = N.oldErrors))
              }
              
              # initialize N.oldZ
              if (!(missing(N.oldZ))) {
                N.args = append(N.args, list(N.oldZ = N.oldZ))
              }
              
              # initialize N.oldSig
              if (!(missing(N.oldSig))) {
                N.args = append(N.args, list(N.oldSig = N.oldSig))
              }
              .Object@N = do.call("new", N.args)
            } 

    
            # initialize Intercept
            if (!missing(Intercept)) {
              if (!(missing(I.values))) {
                stop ("If 'Intercept' is specified, the 'I.values' slot of 'Intercept' cannot be specified outside of the class 'Intercept' object.")
              }
              if (!(missing(I.errors))) {
                stop ("If 'Intercept' is specified, the 'I.errors' slot of 'Intercept' cannot be specified outside of the class 'Intercept' object.")
              }
              if (!(missing(I.z))) {
                stop ("If 'Intercept' is specified, the 'I.z' slot of 'Intercept' cannot be specified outside of the class 'Intercept' object.")
              }
              if (!(missing(I.sig))) {
                stop ("If 'Intercept' is specified, the 'I.sig' slot of 'Intercept' cannot be specified outside of the class 'Intercept' object.")
              }
              if (!(missing(I.oldValues))) {
                stop ("If 'Intercept' is specified, the 'I.oldValues' slot of 'Intercept' cannot be specified outside of the class 'Intercept' object.")
              }
              if (!(missing(I.oldErrors))) {
                stop ("If 'Intercept' is specified, the 'I.oldErrors' slot of 'Intercept' cannot be specified outside of the class 'Intercept' object.")
              }
              if (!(missing(I.oldZ))) {
                stop ("If 'Intercept' is specified, the 'I.oldZ' slot of 'Intercept' cannot be specified outside of the class 'Intercept' object.")
              }
              if (!(missing(I.oldSig))) {
                stop ("If 'Intercept' is specified, the 'I.oldSig' slot of 'Intercept' cannot be specified outside of the class 'Intercept' object.")
              }
              .Object@Intercept <- Intercept(rounds = list(gsub("Round.", "", dimnames(Intercept@I.values)[[3]])),
                                         I.values = Intercept@I.values,
                                         I.errors= Intercept@I.errors,
                                         I.z= Intercept@I.z,
                                         I.sig = Intercept@I.sig,
                                         I.oldValues = Intercept@I.oldValues,
                                         I.oldErrors = Intercept@I.oldErrors,
                                         I.oldZ = Intercept@I.oldZ,
                                         I.oldSig = Intercept@I.oldSig)
            } else {
              I.args = list("Intercept", rounds = .Object@rounds)
              if (!missing(numViews)) {
                I.args = c(I.args, numViews = numViews)
              }
              
              # initialize I.values
              if (!(missing(I.values))) {
                I.args = append(I.args, list(I.values = I.values))
              }
              
              # initialize I.errors
              if (!(missing(I.errors))) {
                I.args = append(I.args, list(I.errors = I.errors))
              }
              
              # initialize I.z
              if (!(missing(I.z))) {
                I.args = append(I.args, list(I.z = I.z))
              }
              
              # initialize I.sig
              if (!(missing(I.sig))) {
                I.args = append(I.args, list(I.sig = I.sig))
              }
              
              # initialize I.oldValues
              if (!(missing(I.oldValues))) {
                I.args = append(I.args, list(I.oldValues = I.oldValues))
              }
              
              # initialize I.oldErrors
              if (!(missing(I.oldErrors))) {
                I.args = append(I.args, list(I.oldErrors = I.oldErrors))
              }
              
              # initialize I.oldZ
              if (!(missing(I.oldZ))) {
                I.args = append(I.args, list(I.oldZ = I.oldZ))
              }
              
              # initialize I.oldSig
              if (!(missing(I.oldSig))) {
                I.args = append(I.args, list(I.oldSig = I.oldSig))
              }
              .Object@Intercept = do.call("new", I.args)
            }
            
            # initialize Shape
            if (!missing(Shape)) {
              if (!(missing(Shape.upFootprintExtend))) {
                stop ("If 'Shape' is specified, the 'Shape.upFootprintExtend' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.downFootprintExtend))) {
                stop ("If 'Shape' is specified, the 'Shape.downFootprintExtend' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.set))) {
                stop ("If 'Shape' is specified, the 'Shape.set' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(shapeParamsUsed))) {
                stop ("If 'Shape' is specified, the 'shapeParamsUsed' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.equivMat))) {
                stop ("If 'Shape' is specified, the 'Shape.equivMat' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.values))) {
                stop ("If 'Shape' is specified, the 'Shape.values' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.errors))) {
                stop ("If 'Shape' is specified, the 'Shape.errors' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.z))) {
                stop ("If 'Shape' is specified, the 'Shape.z' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.sig))) {
                stop ("If 'Shape' is specified, the 'Shape.sig' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.oldValues))) {
                stop ("If 'Shape' is specified, the 'Shape.oldValues' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.oldErrors))) {
                stop ("If 'Shape' is specified, the 'Shape.oldErrors' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.oldZ))) {
                stop ("If 'Shape' is specified, the 'Shape.oldZ' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              if (!(missing(Shape.oldSig))) {
                stop ("If 'Shape' is specified, the 'Shape.oldSig' slot of 'Shape' cannot be specified outside of the class 'Shape' object.")
              }
              .Object@Shape <-  Shape(seedLen = Shape@seedLen,
                                  Shape.upFootprintExtend = Shape@Shape.upFootprintExtend,
                                  Shape.downFootprintExtend = Shape@Shape.downFootprintExtend, 
                                  fS.upFootprintExtend = Shape@fS.upFootprintExtend,
                                  fS.downFootprintExtend = Shape@fS.downFootprintExtend,
                                  Shape.set = Shape@Shape.set,
                                  shapeParamsUsed = Shape@shapeParamsUsed,
                                  Shape.equivMat = Shape@Shape.equivMat,
                                  Shape.values = Shape@Shape.values,
                                  Shape.errors= Shape@Shape.errors,
                                  Shape.z= Shape@Shape.z,
                                  Shape.sig = Shape@Shape.sig,
                                  Shape.oldValues = Shape@Shape.oldValues,
                                  Shape.oldErrors = Shape@Shape.oldErrors,
                                  Shape.oldZ = Shape@Shape.oldZ,
                                  Shape.oldSig = Shape@Shape.oldSig)
            } else {
              Shape.args = list("Shape",
                                seedLen = .Object@seedLen,
                                fS.upFootprintExtend = .Object@upFootprintExtend,
                                fS.downFootprintExtend = .Object@downFootprintExtend,
                                rcSymmetric = .Object@rcSymmetric)
              
              # initialize Shape.upFootprintExtend
              Shape.args = append(Shape.args, list(shapeParamsUsed = shapeParamsUsed))
              # initialize Shape.upFootprintExtend
              if (!(missing(Shape.upFootprintExtend))) {
                Shape.args = append(Shape.args, list(Shape.upFootprintExtend = Shape.upFootprintExtend))
              }
              
              # initialize Shape.downFootprintExtend
              if (!(missing(Shape.downFootprintExtend))) {
                Shape.args = append(Shape.args, list(Shape.downFootprintExtend = Shape.downFootprintExtend))
              }
              
              # initialize Shape.set
              if (!(missing(Shape.set))) {
                Shape.args = append(Shape.args, list(Shape.set = Shape.set))
              }
              
              # initialize Shape.equivMat
              if (!(missing(Shape.equivMat))) {
                Shape.args = append(Shape.args, list(Shape.equivMat = Shape.equivMat))
              }
              
              # initialize Shape.values
              if (!(missing(Shape.values))) {
                Shape.args = append(Shape.args, list(Shape.values = Shape.values))
              }
              
              # initialize Shape.errors
              if (!(missing(Shape.errors))) {
                Shape.args = append(Shape.args, list(Shape.errors = Shape.errors))
              }
              
              # initialize Shape.z
              if (!(missing(Shape.z))) {
                Shape.args = append(Shape.args, list(Shape.z = Shape.z))
              }
              
              # initialize Shape.sig
              if (!(missing(Shape.sig))) {
                Shape.args = c(Shape.args, Shape.sig = Shape.sig)
              }
              
              # initialize Shape.oldValues
              if (!(missing(Shape.oldValues))) {
                Shape.args = append(Shape.args, list(Shape.oldValues = Shape.oldValues))
              }
              
              # initialize Shape.oldErrors
              if (!(missing(Shape.oldErrors))) {
                Shape.args = append(Shape.args, list(Shape.oldErrors = Shape.oldErrors))
                
              }
              
              # initialize Shape.oldZ
              if (!(missing(Shape.oldZ))) {
                Shape.args = append(Shape.args, list(Shape.oldZ = Shape.oldZ))
              
              }
              
              # initialize Shape.oldSig
              if (!(missing(Shape.oldSig))) {
                Shape.args = append(Shape.args, list(Shape.oldSig = Shape.oldSig))
             
              }
              .Object@Shape = do.call("new", Shape.args)
            } 
            
            validObject(.Object)
            .Object
          })



#' Show method for class 'featureSet'
#'
#' Defines a show method for class 'featureSet'.
#'
#' @param object object of class 'featureSet'.
#' @export
setMethod(
  f = "show",
  signature = "featureSet",
  definition = function(object) {
    cat('An object of class "', class(object), '"\n', sep = "")
    cat("\n")
    cat('Slot "seedLen": ', object@seedLen, '\n') 
    cat("\n")
    cat('Slot "upFootprintExtend": ', object@upFootprintExtend, '\n') 
    cat("\n")
    cat('Slot "downFootprintExtend": ', object@downFootprintExtend, '\n')
    cat("\n")
    cat('Slot "numViews": ', object@numViews, '\n')
    cat("\n")
    cat('Slot "rounds": ', object@rounds[[1]], '\n')
    cat("\n")
    cat('Slot "shapeParamsUsed": ', object@shapeParamsUsed[[1]], '\n')
    cat("\n")
    cat('Slot "rcSymmetric": ', object@rcSymmetric, '\n')
    cat("\n")
    cat('Slot "N": \n') 
    show(object@N)
    cat("\n")
    cat('Slot "Intercept": \n') 
    show(object@Intercept)
    cat("\n")
    cat('Slot "Shape": \n') 
    show(object@Shape)
    invisible(NULL)
  }
)


#' Print method for \linkS4class{featureSet}
#' 
#' Defines a print method for \linkS4class{featureSet}.
#' 
#' @param x Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "print",
  signature = "featureSet",
  definition = function(x) {
    show(x)
  }
)

#' Get Feature Design for \linkS4class{featureSet}
#' 
#' Defines a print design method for class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getFeatureDesign",
  signature = "featureSet",
  definition = function(object) {
    cat("Feature design for object of class '", class(object), "'\n", sep = "")
    cat("\n")
    cat('seedLen: ', object@seedLen, '\n')
    cat('upFootprintExtend: ', object@upFootprintExtend, '\n')
    cat('downFootprintExtend: ', object@downFootprintExtend, '\n')
    cat('rcSymmetric: ', object@rcSymmetric, '\n')
    cat("\n")
    cat('Slot "N":', '\n')
    cat('N.upFootprintExtend: ', object@N@N.upFootprintExtend, '\n')
    cat('N.downFootprintExtend: ', object@N@N.downFootprintExtend, '\n')
    cat('N.set: ', object@N@N.set, '\n')
    cat('Number of previous iterations: ', dim(object@N@N.oldValues)[3], '\n')
    cat("\n")
    cat('Slot "Intercept":', '\n')
    cat("Number of Views per Strand of DNA: ", ncol(object@Intercept@I.values), '\n', sep = "")
    cat("Number of Rounds: ", dim(object@Intercept@I.values)[3]," (", paste(gsub("Round.", "", dimnames(object@Intercept@I.values)[[3]]), collapse = ", "),")\n", sep = "")
    cat('Number of previous iterations: ', dim(object@Intercept@I.oldValues)[4], '\n', sep = "")
    cat("\n")
    cat('Slot "Shape":', '\n')
    if (nrow(object@Shape@Shape.values) > 0) {
      cat('"ShapeParamsUsed": ', object@Shape@shapeParamsUsed[[1]])
      cat('Shape.upFootprintExtend: ', object@Shape@Shape.upFootprintExtend, '\n')
      cat('Shape.downFootprintExtend: ', object@Shape@Shape.downFootprintExtend, '\n')
      cat('Shape.set: ', object@Shape@Shape.set, '\n')
      cat('Number of previous iterations: ', dim(object@Shape@Shape.oldValues)[3], '\n')
    } else {
      cat('"ShapeParamsUsed": NONE')
    }
    cat("\n")
    invisible(NULL)
  }
)


#' Summary Method for \linkS4class{featureSet}
#' 
#' Defines a summary method for class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "summary",
  signature = "featureSet",
  definition = function(object) {
    cat("An object of class ", class(object), "\n", sep = "")
    cat("\n")
    shapeThing = ""
    if (length(object@shapeParamsUsed[[1]]) > 0) {
      shapeThing = paste(length(object@shapeParamsUsed[[1]]),
            " shape parameter(s) (",
            paste(object@shapeParamsUsed[[1]], collapse = ", "),
            ") at ",
            length(object@Shape@Shape.set),
            " positions (Shape.set = ",
            paste(object@Shape@Shape.set, collapse = ", "),
            "), ",
            sep = "")
    }
    
    cat("Fits a feature model of length ",
        length(object@seedLen+object@upFootprintExtend+object@downFootprintExtend),
        " with mononucleotides at ",
        length(object@N@N.set),
        " positions (N.set = ",
        paste(object@N@N.set, collapse = ", "),
        "), ",
        shapeThing,
        object@numViews,
        " view(s) per strand of DNA, ",
        length(object@rounds[[1]]),
        " round(s) of data, and reverse complement symmetry = ",
        object@rcSymmetric,
        ".\n",
        sep = "")
    cat("Includes the following feature sub-classes: \n")
    cat("\n")
    cat("An object of class ", class(object@N), "\n", sep = "")
    cat("Fits ",
        length(object@N@N.set),
        " nucleotides at positions ",
        paste(object@N@N.set, collapse = ", "),
        " for a feature model of length ",
        (object@N@fS.upFootprintExtend+object@N@fS.downFootprintExtend+object@N@seedLen),
        ".\n",
        sep = "")
    if (!all(object@N@N.equivMat == 0)) {
      if (all(object@N@N.equivMat[,1:(object@N@fpLen %/%2)] == 0)) {
        halfMatCheck = object@N@N.equivMat[1:(object@N@fpLen %/%2+object@N@fpLen %%2),
                                           (object@N@fpLen %/%2+1):object@N@fpLen]
        halfMatCheck = halfMatCheck[,ncol(halfMatCheck):1]
        if (all(diag(halfMatCheck) == -1)) {
          cat("Nucleotide features are reverse complement symmetric.\n")
        }
      }
    }
    cat("\n")
    cat("Nucleotide beta values:\n")
    print(getValues(object@N))
    cat("\n")
    cat("Nucleotide beta errors:\n")
    print(getErrors(object@N))
    cat("\n")
    cat("Nucleotide beta z-scores:\n")
    print(getZ(object@N))
    cat("\n")
    cat("Nucleotide beta p-values:\n")
    print(getSig(object@N))
    cat("\n")
    cat("\n")
    cat("An object of class ", class(object@Intercept), "\n", sep = "")
    cat("Fits up to ",
        2*ncol(object@Intercept@I.values),
        " views and ",
        dim(object@Intercept@I.values)[3],
        " rounds.\n",
        sep = "")
    cat("\n")
    cat("Intercept beta values:\n")
    for (i in 1:dim(object@Intercept@I.values)[3]) {
      cat(paste(dimnames(object@Intercept@I.values)[[3]][i], ":\n", sep = ""))
      print(object@Intercept@I.values[,,i])
      cat("\n")
    }
    cat("Intercept beta errors:\n")
    for (i in 1:dim(object@Intercept@I.values)[3]) {
      cat(paste(dimnames(object@Intercept@I.errors)[[3]][i], ":\n", sep = ""))
      print(object@Intercept@I.errors[,,i])
      cat("\n")
    }
    cat("Intercept beta z-scores:\n")
    for (i in 1:dim(object@Intercept@I.values)[3]) {
      cat(paste(dimnames(object@Intercept@I.z)[[3]][i], ":\n", sep = ""))
      print(object@Intercept@I.z[,,i])
      cat("\n")
    }
    cat("Intercept beta p-values:\n")
    for (i in 1:dim(object@Intercept@I.values)[3]) {
      cat(paste(dimnames(object@Intercept@I.sig)[[3]][i], ":\n", sep = ""))
      print(object@Intercept@I.sig[,,i])
      cat("\n")
    }
    cat("\n")
    cat("\n")
    cat("An object of class ", class(object@Shape), "\n", sep = "")
    cat("Fits shape coefficients at ",
        length(object@Shape@Shape.set),
        " positions (",
        paste(object@Shape@Shape.set, collapse = ", "),
        ") for ",
        nrow(object@Shape@Shape.values),
        " kinds of shape parameter(s) in a feature model of length ",
        (object@Shape@fpLen),
        ".\n",
        sep = "")
    
    if (nrow(object@Shape@Shape.values) > 0) {
      if (!all(object@Shape@Shape.equivMat == 0)) {
        if (all(object@Shape@Shape.equivMat[,1:(object@Shape@fpLen %/%2)] == 0)) {
          halfMatCheck = object@Shape@Shape.equivMat[1:(object@Shape@fpLen %/%2+object@Shape@fpLen %%2),
                                               (object@Shape@fpLen %/%2+1):object@Shape@fpLen]
          halfMatCheck = halfMatCheck[,ncol(halfMatCheck):1]
          if (all(diag(halfMatCheck) == -1)) {
            cat("Shape features are reverse complement symmetric.\n")
          }
        }
      }
      cat("\n")
      cat("Shape beta values:\n")
      print(getValues(object@Shape))
      cat("\n")
      cat("Shape beta errors:\n")
      print(getErrors(object@Shape))
      cat("\n")
      cat("Shape beta z-scores:\n")
      print(getZ(object@Shape))
      cat("\n")
      cat("Shape beta p-values:\n")
      print(getSig(object@Shape))
      cat("\n")
    }
    invisible(NULL)
  }
)


#' Length of \linkS4class{featureSet}
#'
#' Defines 'length' method for object of class \linkS4class{featureSet}.
#' 
#' The length of \linkS4class{featureSet} object is defined as the footprint length of the model.
#' The footprint length is given by \code{fpLen} = \code{seedLen}+\code{upFootprintExtend}+\code{downFootprintExtend}.
#' 
#' @param x Object of class \linkS4class{featureSet}
#' @export
setMethod(
  f = "length",
  signature = "featureSet",
  definition = function(x) {
    length = x@seedLen+x@upFootprintExtend+x@downFootprintExtend
    return (length)
  }
)


#' Get Seed Length for \linkS4class{featureSet}
#' 
#' Gets seedLen slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getSeedLen",
  signature = "featureSet",
  definition = function(object) {
    return(object@seedLen)
  }
)

#' Get Length Upstream Footprint Extension
#'
#' Gets 'upFootprintExtend' slot from object of class \linkS4class{featureSet}.
#'
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getUpFootprintExtend",
  signature = "featureSet",
  definition = function(object) {
    return(object@upFootprintExtend)
  }
)

#' Get Length Downstream Footprint Extension
#'
#' Gets 'downFootprintExtend' slot from object of class \linkS4class{featureSet}.
#'
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getDownFootprintExtend",
  signature = "featureSet",
  definition = function(object) {
    return(object@downFootprintExtend)
  }
)

#' Get Number of Views for \linkS4class{featureSet}
#' 
#' Gets 'numViews' slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getNumViews",
  signature = "featureSet",
  definition = function(object) {
    return(object@numViews)
  }
)

#' Get Reverse Complement Symmetric (logical)
#' 
#' Gets 'rcSymmetric' slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getRcSymmetric",
  signature = "featureSet",
  definition = function(object) {
    return(object@rcSymmetric)
  }
)

#' Get Rounds
#' 
#' Gets 'rounds' slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getRounds",
  signature = "featureSet",
  definition = function(object) {
    return(object@rounds)
  }
)


#' Get Shape Parameters Used 
#' 
#' Gets 'shapeParamsUsed' slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getShapeParamsUsed",
  signature = "featureSet",
  definition = function(object) {
    return(object@shapeParamsUsed)
  }
)

#' Get \linkS4class{N} for \linkS4class{featureSet}
#' 
#' Gets \linkS4class{N} slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getN",
  signature = "featureSet",
  definition = function(object) {
    return(object@N)
  }
)

#' Get \linkS4class{Intercept} for \linkS4class{featureSet}
#' 
#' Gets \linkS4class{Intercept} slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getIntercept",
  signature = "featureSet",
  definition = function(object) {
    return(object@Intercept)
  }
)

#' Get \linkS4class{Shape} for \linkS4class{featureSet}
#' 
#' Gets \linkS4class{Shape} slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getShape",
  signature = "featureSet",
  definition = function(object) {
    return(object@Shape)
  }
)

#' Get 'GlmFits' for class \linkS4class{featureSet}
#' 
#' Gets 'GlmFits' slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getGlmFits",
  signature = "featureSet",
  definition = function(object) {
    return(object@glmFits)
  }
)

#' Set 'GlmFits' for class \linkS4class{featureSet}
#' 
#' Sets 'GlmFits' slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setReplaceMethod(
  f = "setGlmFits",
  signature = "featureSet",
  definition = function(object, value) {
    stop("GLM fit summaries must be added using addNewBetas function.")
    return(object)
  }
)


#' Get Design Matrix Summary for \linkS4class{featureSet}
#' 
#' Gets 'designMatrixSummary' slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setMethod(
  f = "getDesignMatrixSummary",
  signature = "featureSet",
  definition = function(object) {
    return(object@designMatrixSummary)
  }
)

#' Set Design Matrix Summary for \linkS4class{featureSet}
#' 
#' Sets 'designMatrixSummary' slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
setReplaceMethod(
  f = "setDesignMatrixSummary",
  signature = "featureSet",
  definition = function(object, value) {
    stop("Design Matrix summaries must be added using addNewBetas function.")
    return(object)
  }
)

#' Set Seed Length for \linkS4class{featureSet}
#' 
#' Sets seedLen slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @param value Seed length for 'model' object to which \linkS4class{featureSet} object belongs.
#' @export
setReplaceMethod(
  f = "setSeedLen",
  signature = "featureSet",
  definition = function(object, value) {
    stop('Changes to the seed length cannot be made at the level of "featureSet" class objects.  "seedLen" slot should be changed in the full "model" object. Note that any changes in the seed length must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)


#' Set Length Upstream Footprint Extension
#'
#' Sets 'upFootprintExtend' slot from object of class \linkS4class{featureSet}.
#'
#' @param object Object of class \linkS4class{featureSet}.
#' @param value upFootprintExtend for 'model' object to which \linkS4class{featureSet} object belongs.
#' @export
setReplaceMethod(
  f = "setUpFootprintExtend",
  signature = "featureSet",
  definition = function(object, value) {
    stop("Changes to upFootprintExtend slot cannot be made at the level of 'featureSet' class objects.  'upFootprintExtend' slot should be changed in the full 'model' object. Note that any changes in the 'upFootprintExtend' slot must be made to 'model' class objects before any beta information is added.")
    validObject(object)
    return(object)
  }
)

#' Set Length Downstream Footprint Extension
#'
#' Sets 'downFootprintExtend' slot from object of class \linkS4class{featureSet}.
#'
#' @param object Object of class \linkS4class{featureSet}.
#' @param value downFootprintExtend for 'model' object to which \linkS4class{featureSet} object belongs.
#' @export
setReplaceMethod(
  f = "setDownFootprintExtend",
  signature = "featureSet",
  definition = function(object, value) {
    stop("Changes to downFootprintExtend slot cannot be made at the level of 'featureSet' class objects. 'downFootprintExtend' slot should be changed in the full 'model' object. Note that any changes in the 'downFootprintExtend' slot must be made to 'model' class objects before any beta information is added.")
    validObject(object)
    return(object)
  }
)


#' Set Number of Views for \linkS4class{featureSet}
#' 
#' Sets numViews slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @param value NumViews for \linkS4class{model} object to which \linkS4class{featureSet} object belongs
#' @export
setReplaceMethod(
  f = "setNumViews",
  signature = "featureSet",
  definition = function(object, value) {
    stop('Changes to numViews slot cannot be made at the level of "featureSet" class objects. "NumViews" slot should be changed in the full "model" object.  Note that any changes in the number of views slot must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)

#' Set Reverse Complement Symmetry (logical)
#' 
#' Sets rcSymmetric slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @param value RcSymmetric for \linkS4class{model} object to which \linkS4class{featureSet} object belongs.
#' @export
setReplaceMethod(
  f = "setRcSymmetric",
  signature = "featureSet",
  definition = function(object, value) {
    stop('Changes to rcSymmetric slot cannot be made at the level of "featureSet" class objects. "rcSymmetric" slot should be changed in the full "model" object.  Note that any changes in the rcSymmetric slot must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)

#' Set Rounds for \linkS4class{featureSet}
#' 
#' Sets rounds slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @param value Rounds for \linkS4class{model} object to which \linkS4class{featureSet} object belongs.
#' @export
setReplaceMethod(
  f = "setRounds",
  signature = "featureSet",
  definition = function(object, value) {
    stop('Changes to rounds slot cannot be made at the level of "featureSet" class objects. "rounds" slot should be changed in the full "model" object.  Note that any changes in the rounds slot must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)

#' Set Shape Parameters Used for \linkS4class{featureSet}
#' 
#' Sets shapeParamsUsed slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @param value ShapeParamsUsed for \linkS4class{model} object to which \linkS4class{featureSet} object belongs.
#' @export
setReplaceMethod(
  f = "setShapeParamsUsed",
  signature = "featureSet",
  definition = function(object, value) {
    stop('Changes to shapeParamsUsed slot cannot be made at the level of "featureSet" class objects. "shapeParamsUsed" slot should be changed in the full "model" object.  Note that any changes in the shapeParamsUsed slot must be made to "model" class objects before any beta information is added.')
    validObject(object)
    return(object)
  }
)


#' Set \linkS4class{N} for class \linkS4class{featureSet}
#' 
#' Sets \linkS4class{N} slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @param value Object of class \linkS4class{N}.
#' @export
setReplaceMethod(
  f = "setN",
  signature = "featureSet",
  definition = function(object, value) {
    oldN = getN(object)
    if (oldN@seedLen != value@seedLen) {
      stop('Current and replacement versions of N must have the same seed length.')
    }
    if (oldN@fS.upFootprintExtend != value@fS.upFootprintExtend) {
      stop('Current and replacement versions of N must have the same value of fS.upFootprintExtend.')
    }
    if (oldN@fS.downFootprintExtend != value@fS.downFootprintExtend) {
      stop('Current and replacement versions of N must have the same value of fS.downFootprintExtend.')
    }
    if (!all(c(ncol(getValues(value)),ncol(getErrors(value)), ncol(getZ(value)), ncol(getSig(value)), ncol(getOldValues(value)),
               ncol(getOldErrors(value)), ncol(getOldZ(value)), ncol(getOldSig(value))) == ncol(getValues(oldN)))) {
      stop("All beta related slots in new 'N' must have a number of columns equal to the footprint length specified by seedLen, upFootprintExtend, and downFootprintExtend.")
    }
    
    
    setValues(value) = getValues(value)
    setErrors(value) = getErrors(value)
    setZ(value) = getZ(value)
    setSig(value) = getSig(value)
    setOldValues(value) = getOldValues(value)
    setOldErrors(value) = getOldErrors(value)
    setOldZ(value) = getOldZ(value)
    setOldSig(value) = getOldSig(value)
    setEquivMat(value) = getEquivMat(value)
    validObject(value)
    object@N <- value
    validObject(object)
    return(object)
  }
)


#' Set \linkS4class{Intercept} for class \linkS4class{featureSet}
#' 
#' Sets \linkS4class{Intercept} slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @param value Object of class \linkS4class{Intercept}.
#' @export
setReplaceMethod(
  f = "setIntercept",
  signature = "featureSet",
  definition = function(object, value) {
    oldI = getIntercept(object)
    if (!all(c(ncol(getValues(value)),ncol(getErrors(value)), ncol(getZ(value)), ncol(getSig(value)), ncol(getOldValues(value)),
               ncol(getOldErrors(value)), ncol(getOldZ(value)), ncol(getOldSig(value))) == ncol(getValues(oldI)))) {
      stop("All beta related slots in new 'Intercept' must have a number of columns equal to the footprint length specified by seedLen, upFootprintExtend, and downFootprintExtend.")
    }
    if (!all(c(dim(getValues(value))[3],dim(getErrors(value))[3], dim(getZ(value))[3], dim(getSig(value))[3], dim(getOldValues(value))[3],
               dim(getOldErrors(value))[3], dim(getOldZ(value))[3], dim(getOldSig(value))[3]) == dim(getValues(oldI))[3])) {
      stop("All beta related slots in new 'Intercept' must have the same number of rounds (dim[3]).")
    }
    
    
    setValues(value) = getValues(value)
    setErrors(value) = getErrors(value)
    setZ(value) = getZ(value)
    setSig(value) = getSig(value)
    setOldValues(value) = getOldValues(value)
    setOldErrors(value) = getOldErrors(value)
    setOldZ(value) = getOldZ(value)
    setOldSig(value) = getOldSig(value)
    validObject(value)
    object@Intercept <- value
    validObject(object)
    return(object)
  }
)


#' Set \linkS4class{Shape} for class \linkS4class{featureSet}
#' 
#' Sets \linkS4class{Shape} slot from object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @param value Object of class \linkS4class{Shape}.
#' @export
setReplaceMethod(
  f = "setShape",
  signature = "featureSet",
  definition = function(object, value) {
    oldS = getShape(object)
    if (oldS@seedLen != value@seedLen) {
      stop('Current and replacement versions of Shape slot must have the same seed length.')
    }
    if (oldS@fS.upFootprintExtend != value@fS.upFootprintExtend) {
      stop('Current and replacement versions of Shape slot must have the same value of fS.upFootprintExtend.')
    }
    if (oldS@fS.downFootprintExtend != value@fS.downFootprintExtend) {
      stop('Current and replacement versions of Shape slot must have the same value of fS.downFootprintExtend.')
    }
    if (!all(c(ncol(getValues(value)),ncol(getErrors(value)), ncol(getZ(value)), ncol(getSig(value)), ncol(getOldValues(value)),
               ncol(getOldErrors(value)), ncol(getOldZ(value)), ncol(getOldSig(value))) == ncol(getValues(oldS)))) {
      stop("All beta related slots in new 'Shape' must have a number of columns equal to the footprint length specified by seedLen, upFootprintExtend, and downFootprintExtend.")
    }
    if (!all(c(nrow(getValues(value)),nrow(getErrors(value)), nrow(getZ(value)), nrow(getSig(value)), nrow(getOldValues(value)),
               nrow(getOldErrors(value)), nrow(getOldZ(value)), nrow(getOldSig(value))) == nrow(getValues(oldS)))) {
      stop("All beta related slots in new 'Shape' must have a number of rows equal to the footprint length specified by seedLen, upFootprintExtend, and downFootprintExtend.")
    }
    
    
    setValues(value) = getValues(value)
    setErrors(value) = getErrors(value)
    setZ(value) = getZ(value)
    setSig(value) = getSig(value)
    setOldValues(value) = getOldValues(value)
    setOldErrors(value) = getOldErrors(value)
    setOldZ(value) = getOldZ(value)
    setOldSig(value) = getOldSig(value)
    setEquivMat(value) = getEquivMat(value)
    validObject(value)
    object@Shape <- value
    validObject(object)
    return(object)
  }
)


















