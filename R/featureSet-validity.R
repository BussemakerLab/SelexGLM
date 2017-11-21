
#' Validity Function for \linkS4class{featureSet}
#' 
#' Validity Function for object of class \linkS4class{featureSet}.
#' 
#' @param object Object of class \linkS4class{featureSet}.
#' @export
validFeatureSet = function(object) {
  validObject(object@Shape)
  validObject(object@Intercept)
  validObject(object@N)
  
  # Check seedLen parameter
  if (object@seedLen <= 0) {
    return ('seedLen slot must be a positive value.')
  }
  if (length(unique(c(object@N@seedLen, object@Shape@seedLen, object@seedLen))) != 1) {
    return("seedLen slots in 'N', 'Shape', and 'featureSet' must all be the same.")
  }
  
  # Check upFootprintExtend parameter
  if (length(unique(c(object@N@fS.upFootprintExtend, object@Shape@fS.upFootprintExtend, object@upFootprintExtend))) != 1) {
    return("Must have featureSet@upFootprintExtend == featureSet@N@fS.upFootprintExtend == featureSet@Shape@fS.upFootprintExtend.")
  }
  
  # Check downFootprintExtend parameter
  if (length(unique(c(object@N@fS.downFootprintExtend, object@Shape@fS.downFootprintExtend, object@downFootprintExtend))) != 1) {
    return("Must have featureSet@downFootprintExtend == featureSet@N@fS.downFootprintExtend == featureSet@Shape@fS.downFootprintExtend.")
  }
  
  
  # Check numViews
  if (object@numViews <= 0) {
    return ('featureSet@numViews slot must be a positive value.')
  }
  
  if (object@numViews != ncol(object@Intercept@I.values)) {
    return("Must have featureSet@numViews == ncol(featureSet@Intercept@I.values).")
  }
  
  # Check rounds
  if (!is.list(object@rounds)) {
    return("'featureSet@rounds' must be a list object of length 1.")
  }
  if (length(object@rounds) != 1) {
    return("'featureSet@rounds' must be a list object of length 1.")
  }
  if (length(object@rounds[[1]]) < 1) {
    return("'featureSet@rounds[[1]]' must contain at least one round variable (length(featureSet@rounds[[1]]) >= 1).")
  }
  Irounds = as.numeric(gsub("Round.", "", dimnames(object@Intercept@I.values)[[3]]))
  Irounds = Irounds[order(Irounds)]
  if (length(object@rounds[[1]]) != length(Irounds)) {
    return("featureSet@rounds must be equivalent to dimnames(featureSet@Intercept@I.values)[[3]].")
  } else if (!all(object@rounds[[1]][order(object@rounds[[1]])] == Irounds[order(Irounds)])) {
    return("featureSet@rounds must be equivalent to dimnames(featureSet@Intercept@I.values)[[3]].")
  }
  
  # Check shapeParamsUsed
  if (!is.list(object@shapeParamsUsed)) {
    return("'featureSet@shapeParamsUsed' must be a list object of length 1.")
  }
  if (length(object@shapeParamsUsed) != 1) {
    return("'featureSet@shapeParamsUsed' must be a list object of length 1.")
  }
  if (length(object@shapeParamsUsed[[1]]) > 0) {
    shapeParams = object@shapeParamsUsed[[1]]
    shapeParams = shapeParams[order(shapeParams)]
    if (!all(shapeParams %in% c("HelT", "MGW", "ProT", "Roll"))) {
      return(paste("shapeParamsUsed[[1]] = c(", paste(shapeParams, collapse = ", "), ") is not valid. All shape parameters must be in c('HelT', 'MGW', 'ProT', 'Roll').", sep = ""))
    }
    if (length(object@shapeParamsUsed[[1]]) != length(object@Shape@shapeParamsUsed[[1]])) {
      return(paste("featureSet@shapeParamsUsed[[1]] (",
                   paste(object@shapeParamsUsed[[1]],
                         collapse = ", "),
                   ") != featureSet@Shape@shapeParamsUsed[[1]](",
                   paste(object@Shape@shapeParamsUsed[[1]], collapse = ", "),
                   ").",
                   sep = ""))
    } else if (!all(object@shapeParamsUsed[[1]] == object@Shape@shapeParamsUsed[[1]])) {
      return(paste("featureSet@ShapeParamsUsed[[1]] (",
                   paste(object@shapeParamsUsed[[1]],
                         collapse = ", "),
                   ") != featureSet@Shape@ShapeParamsUsed[[1]] (",
                   paste(object@Shape@shapeParamsUsed[[1]],
                         collapse = ", "),
                   ").",
                   sep = ""))
    }
  } else if (length(object@Shape@shapeParamsUsed[[1]]) != 0) {
    return(paste("featureSet@ShapeParamsUsed[[1]] (",
                 paste(object@shapeParamsUsed[[1]],
                       collapse = ", "),
                 ") != featureSet@Shape@ShapeParamsUsed[[1]] (",
                 paste(object@Shape@shapeParamsUsed[[1]],
                       collapse = ", "),
                 ").",
                 sep = ""))
  }

  # Check rcSymmetric
  if (object@rcSymmetric == TRUE) {
    NmatTest = buildSymmetricEquivalenceMatrix(object@N)
    SmatTest = buildSymmetricEquivalenceMatrix(object@Shape)
    if (!all(object@N@N.equivMat == NmatTest)) {
      return("featureSet@N@N.equivMat is not consistent with a reverse complement symmetric featureSet.")
    }
    if (!all(object@Shape@Shape.equivMat == SmatTest)) {
      return("featureSet@Shape@Shape.equivMat is not consistent with a reverse complement symmetric featureSet.")
    }
    
    # TESTING Nucleotide specific affinity
    if (!testNBetaSymmetry(object@N@N.values, object@N@N.set)) {
      return("featureSet@N@N.values is not consistent with a reverse complement symmetric featureSet.")
    }
    if (!testNBetaSymmetry(object@N@N.errors, object@N@N.set)) {
      return("featureSet@N@N.errors is not consistent with a reverse complement symmetric featureSet.")
    }
    
    if (!testNBetaSymmetry(object@N@N.z, object@N@N.set)) {
      return("featureSet@N@N.z is not consistent with a reverse complement symmetric featureSet.")
    }
    
    if (!testNBetaSymmetry(object@N@N.sig, object@N@N.set)) {
      return("featureSet@N@N.sig is not consistent with a reverse complement symmetric featureSet.")
    }
    
    # TESTING Shape specific affinity
    if (!testShapeBetaSymmetry(object@Shape@Shape.values, object@Shape@Shape.set)) {
      return("featureSet@Shape@Shape.values is not consistent with a reverse complement symmetric featureSet.")
    }
    if (!testShapeBetaSymmetry(object@Shape@Shape.errors, object@Shape@Shape.set)) {
      return("featureSet@Shape@Shape.errors is not consistent with a reverse complement symmetric featureSet.")
    }
    
    if (!testShapeBetaSymmetry(object@Shape@Shape.z, object@Shape@Shape.set)) {
      return("featureSet@Shape@Shape.z is not consistent with a reverse complement symmetric featureSet.")
    }
    
    if (!testShapeBetaSymmetry(object@Shape@Shape.sig, object@Shape@Shape.set)) {
      return("featureSet@Shape@Shape.sig is not consistent with a reverse complement symmetric featureSet.")
    }
    
    # TESTING View specific affinity
    if (!testViewBetaSymmetry(object@Intercept@I.values)) {
      return("featureSet@Intercept@I.values is not consistent with a reverse complement symmetric featureSet.")
    }
    if (!testViewBetaSymmetry(object@Intercept@I.errors)) {
      return("featureSet@Intercept@I.errors is not consistent with a reverse complement symmetric featureSet.")
    }
    
    if (!testViewBetaSymmetry(object@Intercept@I.z)) {
      return("featureSet@Intercept@I.z is not consistent with a reverse complement symmetric featureSet.")
    }
    
    if (!testViewBetaSymmetry(object@Intercept@I.sig)) {
      return("featureSet@Intercept@I.sig is not consistent with a reverse complement symmetric featureSet.")
    }
    
  }
  
  return (TRUE)
}

#' Test Symmetry of Nucleotide Betas
#' 
#' Tests reverse complement symmetry of nucleotide betas for nucleotides in N.set
#' 
#' @param mat Matrix of nucleotide beta values from \linkS4class{N}.
#' @param Nset Set of nucleotide values included in model (i.e. N.set from \linkS4class{N}).
testNBetaSymmetry = function(mat, NSet) {
  fpLen = ncol(mat)
  for (i in 1:fpLen) {
    rc.i = fpLen-i+1
    if (!(i %in% NSet)) {
      mat[,i] = 0
      mat[,rc.i] = 0
    }
  }
  Nm.h1 = mat[1:4, 1:ceiling(fpLen/2)]
  Nm.h2 = mat[4:1, fpLen:(floor(fpLen/2)+1)]
  epsilon <- 1e-10
  if (!all(abs(Nm.h1- Nm.h2) < epsilon, na.rm = TRUE)){
    return (FALSE)
  } else {
    return (TRUE)
  }
}

#' Test Symmetry of Shape Betas
#' 
#' Tests reverse complement symmetry of shape betas for nucleotides in Shape.set
#' 
#' @param mat Matrix of shape beta values from \linkS4class{Shape}.
#' @param Sset Set of shape values included in model (i.e. Shape.set from \linkS4class{Shape}).
testShapeBetaSymmetry = function(mat, SSet) {
  fpLen = ncol(mat)
  for (i in 1:fpLen) {
    rc.i = fpLen-i+1
    if (!(i %in% SSet)) {
      mat[,i] = 0
    }
  }
  Sm.h1 = mat[, 1:ceiling(fpLen/2)]
  Sm.h2 = mat[, fpLen:(floor(fpLen/2)+1)]
  HelTtest = grep("Shape.HelTA", rownames(mat))
  if (length(HelTtest) == 1) {
    Sm.h2[c(HelTtest,HelTtest+1), 1:(fpLen %/% 2)] = Sm.h2[c(HelTtest+1,HelTtest), 1:(fpLen %/% 2)]
  }
  Rolltest = grep("Shape.RollA", rownames(mat))
  if (length(Rolltest) == 1) {
    Sm.h2[c(Rolltest,Rolltest+1), 1:(fpLen %/% 2)] = Sm.h2[c(Rolltest+1,Rolltest), 1:(fpLen %/% 2)]
  }
  epsilon <- 1e-10
  if (!all(abs(Sm.h1- Sm.h2) < epsilon, na.rm = TRUE)){
    return (FALSE)
  } else {
    return (TRUE)
  }
}

#' Test Symmetry of View Betas
#' 
#' Tests reverse complement symmetry of view betas.
#' 
#' @param arr Array of beta values from \linkS4class{Intercept}.
testViewBetaSymmetry = function(arr) {
  if (all(arr["Strand.R",,]  == 0, na.rm = TRUE)) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}

