
#' Validity Function for \linkS4class{model}
#' 
#' Validity Function for object of class \linkS4class{model}.
#' 
#' @param object Object of class \linkS4class{model}.
#' @export
validModel = function(object) {
  if (object@varRegLen <= 0) {
    return("'varRegLen' must be a positive integer.")
  }
  if (length(grep(paste("[ACGT]{", nchar(object@leftFixedSeq), "}", sep = ""), object@leftFixedSeq)) != 1) {
    return("'leftFixedSeq' must be a character object composed of A,C,G, and T only.")
  }
  if (length(grep(paste("[ACGT]{", nchar(object@rightFixedSeq), "}", sep = ""), object@rightFixedSeq)) != 1) {
    return("'rightFixedSeq' must be a character object composed of A,C,G, and T only.")
  }
  if (object@seedLen <= 0) {
    return("'seedLen' must be a positive integer.")
  }
  if (object@seedLen != countDNAseqLength(object@consensusSeq)) {
    return("'seedLen' must be the same length as object@consensusSeq.")
  }

  if (object@leftFixedSeqOverlap < 0) {
    return("'leftFixedSeqOverlap' must be a non-negative integer.")
  }
  if (object@leftFixedSeqOverlap > nchar(object@leftFixedSeq)) {
    return("'leftFixedSeqOverlap' must be < nchar('leftFixedSeq').")
  }
  
  if (object@rightFixedSeqOverlap < 0) {
    return("'rightFixedSeqOverlap' must be a non-negative integer.")
  }
  if (object@rightFixedSeqOverlap > nchar(object@rightFixedSeq)) {
    return("'rightFixedSeqOverlap' must be < nchar('rightFixedSeq').")
  }
  
  if (object@numViews != object@leftFixedSeqOverlap+object@rightFixedSeqOverlap+object@varRegLen-object@fpLen+1) {
    return("'numViews' must equal 'lefFixedSeqOverlap'+'rightFixedSeqOverlap'+'varRegLen'-'fpLen'+1 ")
  }
  if (length(object@rounds) != 1) {
    return("'round' must be a list of length one. 'round[[1]]' should be a vector with lenght >= 1 containing the rounds of data to be analyzed.")
  }
  if ((object@confidenceLevel <= 0) | (object@confidenceLevel >= 1)) {
    return("'confidenceLevel' must be in the open interval from 0 to 1.")
  }
  if ((object@minAffinity < 0) | (object@minAffinity > 1)) {
    return("'minAffinity' must be in the closed interval from 0 to 1.")
  }
  if (object@missingValueSuppression < 0) {
    return("'missingValueSuppression' must be >= 0.")
  }
  if ((object@minSeedValue < 0) | (object@minSeedValue >= 1)) {
    return("'minSeedValue' must be in the half-open interval [0, 1).")
  }
  if (object@upFootprintExtend < 0) {
    return("'upFootprintExtend' must be >= 0.")
  }
  
  if (all(c(object@features@N@N.upFootprintExtend, object@features@Shape@Shape.upFootprintExtend) < object@upFootprintExtend)) {
    if ((object@useFixedValuesOffset.N == FALSE) & (object@useFixedValuesOffset.Shape)) {
      return("If 'useFixedValuesOffset.N' == 'useFixedValuesOffset.Shape' == FALSE, either Shape.upFootprintExtend or N.upFootprintExtend must be equal to upFootprintExtend.")
    }
  }
  
  if (all(c(object@features@N@N.downFootprintExtend, object@features@Shape@Shape.downFootprintExtend) < object@downFootprintExtend)) {
    if ((object@useFixedValuesOffset.N == FALSE) & (object@useFixedValuesOffset.Shape)) {
      return("If 'useFixedValuesOffset.N' == 'useFixedValuesOffset.Shape' == FALSE, either Shape.downFootprintExtend or N.downFootprintExtend must be equal to downFootprintExtend.")
    }
  }
  
  if (object@downFootprintExtend < 0) {
    return("'downFootprintExtend' must be >= 0.")
  }
  if ((object@rcSymmetric == TRUE) & (object@includeDNAstrand == TRUE)) {
    return("'rcSymmetric' and 'includeDNAstrand' cannot both be true for the same 'model' object.")
  }
  if (object@exUpBases < 0) {
    return("'exUpBases' must be >= 0.")
  }
  if (object@exDownBases < 0) {
    return("'exDownBases' must be >= 0.")
  }
  if (object@fpLen > object@varRegLen+object@leftFixedSeqOverlap+object@rightFixedSeqOverlap) {
    return("'fpLen' must be shorter than the DNA strands specified by 'varRegLen', 'leftFixedSeq', and 'rightFixedSeq'.")
  } 
  if (object@fpLen != object@seedLen+object@upFootprintExtend+object@downFootprintExtend) {
    return("'fpLen' must equal 'upFootprintExtend'+'downFootprintExtend'+'seedLen'.")
  }
  if (length(object@shapeParamsUsed[[1]]) > 0) {
    if (length(object@shapeTable) == 0) {
      return("'shapeTable' must be supplied if shapeParamsUsed[[1]] != c(character(0)).")
    }
    for (sP in object@shapeParamsUsed[[1]]) {
      if (sP == "MGW") {
        if (!("MGW" %in% colnames(object@shapeTable))) {
          return ("'MGW' is in shapeParamsUsed[[1]] but is not in the shapeTable supplied.")
        }
      } else if (sP == "ProT") {
        if (!("ProT" %in% colnames(object@shapeTable))) {
          return ("'ProT' is in shapeParamsUsed[[1]] but is not in the shapeTable supplied.")
        }
      } else if (sP == "HelT") {
        if (!("HelTA" %in% colnames(object@shapeTable))) {
          return ("'HelT' is in shapeParamsUsed[[1]] but 'HelTA' is not in the shapeTable supplied.")
        }
        if (!("HelTB" %in% colnames(object@shapeTable))) {
          return ("'HelT' is in shapeParamsUsed[[1]] but 'HelTB' is not in the shapeTable supplied.")
        }
      } else if (sP == "Roll") {
        if (!("RollA" %in% colnames(object@shapeTable))) {
          return ("'Roll' is in shapeParamsUsed[[1]] but 'RollA' is not in the shapeTable supplied.")
        }
        if (!("RollB" %in% colnames(object@shapeTable))) {
          return ("'Roll' is in shapeParamsUsed[[1]] but 'RollB' is not in the shapeTable supplied.")
        }
      } else {
        return (paste("Shape parameter ", sP, " is not in c('HelT', 'MGW', 'ProT', 'Roll').", sep = ""))
      }
    }
  }
  
  if (length(object@shapeTable) > 0) {
    if (!all(colnames(object@shapeTable) %in% c("Sequence", "MGW", "ProT", "HelTA", "HelTB", "RollA", "RollB"))){
      return("'shapeTable' columns must be named Sequence, MGW, ProT, HelTA, HelTB, RollA, or RollB.")
    }
  }
  if (object@features@seedLen != object@seedLen) {
    return("'seedLen' slots in 'features' and 'model' must have the same value.")
  }
  if (object@features@upFootprintExtend != object@upFootprintExtend) {
    return("'upFootprintExtend' slots in 'features' and 'model' must have the same value.")
  }
  if (object@features@downFootprintExtend != object@downFootprintExtend) {
    return("'downFootprintExtend' slots in 'features' and 'model' must have the same value.")
  }
  if (object@features@rcSymmetric != object@rcSymmetric) {
    return("'rcSymmetric' slots in 'features' and 'model' must have the same value.")
  }
  if (object@features@numViews != object@numViews) {
    return("'numViews' slots in 'features' and 'model' must have the same value.")
  }
  if (length(object@rounds[[1]]) == length(object@features@rounds[[1]])) {
    if (!all(object@rounds[[1]] == object@features@rounds[[1]])) {
      return("'rounds' slots in 'features' and 'model' must have the same value.")
    }
  } else {
    return("'rounds' slots in 'features' and 'model' must have the same value.")
  }
  
  # Check shapeParamsUsed
  if (!is.list(object@shapeParamsUsed)) {
    return("'model@shapeParamsUsed' must be a list object of length 1.")
  }
  if (length(object@shapeParamsUsed) != 1) {
    return("'model@shapeParamsUsed' must be a list object of length 1.")
  }
  if (length(object@shapeParamsUsed[[1]]) > 0) {
    shapeParams = object@shapeParamsUsed[[1]]
    shapeParams = shapeParams[order(shapeParams)]
    if (!all(shapeParams %in% c("HelT", "MGW", "ProT", "Roll"))) {
      return(paste("shapeParamsUsed[[1]] = c(", paste(shapeParams, collapse = ", "), ") is not valid. All shape parameters must be in c('HelT', 'MGW', 'ProT', 'Roll').", sep = ""))
    }
    if (length(object@shapeParamsUsed[[1]]) != length(object@features@shapeParamsUsed[[1]])) {
      return(paste("model@shapeParamsUsed[[1]] (",
                   paste(object@shapeParamsUsed[[1]],
                         collapse = ", "),
                   ") != model@features@shapeParamsUsed[[1]](",
                   paste(object@features@shapeParamsUsed[[1]], collapse = ", "),
                   ").",
                   sep = ""))
    } else if (!all(object@shapeParamsUsed[[1]] == object@features@shapeParamsUsed[[1]])) {
      return(paste("model@ShapeParamsUsed[[1]] (",
                   paste(object@shapeParamsUsed[[1]],
                         collapse = ", "),
                   ") != model@features@ShapeParamsUsed[[1]] (",
                   paste(object@features@shapeParamsUsed[[1]],
                         collapse = ", "),
                   ").",
                   sep = ""))
    }
  } else if (length(object@features@shapeParamsUsed[[1]]) != 0) {
    return(paste("model@ShapeParamsUsed[[1]] (",
                 paste(object@shapeParamsUsed[[1]],
                       collapse = ", "),
                 ") != model@features@ShapeParamsUsed[[1]] (",
                 paste(object@features@shapeParamsUsed[[1]],
                       collapse = ", "),
                 ").",
                 sep = ""))
  }
  
  if ((object@useFixedValuesOffset.Shape == TRUE) & (object@includeShape == FALSE)) {
    return("'useFixedValuesOffset.Shape' must be used with 'includeShape' = TRUE.")
  }
  
  validObject(object@features)
  return(TRUE)
}

