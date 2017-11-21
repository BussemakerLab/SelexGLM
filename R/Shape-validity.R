
#' Validity Function for \linkS4class{Shape}
#' 
#' Validity Function for object of class \linkS4class{Shape}.
#' 
#' @param object Object of class \linkS4class{Shape}.
#' @export
validShape = function(object) {
  if (object@seedLen <= 0) {
    return('seedLen must be positive.')
  }
  if ((object@Shape.upFootprintExtend < 0) | (object@Shape.upFootprintExtend > object@fS.upFootprintExtend)) {
    return('Must choose Shape.upFootprintExtend such that 0 <= Shape.upFootprintExtend <= fS.upFootprintExtend.')
  }
  if (0 %in% object@Shape.set) {
    if (length(object@Shape.set) != 1) {
      return('If 0 is in Shape.set, it must be the only element of Shape.set.')
    }
    if (object@Shape.downFootprintExtend != 0) {
      return('Shape.downFootprintExtend must be equal to 0 for Shape.set = c(0).')
    }
    if (object@Shape.upFootprintExtend != 0) {
      return('Shape.upFootprintExtend must be equal to 0 for Shape.set = c(0).')
    }
  }

  if ((object@Shape.downFootprintExtend < 0) | (object@Shape.downFootprintExtend > object@fS.downFootprintExtend)) {
    return('Must choose Shape.downFootprintExtend such that 0 <= Shape.downFootprintExtend <= fS.downFootprintExtend.')
  }


  if (min(object@Shape.set) <= object@fS.upFootprintExtend) {
    if ((length(object@Shape.set) == 1) & (object@Shape.set[1] == 0)) {
      if (object@Shape.upFootprintExtend != 0)  {
        return('Shape.upFootprintExtend and Shape.set are inconsistent.')
      }
    } else {
      testMinSet = object@fS.upFootprintExtend-object@Shape.upFootprintExtend+1
      if (testMinSet != min(object@Shape.set)) {
        return('Shape.upFootprintExtend and Shape.set are inconsistent.')
      }
    }

  }

  if (max(object@Shape.set) > object@fS.upFootprintExtend+object@seedLen) {
    if ((length(object@Shape.set) == 1) & (object@Shape.set[1] == 0)) {
      if (object@Shape.downFootprintExtend != 0)  {
        return('Shape.downFootprintExtend and Shape.set are inconsistent.')
      }
    } else {
      testMaxSet = object@fS.upFootprintExtend+object@seedLen+object@Shape.downFootprintExtend
      if (testMaxSet != max(object@Shape.set)) {
        return('Shape.downFootprintExtend and Shape.set are inconsistent.')
      }
    }
  }
  if (object@fS.upFootprintExtend+object@seedLen+object@Shape.downFootprintExtend < max(object@Shape.set)) {
    return('Shape.downFootprintExtend and Shape.set are inconsistent.')
  }
  if (object@fpLen != object@fS.upFootprintExtend+object@fS.downFootprintExtend+object@seedLen) {
    return('fpLen and (seedLen, fS.upFootprintExtend, fS.downFootprintExtend) are inconsistent.')
  }
  if (max(object@Shape.set) > object@fpLen) {
    return('fpLen and Shape.set are inconsistent.')
  }

  if ((nrow(object@Shape.equivMat) != object@fpLen) | (ncol(object@Shape.equivMat) != object@fpLen)) {
    return('fpLen and Shape.equivMat dimensionality are inconsistent.')
  }
  ##############
  ##################
  if (!all((object@Shape.equivMat >= -1) | (object@Shape.equivMat <= 1))) {
    return('Shape.equivMat entries do not fall within allowed range from -1 to 1.')
  }
  
  if (length(object@shapeParamsUsed[[1]]) > 0) {
    shapeParams = object@shapeParamsUsed[[1]]
    shapeParams = shapeParams[order(shapeParams)]
    if (!all(shapeParams %in% c("HelT", "MGW", "ProT", "Roll"))) {
      return(paste("shapeParamsUsed[[1]] = c(", paste(shapeParams, collapse = ", "), ") is not valid. All shape parameters must be in c('HelT', 'MGW', 'ProT', 'Roll').", sep = ""))
    }
    rn = c(character(0))
    if ("HelT" %in% shapeParams) {
      rn = c(rn, "Shape.HelTA", "Shape.HelTB")
    }
    if ("MGW" %in% shapeParams) {
      rn = c(rn, "Shape.MGW")
    }
    if ("ProT" %in% shapeParams) {
      rn = c(rn, "Shape.ProT")
    }
    if ("Roll" %in% shapeParams) {
      rn = c(rn, "Shape.RollA", "Shape.RollB")
    }
    
    if (!all(rownames(object@Shape.values) == rn)) {
      return(paste("rownames('Shape.values') = c(", paste(rownames(object@Shape.values), collapse = ", "), ") is not consistent with shapeParamsUsed: ", paste(shapeParams, collapse = ", "), sep = ""))
    }
    if (!all(rownames(object@Shape.errors) == rn)) {
      return(paste("rownames('Shape.errors') = c(", paste(rownames(object@Shape.errors), collapse = ", "), ") is not consistent with shapeParamsUsed: ", paste(shapeParams, collapse = ", "), sep = ""))
    }
    if (!all(rownames(object@Shape.z) == rn)) {
      return(paste("rownames('Shape.z') = c(", paste(rownames(object@Shape.z), collapse = ", "), ") is not consistent with shapeParamsUsed: ", paste(shapeParams, collapse = ", "), sep = ""))
    }
    if (!all(rownames(object@Shape.sig) == rn)) {
      return(paste("rownames('Shape.sig') = c(", paste(rownames(object@Shape.sig), collapse = ", "), ") is not consistent with shapeParamsUsed: ", paste(shapeParams, collapse = ", "), sep = ""))
    }
    if (!all(rownames(object@Shape.oldValues) == rn)) {
      return(paste("rownames('Shape.oldValues') = c(", paste(rownames(object@Shape.oldValues), collapse = ", "), ") is not consistent with shapeParamsUsed: ", paste(shapeParams, collapse = ", "), sep = ""))
    }
    if (!all(rownames(object@Shape.oldErrors) == rn)) {
      return(paste("rownames('Shape.oldErrors') = c(", paste(rownames(object@Shape.oldErrors), collapse = ", "), ") is not consistent with shapeParamsUsed: ", paste(shapeParams, collapse = ", "), sep = ""))
    }
    if (!all(rownames(object@Shape.oldZ) == rn)) {
      return(paste("rownames('Shape.oldZ') = c(", paste(rownames(object@Shape.oldZ), collapse = ", "), ") is not consistent with shapeParamsUsed: ", paste(shapeParams, collapse = ", "), sep = ""))
    }
    if (!all(rownames(object@Shape.oldSig) == rn)) {
      return(paste("rownames('Shape.oldSig') = c(", paste(rownames(object@Shape.oldSig), collapse = ", "), ") is not consistent with shapeParamsUsed: ", paste(shapeParams, collapse = ", "), sep = ""))
    }
  } else {
    if (!(all(c(nrow(object@Shape.values),
                nrow(object@Shape.errors),
                nrow(object@Shape.z),
                nrow(object@Shape.sig),
                nrow(object@Shape.oldValues),
                nrow(object@Shape.oldErrors),
                nrow(object@Shape.oldZ),
                nrow(object@Shape.oldSig)) == 0))) {
      return('Dimensions of of one of the following slots is inconsistent with shapeParamsUsed[[1]] = c(character(0)): Shape.values, Shape.errors, Shape.z, Shape.sig, Shape.oldValues, Shape.oldErrors, Shape.oldZ, Shape.oldSig.')
    }
    if ((length(object@Shape.set) != 1) | (object@Shape.set[1] != 0)) {
      return(paste("Shape.set = c(", paste(object@Shape.set, collapse = ", "), ") is not consistent with shapeParamsUsed[[1]] = c(character(0)). Shape.set must be c(numeric(0)).", sep = ""))
    }
    if (object@Shape.downFootprintExtend != 0) {
      return(paste("Shape.downFootprintExtend = ",object@Shape.downFootprintExtend, " is not consistent with shapeParamsUsed[[1]] = c(character(0)). Shape.downFootprintExted must equal 0.", sep = ""))
    }
    if (object@Shape.upFootprintExtend != 0) {
      return(paste("Shape.upFootprintExtend = ",object@Shape.upFootprintExtend, " is not consistent with shapeParamsUsed[[1]] = c(character(0)). Shape.upFootprintExtend must equal 0.", sep = ""))
    }
  }
  
  if (!(all(c(ncol(object@Shape.values),
              ncol(object@Shape.errors),
              ncol(object@Shape.z),
              ncol(object@Shape.sig),
              ncol(object@Shape.oldValues),
              ncol(object@Shape.oldErrors),
              ncol(object@Shape.oldZ),
              ncol(object@Shape.oldSig)) == object@fpLen))) {
    return('Dimensions of of one of the following slots is inconsistent with fpLen: Shape.values, Shape.errors, Shape.z, Shape.sig, Shape.oldValues, Shape.oldErrors, Shape.oldZ, Shape.oldSig.')
  }
  
  if (length(unique(c(dim(object@Shape.oldValues)[3],
                      dim(object@Shape.oldErrors)[3],
                      dim(object@Shape.oldZ)[3],
                      dim(object@Shape.oldSig)[3]))) != 1) {
    return('Dimensions of of one of the following slots is inconsistent with the others: Shape.oldValues, Shape.oldErrors, Shape.oldZ, Shape.oldSig.')
  }

  return(TRUE)
  
}

