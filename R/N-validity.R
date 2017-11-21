
#' Validity Function for \linkS4class{N}
#' 
#' Validity Function for object of class \linkS4class{N}.
#' 
#' @param object Object of class \linkS4class{N}.
#' @export 
N.validity = function(object) {
  if (object@seedLen <= 0) {
    return('seedLen must be positive.')
  }
  if ((object@N.upFootprintExtend < 0) | (object@N.upFootprintExtend > object@fS.upFootprintExtend)) {
    return('Must choose N.upFootprintExtend such that 0 <= N.upFootprintExtend <= fS.upFootprintExtend.')
  }
  if (0 %in% object@N.set) {
    if (length(object@N.set) != 1) {
      return('If 0 is in N.set, it must be the only element of N.set.')
    }
    if (object@N.downFootprintExtend != 0) {
      return('N.downFootprintExtend must be equal to 0 for N.set = c(0).')
    }
    if (object@N.upFootprintExtend != 0) {
      return('N.upFootprintExtend must be equal to 0 for N.set = c(0).')
    }
  }
  if (min(object@N.set) <= object@fS.upFootprintExtend) {
    if ((length(object@N.set) == 1) & (object@N.set[1] == 0)) {
      if (object@N.upFootprintExtend != 0)  {
        return('N.upFootprintExtend and N.set are inconsistent.')
      }
    } else {
      testMinSet = object@fS.upFootprintExtend-object@N.upFootprintExtend+1
      if (testMinSet != min(object@N.set)) {
        return('N.upFootprintExtend and N.set are inconsistent.')
      }
    }

  }


  if ((object@N.downFootprintExtend < 0) | (object@N.downFootprintExtend > object@fS.downFootprintExtend)) {
    return('Must choose N.downFootprintExtend such that 0 <= N.downFootprintExtend <= fS.downFootprintExtend.')
  }
  if (max(object@N.set) > object@fS.upFootprintExtend+object@seedLen) {
    if ((length(object@N.set) == 1) & (object@N.set[1] == 0)) {
      if (object@N.downFootprintExtend != 0)  {
        return('N.downFootprintExtend and N.set are inconsistent.')
      }
    } else {
      testMaxSet = object@fS.upFootprintExtend+object@seedLen+object@N.downFootprintExtend
      if (testMaxSet != max(object@N.set)) {
        return('N.downFootprintExtend and N.set are inconsistent.')
      }
    }
  }

  if (object@fpLen != object@fS.upFootprintExtend+object@fS.downFootprintExtend+object@seedLen) {
    return('fpLen and (seedLen, fS.upFootprintExtend, fS.downFootprintExtend) are inconsistent.')
  }
  if (max(object@N.set) > object@fpLen) {
    return('fpLen and N.set are inconsistent.')
  }
  if ((nrow(object@N.equivMat) != object@fpLen) | (ncol(object@N.equivMat) != object@fpLen)) {
    return('fpLen and N.equivMat dimensionality are inconsistent.')
  }
  if (!all((object@N.equivMat >= -1) | (object@N.equivMat <= 1))) {
    return('N.equivMat entries do not fall within allowed range from -1 to 1.')
  }
  if (!(all(c(ncol(object@N.values),
              ncol(object@N.errors),
              ncol(object@N.z),
              ncol(object@N.sig),
              ncol(object@N.oldValues),
              ncol(object@N.oldErrors),
              ncol(object@N.oldZ),
              ncol(object@N.oldSig)) == object@fpLen))) {
    return('Dimensions of of one of the following slots is inconsistent with fpLen: N.values, N.errors, N.z, N.sig, N.oldValues, N.oldErrors, N.oldZ, N.oldSig.')
  }
  if (!(all(c(nrow(object@N.values),
              nrow(object@N.errors),
              nrow(object@N.z),
              nrow(object@N.sig),
              nrow(object@N.oldValues),
              nrow(object@N.oldErrors),
              nrow(object@N.oldZ),
              nrow(object@N.oldSig)) == 4))) {
    return('Dimensions of of one of the following slots is inconsistent with fpLen: N.values, N.errors, N.z, N.sig, N.oldValues, N.oldErrors, N.oldZ, N.oldSig.')
  }
  if (length(unique(c(dim(object@N.oldValues)[3],
              dim(object@N.oldErrors)[3],
              dim(object@N.oldZ)[3],
              dim(object@N.oldSig)[3]))) != 1) {
    return('Dimensions of of one of the following slots is inconsistent with the others: N.oldValues, N.oldErrors, N.oldZ, N.oldSig.')
  }
  
  return(TRUE)
}





