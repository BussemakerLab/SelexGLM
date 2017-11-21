
#' Validity Function for \linkS4class{Intercept}
#' 
#' Validity Function for object of class \linkS4class{Intercept}.
#' 
#' @param object Object of class \linkS4class{Intercept}.
#' @export
validIntercept = function(object) {
  if (length(unique(c(2,
                      dim(object@I.oldValues)[1],
                   dim(object@I.oldErrors)[1],
                   dim(object@I.oldZ)[1],
                   dim(object@I.oldSig)[1],
                   dim(object@I.values)[1],
                   dim(object@I.errors)[1],
                   dim(object@I.z)[1],
                   dim(object@I.sig)[1]))) != 1) {
    return('I.values, I.errors, I.z, I.sig, I.oldValues, I.oldErrors, I.oldZ, and I.oldSig must all have 2 rows.')
  }

  if (length(unique(c(dim(object@I.oldValues)[2],
                      dim(object@I.oldErrors)[2],
                      dim(object@I.oldZ)[2],
                      dim(object@I.oldSig)[2],
                      dim(object@I.values)[2],
                      dim(object@I.errors)[2],
                      dim(object@I.z)[2],
                      dim(object@I.sig)[2]))) != 1) {
    return('I.values, I.errors, I.z, I.sig, I.oldValues, I.oldErrors, I.oldZ, and I.oldSig must all have the same number of columns (i.e. views).')
  }

  if (length(unique(c(dim(object@I.oldValues)[3],
                      dim(object@I.oldErrors)[3],
                      dim(object@I.oldZ)[3],
                      dim(object@I.oldSig)[3],
                      dim(object@I.values)[3],
                      dim(object@I.errors)[3],
                      dim(object@I.z)[3],
                      dim(object@I.sig)[3]))) != 1) {
    return('I.values, I.errors, I.z, I.sig, I.oldValues, I.oldErrors, I.oldZ, and I.oldSig must all have the same length of dimension 3 (i.e. number of rounds).')
  }

  if (length(unique(c(dim(object@I.oldValues)[4],
                      dim(object@I.oldErrors)[4],
                      dim(object@I.oldZ)[4],
                      dim(object@I.oldSig)[4]))) != 1) {
    return('I.oldValues, I.oldErrors, I.oldZ, and I.oldSig must all have the same length of dimension 4 (i.e. number of previous iterations).')
  }

  return (TRUE)
}

