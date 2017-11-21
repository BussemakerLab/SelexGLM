#' @include N-class.R
NULL

#' Get Updated Regression Formula
#' 
#' Uses design matrix in \code{design} to build regression formula with independent nucleotide features.
#' Additionally, removes features with the same value for every observation in \code{design}.
#' 
#' @param design Design Matrix for regression. 
#' @param model Object of class \linkS4class{model} used to score binding windows and direct analysis. 
#' @return Regression formula of class character. 
#' @export
updatedRegressionFormula = function(design, model) {
  regForm = model@regressionFormula
  if (model@useFixedValuesOffset.N == TRUE) {
    if (all(design$fixedNddG == 0)) {
      fvnVar = "[+]offset(fixedNddG)[+]"
      regForm = gsub(fvnVar, "+", regForm)
      fvnVar2 = "[+]offset(fixedNddG)$"
      regForm = gsub(fvnVar2, "", regForm)
    }
  }
  
  if (model@useFixedValuesOffset.Shape == TRUE) {
    if (all(design$fixedSddG == 0)) {
      fvsVar = "[+]offset(fixedSddG)[+]"
      regForm = gsub(fvsVar, "+", regForm)
      fvsVar2 = "[+]offset(fixedSddG)$"
      regForm = gsub(fvsVar2, "", regForm)
    }
  }
  if (!all(model@features@N@N.set == 0)) {
    sumN = getDesignMatrix(model@features@N, design)
    Ns = c("N.A", "N.C", "N.G", "N.T")
    maxNs = Ns[apply(sumN, 1, which.max)]
    names(maxNs) = rownames(sumN)
    occN = sumN
    occN[occN > 0] = 1
    occLevels = apply(occN, 1, sum)
    names(occLevels) = rownames(sumN)
    for (i in as.numeric(rownames(sumN))) {
      index = as.character(i)
      if (occLevels[names(occLevels) == index] < 2) {
        bs = Ns
      } else {
        bs = c(maxNs[names(maxNs) == index])
        bsZ = Ns[occN[index,] == 0]
        bs = c(bs, bsZ)
      }
      bs = paste(bs, i, sep = "")
      for (b in bs) {
        baseVar = paste("[+]", b,"[+]", sep = "")
        regForm = gsub(baseVar, "+", regForm)
        baseVar2 = paste("[+]", b,"$", sep = "")
        regForm = gsub(baseVar2, "", regForm)
      }
    }
  }
  sumR = getRoundDesign(model@features@Intercept, design) 
  if ((ncol(sumR) == 1) & (nrow(sumR) == 1)) {
    uniqueRound = colnames(sumR)[1]
    r = paste("Round.", uniqueRound, sep = "")
    roundVar = paste("[+]", r,"[+]", sep = "")
    regForm = gsub(roundVar, "+", regForm)
    roundVar2 = paste("[+]", r,"$", sep = "")
    regForm = gsub(roundVar2, "", regForm)
  } else {
    uniqueRounds = colnames(sumR)
    maxR = which.max(sumR[1,])
    r = paste("Round.", names(maxR)[1], sep = "")
    roundVar = paste("[+]", r,"[+]", sep = "")
    regForm = gsub(roundVar, "+", regForm)
    roundVar2 = paste("[+]", r,"$", sep = "")
    regForm = gsub(roundVar2, "", regForm)
    zeroRounds = uniqueRounds[sumR[1,] == 0]
    if (length(zeroRounds) > 0) {
      zeroRounds = paste("Round.", zeroRounds , sep = "")
        for (r in zeroRounds) {
          roundVar = paste("[+]", r,"[+]", sep = "")
          regForm = gsub(roundVar, "+", regForm)
          roundVar2 = paste("[+]", r,"$", sep = "")
          regForm = gsub(roundVar2, "", regForm)
        }
    }
    
  }
  sumI = getDesignMatrix(model@features@Intercept, design)
  
  if (model@includeView == TRUE) {
    valuesI = reshape2::melt(sumI, varnames = c("Strand", "View"),
                             na.rm = FALSE, as.is = FALSE, value.name = "count")
    maxVindex = which.max(valuesI$count)
    zeroVindex = c(1:nrow(valuesI))[valuesI$count == 0]
    vs = c(maxVindex, zeroVindex)
    for (v in vs) {
      vLabel = paste(valuesI$Strand[v], gsub("View.", "", valuesI$View[v]), sep = "")
      viewVar = paste("[+]", vLabel,"[+]", sep = "")
      regForm = gsub(viewVar, "+", regForm)
      viewVar2 = paste("[+]", r,"$", sep = "")
      regForm = gsub(viewVar2, "", regForm)
    }
  } else if (model@includeDNAstrand == TRUE) {
    sumS = matrix(0, nrow = 2, ncol = 1, dimnames = list(rownames(sumI),
                                                         NULL))
    for (st in rownames(sumI)) {
      sumS[st,1] = sum(sumI[st,])
    }
    
    maxSt = which.max(sumS[,1])
    st = rownames(sumS)[maxSt]
    strandVar = paste("[+]", st,"[+]", sep = "")
    regForm = gsub(strandVar, "+", regForm)
    strandVar2 = paste("[+]", st,"$", sep = "")
    regForm = gsub(strandVar2, "", regForm)
  }
  return (regForm)
}



