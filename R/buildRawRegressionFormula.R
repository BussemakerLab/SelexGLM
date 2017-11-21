
#' Build Raw Regression Formula
#' 
#' Builds a raw regression formula for an object of class \linkS4class{model} based on its feature design.
#' 
#' This regression formula must be updated using the data's design matrix to exclude features
#' in the model that do not appear in the aligned data and/or a single mono-nucleotide at each position
#' such that the data is not overdescribed.
#' 
#' @param model Object of class \linkS4class{model}.
#' @return Raw regression formula containing features of \linkS4class{model} object.
#' @export
buildRawRegressionFormula = function(model) {
  regForm = "ObservedCount ~ offset(logProb)+"
  if ((model@useFixedValuesOffset.N == TRUE) & (!all(c(1:model@fpLen) %in% model@features@N@N.set))) {
    regForm = paste(regForm,"offset(fixedNddG)+", sep = "")
  }
  if ((model@useFixedValuesOffset.Shape == TRUE) & (!all(c(1:model@fpLen) %in% model@features@Shape@Shape.set))) {
    regForm = paste(regForm,"offset(fixedSddG)+", sep = "")
  }
  roundVars = paste(paste("Round.", model@rounds[[1]], sep = ""), collapse = "+")
  regForm = paste(regForm, "+", roundVars, "+", sep = "")
  if (length(model@features@N@N.set[model@features@N@N.set > 0]) > 0) {
    if (all(model@features@N@N.equivMat==0)) {
      Params = expand.grid(c("N.A", "N.C", "N.G", "N.T"), model@features@N@N.set[model@features@N@N.set > 0])
      
    } else {
      N.set = model@features@N@N.set[model@features@N@N.set > 0]
      N.equivMat = model@features@N@N.equivMat
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
      
      
      Params = expand.grid(c("N.A", "N.C", "N.G", "N.T"), N.included)
      Params.selfSym = expand.grid(c("N.A", "N.C"), N.selfSym)
      Params = rbind(Params, Params.selfSym)
      
    }
    regVars = paste(as.character(Params$Var1), as.character(Params$Var2), sep = "")
    regForm = paste(regForm, paste(regVars, collapse="+"), sep="")
  }
  
  
  if (model@includeView == TRUE){
    if (model@rcSymmetric == FALSE) {
      Params = expand.grid(paste("Strand.", c("F", "R"), sep = ""), c(1:model@numViews))
      regVars = paste(as.character(Params$Var1), as.character(Params$Var2), sep = "")
      regForm = paste(regForm, "+",paste(regVars, collapse="+"), sep="")
    } else {
      regVars = paste("Strand.F", c(1:model@numViews), sep = "")
      regForm = paste(regForm, "+",paste(regVars, collapse="+"), sep="")
    }
  } else if (model@includeDNAstrand == TRUE) {
    regVars = paste("Strand.", c("F", "R"), sep = "")
    regForm = paste(regForm, "+",paste(regVars, collapse="+"), sep="" )
  }
  if (model@includeShape == TRUE) {
    if (length(model@features@Shape@Shape.set[model@features@Shape@Shape.set > 0]) > 0) {
      if (all(model@features@Shape@Shape.equivMat==0)) {
        
        Params = expand.grid(rownames(model@features@Shape@Shape.values), model@features@Shape@Shape.set)
        
      } else {
        Shape.set = model@features@Shape@Shape.set
        Shape.equivMat = model@features@Shape@Shape.equivMat
        Shape.eMsub = Shape.equivMat[Shape.set,]
        Shape.eMsub = Shape.eMsub[,Shape.set]
        rownames(Shape.eMsub) = Shape.set
        colnames(Shape.eMsub) = Shape.set
        
        Shape.included = c(numeric(0))
        Shape.selfSym= c(numeric(0))
        for (i in Shape.set) {
          if (length(which(Shape.eMsub[,as.character(i)] == 0)) == nrow(Shape.eMsub)) {
            Shape.included = c(Shape.included, i)
          } else if (length(which(Shape.eMsub[,as.character(i)] == 0)) == (nrow(Shape.eMsub)-1)) {
            if (Shape.eMsub[as.character(i), as.character(i)] == -1) {
              Shape.selfSym = c(Shape.selfSym, i)
            }
          }
        }
        
        Params = expand.grid(rownames(model@features@Shape@Shape.values), c(Shape.included))
        if (length(Shape.selfSym) > 0) {
          sVals = rownames(model@features@Shape@Shape.values)
          sVals = sVals[!grepl("B", sVals)]
          Params2 = expand.grid(sVals, c(Shape.selfSym))
          Params = rbind(Params, Params2)
        }
      }
      regVars = paste(as.character(Params$Var1), as.character(Params$Var2), sep = "")
      regForm = paste(regForm, "+",paste(regVars, collapse="+"), sep="")
    }
  }
  
  
  
  regForm = gsub("[+]{2}", "+", regForm)
  regForm = gsub("[+]$", "", regForm)
  
  return(regForm)
}