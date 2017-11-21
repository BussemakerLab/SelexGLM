#' @include topModelMatch.R
#' @include addDesignMatrix.R
#' @include updatedRegressionFormula.R
NULL

#' Iterates Selex GLM Algorithm
#' 
#' Loops through SelexGLM algorithm until stable model is reached. 
#' 
#' @param model Object of class \linkS4class{model}.
#' @param initalData Table of initial probe counts.
#' @export
iterateSelexGLM = function(model, initialData) {
  data = topModelMatch(initialData, model)
  data = addDesignMatrix(data, model)
  designMatrixSummary = getDesignMatrix(model, data)
  if (model@verbose == TRUE) {
    print("\n")
    print("Round summary: ")
    print (designMatrixSummary$Round)
    print("\n")
    print("Mono-nucleotide summary: ")
    print (designMatrixSummary$N)
    print("\n")
    print("View/strand orientation summary: ")
    print (designMatrixSummary$Intercept)
  }
  regressionFormula = updatedRegressionFormula(data, model)
  if (model@verbose == TRUE) {
    print("\n")
    print("Regression Formula: ")
    print (regressionFormula)
  }
  fit = glm(regressionFormula, 
            data=data, 
            family = poisson(link="log"))
  model = addNewBetas(model, data, fit)
  data = topModelMatch(initialData, model)
  data = addDesignMatrix(data, model)
  designMatrixSummary.v2 = getDesignMatrix(model, data)
  if ((all(designMatrixSummary.v2$N == designMatrixSummary$N)) & (all(designMatrixSummary.v2$Round == designMatrixSummary$Round)) & (all(designMatrixSummary.v2$Intercept == designMatrixSummary$Intercept)))  {
    print ("Stability Reached")
  }
  for (i in 2:20) {
    if (data.nrow == nrow(data)) {
      break
    }
    data.nrow = nrow(data)
    print (paste("i =",i))
    
    designMatrixSummary = getDesignMatrix(model, data)
    if (model@verbose == TRUE) {
      print("\n")
      print("Round summary: ")
      print (designMatrixSummary$Round)
      print("\n")
      print("Mono-nucleotide summary: ")
      print (designMatrixSummary$N)
      print("\n")
      print("View/strand orientation summary: ")
      print (designMatrixSummary$Intercept)
    }
    
    regressionFormula = updatedRegressionFormula(data, model)
    if (model@verbose == TRUE) {
      print("\n")
      print("Regression Formula: ")
      print (regressionFormula)
    }
    fit = glm(regressionFormula, 
              data=data, 
              family = poisson(link="log"))
    summary(fit)
    model = addNewBetas(model, data, fit)
    data = topModelMatch(data, model)
    data = addDesignMatrix(data, model)
    if (model@verbose == TRUE) {
      print(paste("Number of Observations in Design Matrix: ",nrow(data), sep = ""))
    }
    designMatrixSummary.v2 = getDesignMatrix(model, data)
    
    if ((all(designMatrixSummary.v2$N == designMatrixSummary$N)) & (all(designMatrixSummary.v2$Round == designMatrixSummary$Round)) & (all(designMatrixSummary.v2$Intercept == designMatrixSummary$Intercept)))  {
      if (model@verbose == TRUE) {
        print (paste("Stability Reached after ", i, " iterations.", sep = ""))
      }
      break
    } else if (nrow(data) == 0) {
      if (model@verbose == TRUE) {
        print ("Algorithm failed to converge: No probes meet the confidence level requirement (Confidence Level:", model@confidenceLevel, ")", sep = "")
      }
    }
  }
  model <- finalizeFeatureBetas(model)
  return (model)
}