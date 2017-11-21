#' @include topModelMatch.R
NULL

#' Add Design Matrix
#' 
#' Uses top binding windows to align probes and add corresponding design matrix to \code{\link{topModelMatch}} for regression based on feature design in \linkS4class{model} object.
#'
#' Builds a design matrix for probes in \code{data} using the features contained in \code{model}. The following features will be added for use in the glm regression:
#' \enumerate{
#'  \item \code{logProb}: log(\code{data$Probability})
#'  \item \code{N.Ax, N.Cx, N.Gx, N.Tx}: A nucleotide at position x is represented by a set of 4 nucleotide variables that are equal to 0 or \code{data$Round} depending on the identity of the nucleotide in question. If reverse complement symmetry is being used, contributions from x and its reverse complement symmetric equivalent nucleotide, x', are both represented by a single block of 4 nucleotides.
#'  \item \code{fixedNddG}: If useFixedValuesOffset.N == TRUE, fixedNddG is the ddG contribution from nucleotides not included in \code{model@features@N@N.set}. 
#'  \item \code{Round.R}: The round of a probe, r, is represented by a value of 1 in column \code{Round.r} and 0 in all other round columns.
#'  \item \code{Strand.[FR]}: If includeDNAstrand == TRUE in \code{model} input, the DNA strand of the top view is represented by (\code{data$Strand.F}, \code{data$Strand.R}) equal to 0 or \code{data$Round} depending on which strand the top view is located on. 
#'  \item \code{Strand.[FR]v}: If includeView == TRUE in \code{model} input, the View and strand of the top view is represented by a set of variables \code{data$Strand.[FR]v}, where v is in [1, (number of views per strand)]. These variables are equal to \code{data$Round} if the top view is the one represented by that variable, or zero if not. 
#'  \item \code{Shape.[shape]x}: A shape parameter at position x is represented by a single variable \code{Shape.[shape]x}, which is equal to the shapeTable lookup value for the shape parameter at position x. A variable of this kind will be included for each shape in \code{model@shapeParamsUsed[[1]]} and each position in \code{model@features@Shape@Shape.set}
#' }
#' 
#' @param data Table output by \code{\link{topModelMatch}} containing probes, probe counts and information about the top binding view for each probe.
#' @param model Object of class \linkS4class{model} used to build design matrix from the output of \code{\link{topModelMatch}}.
#'
#' @return Data frame containing probes with design matrix built for regression input.
#' @export
addDesignMatrix = function(data, model) {
  # Filters out probes that don't meet the confidence level minimum.  
  data = data[data$topMatchConfidence >= model@confidenceLevel,]
  if (nrow(data) == 0) {
    warning('No observations meet the confidence level minimum.')
    return(data)
  }
  
  # Filters out probes that don't meet the affinity minimum.  
  data = data[data$topMatchRelAff >= model@minAffinity,]
  if (nrow(data) == 0) {
    warning('No observations meet the confidence level and view affinity minima.')
    return(data)
  }
  
  data$extendLmer = data$alignedFootprint
  data$logProb = log(data$Probability)
  
  # Adds nucleotide design to 'data'
  data = addDesignMatrix.N(data, model)
  
  # Adds round design to 'data'
  data = addDesignMatrix.Round(data, model)
  
  # Adds view design to 'data'
  data = addDesignMatrix.View(data, model)
  
  # Adds shape design to 'data'
  data = addDesignMatrix.Shape(data, model)
  data$extendLmer = NULL
  rownames(data) = NULL
  
  return (data)
}


#' Adds design matrix for nucleotide features. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
addDesignMatrix.N = function(data, model) {  
  if (model@verbose == TRUE) {
    print ("Adding Nucleotide Design Info")
  }
  
  # Performs the following if at least one nucleotide is included in N.set for fitting:
  if (!all(model@features@N@N.set == 0)) {
    for (i in model@features@N@N.set) {
      # adds N.Ai, N.Ci, N.Gi, N.Ti to 'data' using addNDesign from topModelMatch.R
      data = addNDesign(data, model, i)
    }
    
    # Performs the following if N.equivMat has non-zero entries, such that the same betas describe nucleotides at mutliple positions
    if (!(all(model@features@N@N.equivMat == 0))) {
      equivMat = model@features@N@N.equivMat
      for (i in rev(model@features@N@N.set)) {
        matchCol = equivMat[,i]
        matchIndices = which(matchCol != 0)
        if (length(matchIndices[matchIndices %in% model@features@N@N.set]) > 0) {
          matchIndex = min(matchIndices[matchIndices %in% model@features@N@N.set])
          # If the matched entry has  value of 1, position 'i' and position 'matchIndex' should have identical betas (N.A(i) = N.A(matchIndex), etc.)
          if (matchCol[matchIndex] == 1) {
            N = matchIndex
            N.equiv = i
            # The nucleotide variables at position i are added to the nucleotide variables at position 'matchIndex'.
            data[,paste("N.A", N, sep = "")] = data[,paste("N.A", N, sep = "")]+data[,paste("N.A", N.equiv, sep = "")]
            data[,paste("N.C", N, sep = "")] = data[,paste("N.C", N, sep = "")]+data[,paste("N.C", N.equiv, sep = "")] 
            data[,paste("N.G", N, sep = "")] = data[,paste("N.G", N, sep = "")]+data[,paste("N.G", N.equiv, sep = "")] 
            data[,paste("N.T", N, sep = "")] = data[,paste("N.T", N, sep = "")]+data[,paste("N.T", N.equiv, sep = "")]
            # The nucleotide variables from position i are removed from 'data'.
            data[,paste("N.A", N.equiv, sep = "")] = NULL
            data[,paste("N.C", N.equiv, sep = "")] = NULL
            data[,paste("N.G", N.equiv, sep = "")] = NULL
            data[,paste("N.T", N.equiv, sep = "")] = NULL
            equivMat[i,] = 0
          # If the matched entry has  value of -1, position 'i' and position 'matchIndex' should have reverse complement symmetrized betas (N.A(i) = N.T(matchIndex), etc.)
          } else if (matchCol[matchIndex] == -1) {
            # For the case where i is not equal to matchIndex, symmetrizing proceeds as follows:
            if (matchIndex != i) {
              N = matchIndex
              N.RC = i
              # The nucleotide variables at position i are added to the reverse complement symmetric nucleotide variables at position 'matchIndex'.
              data[,paste("N.A", N, sep = "")] = data[,paste("N.A", N, sep = "")]+data[,paste("N.T", N.RC, sep = "")]
              data[,paste("N.C", N, sep = "")] = data[,paste("N.C", N, sep = "")]+data[,paste("N.G", N.RC, sep = "")] 
              data[,paste("N.G", N, sep = "")] = data[,paste("N.G", N, sep = "")]+data[,paste("N.C", N.RC, sep = "")] 
              data[,paste("N.T", N, sep = "")] = data[,paste("N.T", N, sep = "")]+data[,paste("N.A", N.RC, sep = "")]
              # The nucleotide variables from position i are removed from 'data'.
              data[,paste("N.A", N.RC, sep = "")] = NULL
              data[,paste("N.C", N.RC, sep = "")] = NULL
              data[,paste("N.G", N.RC, sep = "")] = NULL
              data[,paste("N.T", N.RC, sep = "")] = NULL
              equivMat[i,] = 0
            # For the case where i is equal to matchIndex, symmetrizing proceeds as follows:
            } else {
              N = matchIndex
              # The nucleotide variables N.A, N.C at position i are added to N.T, N.G, respectively, at position
              data[,paste("N.A", N, sep = "")] = data[,paste("N.A", N, sep = "")]+data[,paste("N.T", N, sep = "")]
              data[,paste("N.C", N, sep = "")] = data[,paste("N.C", N, sep = "")]+data[,paste("N.G", N, sep = "")]
              # The nucleotide variables N.T, N.G from position i are removed from 'data'.
              data[,paste("N.G", N, sep = "")] = NULL
              data[,paste("N.T", N, sep = "")] = NULL
            }
          }
        }
      }
    }
  }
  
  # Adds fixedNddG variable if useFixedValueOffset.N== TRUE, and N.set does not include all positions in the model footprint.
  if ((model@useFixedValuesOffset.N == TRUE) & (length(model@features@N@N.set) < model@fpLen)) {
    # N.notSet is the set of positions within the model's footprint that are not contained in N.set
    N.notSet = c(1:model@fpLen)[!(c(1:model@fpLen) %in% model@features@N@N.set)]
    # adds N.Ai, N.Ci, N.Gi, N.Ti to 'data' using addNDesign from topModelMatch.R for i in N.notSet
    for (i in N.notSet) {
      data = addNDesign(data, model, i)
    }
    # Adds corresponding ddG value by dotting in beta values for N.notSet entries from N.values.
    N.mat = model@features@N@N.values[,N.notSet]
    rownames(N.mat) = NULL
    colnames(N.mat) = NULL
    flatCM = as.vector(N.mat)
    names(flatCM) = NULL
    matchNum = min(N.notSet)
    aMatchVar = paste("N.A", matchNum, sep = "")
    
    cIndex.start = match(aMatchVar, colnames(data))
    cIndex.end = cIndex.start+4*(length(N.notSet))-1
    rownames(data) = NULL
    ddGhold = as.matrix(data[,cIndex.start:cIndex.end])%*%cbind(flatCM[1:length(flatCM)])
    data$fixedNddG = ddGhold
    
    # Removes Nucleotide variables for N.notSet using cleanNDesign from topModelMatch.R
    for (i in N.notSet) {
      data = cleanNDesign(data, model, i)
    }
  }
  
  return (data)
}

#' Adds design matrix for round features. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
addDesignMatrix.Round = function(data, model) {  
  if (model@verbose == TRUE) {
    print ("Adding Round Design Info")
  }
  # identifies the set of rounds occuring in data after filtering.
  uniqueRounds = model@rounds[[1]]
  
  # Adds Round.r variable for each r in uniqueRounds
  for (uR in uniqueRounds) {
    roundLabel = paste("Round.", uR, sep = "")
    data[,roundLabel] = 0
    data[(data$Round == uR),roundLabel] = 1
  }
  return (data)
}

#' Adds design matrix for view features. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
addDesignMatrix.View = function(data, model) {
  if ((model@includeView == FALSE) & (model@includeDNAstrand == FALSE)) {
    return (data)
  }
  if (model@verbose == TRUE) {
    print ("Adding View Design Info")
  }
  sList = c("F", "R")[1:(1+1*(model@rcSymmetric == FALSE))]
  # If includeView == TRUE, Adds Strand.[FR]v variables for all possible views. 
  if (model@includeView == TRUE) {
    for (s in sList) {
      for (v in 1:model@numViews) {
        vLabel = paste("Strand.", s, v, sep = "")
        data[,vLabel] = 0
        data[((data$topMatchStrand == s) & (data$topMatchView == v)), vLabel] = 1*data$Round[((data$topMatchStrand == s) & (data$topMatchView == v))]
      }
    }
  # If includeDNAstrand == TRUE, Adds Strand.[FR] variables for all possible views. 
  } else if (model@includeDNAstrand == TRUE) {
    for (s in sList) {
      vLabel = paste("Strand.", s, sep = "")
      data[,vLabel] = 0
      data[(data$topMatchStrand == s), vLabel] = 1*data$Round[(data$topMatchStrand == s)]
    }
  }
  
  return (data)
}

#' Adds design matrix for shape features. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
addDesignMatrix.Shape = function(data, model) {
  if (model@includeShape == FALSE) {
    return (data)
  }
  if (model@verbose == TRUE) {
    print ("Adding Shape Design Info")
  }
  if (!all(model@features@Shape@Shape.set == 0)) {
    for (i in model@features@Shape@Shape.set) {
      if (model@verbose == TRUE) {
        print (paste("Probe Index", as.character(i),"of",length(model@features@Shape@Shape.set), sep = " "))
      }
      data = addShapeDesign(data, model, i)
    }
    if (!(all(model@features@Shape@Shape.equivMat == 0))) {
      equivMat = model@features@Shape@Shape.equivMat
      for (i in rev(model@features@Shape@Shape.set)) {
        matchCol = equivMat[,i]
        matchIndices = which(matchCol != 0)
        if (length(matchIndices[matchIndices %in% model@features@Shape@Shape.set]) > 0) {
          matchIndex = min(matchIndices[matchIndices %in% model@features@Shape@Shape.set])
          if (matchCol[matchIndex] == 1) {
            S = matchIndex
            S.equiv = i
            if (S != S.equiv) {
              for (s in rownames(model@features@Shape@Shape.values)) {
                data[,paste(s, S, sep = "")] = data[,paste(s, S, sep = "")]+data[,paste(s, S.equiv, sep = "")]
                data[,paste(s, S.equiv, sep = "")] = NULL
              }
            }
            equivMat[i,] = 0
          } else if (matchCol[matchIndex] == -1) {
            if (matchIndex != i) {
              S = matchIndex
              S.RC= i
              shapeParams = model@shapeParamsUsed[[1]]
              if ("MGW" %in% shapeParams) {
                data[,paste("Shape.MGW", S, sep = "")] = data[,paste("Shape.MGW", S, sep = "")]+data[,paste("Shape.MGW", S.RC, sep = "")]
                data[,paste("Shape.MGW", S.RC, sep = "")] = NULL
              }
              if ("ProT" %in% shapeParams) {
                data[,paste("Shape.ProT", S, sep = "")] = data[,paste("Shape.ProT", S, sep = "")]+data[,paste("Shape.ProT", S.RC, sep = "")]
                data[,paste("Shape.ProT", S.RC, sep = "")] = NULL
              }
              if ("HelT" %in% shapeParams) {
                data[,paste("Shape.HelTA", S, sep = "")] = data[,paste("Shape.HelTA", S, sep = "")]+data[,paste("Shape.HelTB", S.RC, sep = "")]
                data[,paste("Shape.HelTB", S.RC, sep = "")] = NULL
                data[,paste("Shape.HelTB", S, sep = "")] = data[,paste("Shape.HelTB", S, sep = "")]+data[,paste("Shape.HelTA", S.RC, sep = "")]
                data[,paste("Shape.HelTA", S.RC, sep = "")] = NULL
              }
              if ("Roll" %in% shapeParams) {
                data[,paste("Shape.RollA", S, sep = "")] = data[,paste("Shape.RollA", S, sep = "")]+data[,paste("Shape.RollB", S.RC, sep = "")]
                data[,paste("Shape.RollB", S.RC, sep = "")] = NULL
                data[,paste("Shape.RollB", S, sep = "")] = data[,paste("Shape.RollB", S, sep = "")]+data[,paste("Shape.RollA", S.RC, sep = "")]
                data[,paste("Shape.RollA", S.RC, sep = "")] = NULL
              }
              equivMat[i,] = 0
            } else {
              S = i
              shapeParams = model@shapeParamsUsed[[1]]
              if ("HelT" %in% shapeParams) {
                data[,paste("Shape.HelTA", S, sep = "")] = data[,paste("Shape.HelTA", S, sep = "")]+data[,paste("Shape.HelTB", S, sep = "")]
                data[,paste("Shape.HelTB", S, sep = "")] = NULL
              }
              if ("Roll" %in% shapeParams) {
                data[,paste("Shape.RollA", S, sep = "")] = data[,paste("Shape.RollA", S, sep = "")]+data[,paste("Shape.RollB", S, sep = "")]
                data[,paste("Shape.RollB", S, sep = "")] = NULL
              }
            }
          }
        }
      }
    }
  }
  if ((model@useFixedValuesOffset.Shape == TRUE) & (length(model@features@Shape@Shape.set) < model@fpLen) & (length(model@shapeParamsUsed[[1]]) > 0)) {
    Shape.notSet = c(1:model@fpLen)[!(c(1:model@fpLen) %in% model@features@Shape@Shape.set)]
    for (i in Shape.notSet) {
      data = addShapeDesign(data, model, i)
    }
    S.mat = model@features@Shape@Shape.values[,Shape.notSet]
    rownames(S.mat) = NULL
    colnames(S.mat) = NULL
    flatCM = as.vector(S.mat)
    names(flatCM) = NULL
    matchNum = min(Shape.notSet)
    aMatchVar = paste(rownames(model@features@Shape@Shape.values)[1],matchNum, sep = "")
    
    cIndex.start = match(aMatchVar, colnames(data))
    cIndex.end = cIndex.start+nrow(model@features@Shape@Shape.values)*(length(Shape.notSet))-1
    rownames(data) = NULL
    ddGhold = as.matrix(data[,cIndex.start:cIndex.end])%*%cbind(flatCM[1:length(flatCM)])
    data$fixedSddG = ddGhold
    for (i in Shape.notSet) {
      data = cleanShapeDesign(data, model, i)
    }
  }
  
  return (data)
}

