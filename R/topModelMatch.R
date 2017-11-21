#' Finds Top Match View on Probe for Model. 
#' 
#' Selects top binding window for each probe and assigns model affinity score and confidence level for selected window.
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
#' @return A data frame containing probes with top binding windows selected, along with their PSAM scores and confidence levels.
#' @export
topModelMatch = function(data, model) {
  data = data[,c("Probe", "ObservedCount", "Probability", "Round")]
  concatLeft.start = nchar(model@leftFixedSeq)-(model@leftFixedSeqOverlap+model@exUpBases)+1
  concatRight.end = model@rightFixedSeqOverlap+model@exDownBases
  addLeft = stringi::stri_sub(model@leftFixedSeq, concatLeft.start, nchar(model@leftFixedSeq))
  addRight = stringi::stri_sub(model@rightFixedSeq, 1, concatRight.end)
  data$extendLmer = paste(addLeft, data$Probe, addRight, sep = "")
  data$RevComp = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(data$extendLmer)))
  data$Lmer = stringi::stri_sub(data$extendLmer, model@exUpBases+1, nchar(data$extendLmer)-model@exDownBases)
  extendLmer.length = nchar(data$extendLmer[1])
  ddGcols = paste("ddG", 1:((1+1*(model@rcSymmetric == FALSE))*model@numViews), sep = "")
  for (i in 1:length(ddGcols)) {
    data[,ddGcols[i]] = 0
  }
  data = addViewddG(data, model)
  data = addNddG(data, model)
  data = addShapeddG(data, model)
  
  data = ddg2affinity(data, model)
  
  affcols = paste("Affinity", 1:((1+1*(model@rcSymmetric == FALSE))*model@numViews), sep = "")
  if (model@verbose == TRUE){
    print ("Finding Top Affinity Binding Views")
  }
  topMatchFrame = apply(data[,affcols], 1, which.max)
  topMatchRelAff = apply(data[,affcols], 1, max)
  confLevelDenom = apply(data[,affcols]^data$Round, 1, sum)
  topMatchConfidence = ((topMatchRelAff)^data$Round)/confLevelDenom
  data = cleanAffinityInfo(data, model)
  # ALIGNED PROBE
  lenAlignedSequence = model@fpLen+model@exUpBases+model@exDownBases
  data$alignedFootprint = stringi::stri_sub(data$extendLmer, topMatchFrame, (topMatchFrame+lenAlignedSequence-1))
  data$alignedFootprint[topMatchFrame > model@numViews] = 
    stringi::stri_sub(data$RevComp[topMatchFrame > model@numViews], 
             (topMatchFrame[topMatchFrame > model@numViews]-model@numViews),
             (topMatchFrame[topMatchFrame > model@numViews]-model@numViews+lenAlignedSequence-1))
  
  data$topMatchSequence = stringi::stri_sub(data$alignedFootprint, model@exUpBases+1, model@exUpBases+model@fpLen)
  data$topMatchRelAff = topMatchRelAff
  data$topMatchConfidence = topMatchConfidence
  data$topMatchView = topMatchFrame
  data$topMatchStrand = "F"
  data$topMatchStrand[data$topMatchView>model@numViews] = "R"
  data$topMatchView[data$topMatchStrand == "R"] = data$topMatchView[data$topMatchStrand == "R"]-model@numViews
  data$topMatchStrand = factor(data$topMatchStrand)
  data$extendLmer = NULL
  data$RevComp = NULL
  
  rownames(data) = NULL
  return(data)
}

#' Adds nucleotide ddG values for views of DNA probe.
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}. 
addNddG = function(data, model) {
  extendLmer.length = nchar(data$extendLmer[1])
  if (model@verbose == TRUE) {
    print ("Adding Nucleotide Design Info")
  }
  for (i in 1:((1+1*(model@rcSymmetric == FALSE))*(extendLmer.length-model@exUpBases-model@exDownBases))) {
    if (model@verbose == TRUE) {
      print (paste("Probe Index", as.character(i),"of",as.character((1+1*(model@rcSymmetric == FALSE))*(extendLmer.length-model@exUpBases-model@exDownBases)), sep = " "))
    }
    data = addNDesign(data, model, i)
  }
  N.mat = model@features@N@N.values
  rownames(N.mat) = NULL
  colnames(N.mat) = NULL
  flatCM = as.vector(N.mat)
  names(flatCM) = NULL
  if (model@verbose == TRUE) {
    print ("Calculating Nucleotide ddG values")
  }
  for (j in 1:((1+1*(model@rcSymmetric == FALSE))*(model@numViews))) {
    if (model@verbose == TRUE) {
      print (paste("View ", as.character(j),"of",as.character(((1+1*(model@rcSymmetric == FALSE))*(model@numViews))), sep = " "))
    }
    data = addDdg.N.frame(data, model, flatCM, j)
  }
  if (model@verbose == TRUE) {
    print ("Cleaning Nucleotide Design")
  }
  for (i in 1:((1+1*(model@rcSymmetric == FALSE))*(extendLmer.length-model@exUpBases-model@exDownBases))) {
    data = cleanNDesign(data, model, i)
  } 
  return (data)
}



#' Parses probes into nucleotide feature columns to be used in calculating ddG values for binding windows. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
#' @param frame Position in DNA probe.
addNDesign = function(data, model, frame) {
  extendLmer.length = nchar(data$extendLmer[1])
  if ((frame) <= (extendLmer.length-model@exDownBases-model@exUpBases)) {
    seqCol = "extendLmer"
    index = frame+model@exUpBases
  } else {
    seqCol = "RevComp"
    index = frame-(extendLmer.length-model@exDownBases-model@exUpBases)+model@exUpBases
  }
  
  Alab = paste("N.A", frame, sep = "")
  Clab = paste("N.C", frame, sep = "")
  Glab = paste("N.G", frame, sep = "")
  Tlab = paste("N.T", frame, sep = "")
  
  data[,Alab] = 0
  data[(stringi::stri_sub(data[,seqCol], index, index) == "A"), Alab] = 1
  data[,Alab] = data[,Alab]*data$Round
  data[,Clab] = 0
  data[(stringi::stri_sub(data[,seqCol], index, index) == "C"), Clab] = 1
  data[,Clab] = data[,Clab]*data$Round
  data[,Glab] = 0
  data[(stringi::stri_sub(data[,seqCol], index, index) == "G"), Glab] = 1
  data[,Glab] = data[,Glab]*data$Round
  data[,Tlab] = 0
  data[(stringi::stri_sub(data[,seqCol], index, index) == "T"), Tlab] = 1
  data[,Tlab] = data[,Tlab]*data$Round
  
  return(data) 
}

#' Calculates nucleotide ddG values from nucleotide design and beta values for views of DNA probe. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
#' @param flatCM Vector version of N.values. Built from as.vector() of model@features@N@N.values with rownames and colnames set to NULL.
#' @param window Binding view on DNA probe. 
addDdg.N.frame = function(data, model, flatCM, window) {
  ddGlabel = paste("ddG", as.character(window), sep = "")
  if (window <= model@numViews) {
    matchNum = window
  } else {
    if (model@rcSymmetric == TRUE) {
      return (data)
    }
    matchNum = nchar(data$extendLmer[1])-model@exUpBases-model@exDownBases+window-model@numViews
  }
  
  aMatchVar = paste("N.A", matchNum, sep = "")
  
  cIndex.start = match(aMatchVar, colnames(data))
  cIndex.end = cIndex.start+4*(model@fpLen)-1
  rownames(data) = NULL
  ddGhold = as.matrix(data[,cIndex.start:cIndex.end])%*%cbind(flatCM[1:length(flatCM)])
  data[,ddGlabel] = ddGhold/data$Round
  return(data)
}

#' Removes nucleotide feature columns used in calculating ddG values for binding windows. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
#' @param frame Position in DNA probe.
cleanNDesign = function(data, model, frame) {
  if (model@rcSymmetric == TRUE) {
    if (frame > (nchar(data$extendLmer[1])-model@exDownBases-model@exUpBases)) {
      return (data)
    }
  }
  Alab = paste("N.A", frame, sep = "")
  Clab = paste("N.C", frame, sep = "")
  Glab = paste("N.G", frame, sep = "")
  Tlab = paste("N.T", frame, sep = "")
  
  
  data[,Alab] = NULL
  data[,Clab] = NULL
  data[,Glab] = NULL
  data[,Tlab] = NULL
  
  return(data) 
}


#' Adds shape ddG values for views of DNA probe.
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
addShapeddG = function(data, model) {
  if (model@includeShape) {
    extendLmer.length = nchar(data$extendLmer[1])
    if (model@verbose == TRUE) {
      print ("Adding Shape Design Info")
    }
    for (i in 1:((1+1*(model@rcSymmetric == FALSE))*(extendLmer.length-model@exUpBases-model@exDownBases))) {
      if (model@verbose == TRUE) {
        print (paste("Frame Index", as.character(i),"of",as.character(((1+1*(model@rcSymmetric == FALSE))*(extendLmer.length-model@exUpBases-model@exDownBases))), sep = " "))
      }
      data = addShapeDesign(data, model, i)
    }
    
    if (model@verbose == TRUE) {
      print ("Adding Shape ddG values")
    }
    for (j in 1:((1+1*(model@rcSymmetric == FALSE))*(model@numViews))) {
      if (model@verbose == TRUE) {
        print (paste("View", as.character(j),"of",as.character(((1+1*(model@rcSymmetric == FALSE))*model@numViews)), sep = " "))
      }
      data = addDdg.Shape.frame(data, model, j)
    }
    if (model@verbose == TRUE) {
      print ("Cleaning Shape Design")
    }
    for (i in 1:((1+1*(model@rcSymmetric == FALSE))*(extendLmer.length-model@exUpBases-model@exDownBases))) {
      data = cleanShapeDesign(data, model, i)
    } 
  }
  
  return (data)
}

#' Parses probes into shape feature columns to be used in calculating ddG values for binding windows. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
#' @param frame Position in DNA probe. 
addShapeDesign = function(data, model, frame) {
  extendLmer.length = nchar(data$extendLmer[1])
  shapeLengthAdd = model@exUpBases+model@exDownBases
  if ((frame) <= (extendLmer.length-model@exDownBases-model@exUpBases)) {
    seqCol = "extendLmer"
    index = frame
  } else {
    seqCol = "RevComp"
    index = frame-(extendLmer.length-model@exDownBases-model@exUpBases)
  }
  for (s in rownames(model@features@Shape@Shape.values)) {
    sLabel = paste(s, frame, sep = "")
    sCol = gsub("Shape.", "", s)
    data[,sLabel] = model@shapeTable[stringi::stri_sub(data[,seqCol], index, (index+shapeLengthAdd)),sCol]*data$Round
  }
  return(data) 
}



#' Adds window-specific shape ddG contributions for each binding window of the probe.  
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
#' @param window Binding view on DNA probe. 
addDdg.Shape.frame = function(data, model, window) {
  lmer.length = nchar(data$extendLmer[1])-model@exDownBases-model@exUpBases
  Shape.mat = model@features@Shape@Shape.values
  rownames(Shape.mat) = NULL
  colnames(Shape.mat) = NULL
  ddGlabel = paste("ddG", as.character(window), sep = "")
  if (window <= model@numViews) {
    flatCM = as.vector(Shape.mat)
    matchNum = window
    matchShape = rownames(model@features@Shape@Shape.values)[1]
    matchVar = paste(matchShape, matchNum, sep = "")
    cIndex.start = match(matchVar, colnames(data))
    cIndex.end = cIndex.start+length(rownames(model@features@Shape@Shape.values))*(model@fpLen)-1
  } else {
    flatCM = as.vector(Shape.mat)
    matchNum = lmer.length+(window-(model@numViews))
    matchShape = rownames(model@features@Shape@Shape.values)[1]
    matchVar = paste(matchShape, matchNum, sep = "")
    cIndex.start = match(matchVar, colnames(data))
    cIndex.end = cIndex.start+length(rownames(model@features@Shape@Shape.values))*(model@fpLen)-1
  }
  rownames(data) = NULL
  ddGhold = as.matrix(data[,cIndex.start:cIndex.end])%*%cbind(flatCM[1:length(flatCM)])
  data[,ddGlabel] = data[,ddGlabel]+ddGhold/data$Round
  return(data)
}

#' Removes shape feature columns used in calculating ddG values for binding windows. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
#' @param frame Position in DNA probe.
cleanShapeDesign = function(data, model, frame) {
  for (s in rownames(model@features@Shape@Shape.values)) {
    sLabel = paste(s, frame, sep = "")
    data[,sLabel] = NULL
  }
  return(data) 
}



#' Adds intercept type, view-specific, strand-specific, and round-specific ddG values. 
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}. 
addViewddG = function(data, model) {
  if (model@verbose == TRUE) {
    print ("Adding Intercept ddG Values")
  }
  Intercept = getIntercept(model)
  ddGviews = c(1:((1+1*(model@rcSymmetric == FALSE))*model@numViews))
  for (v in ddGviews) {
    ddGlabel = paste("ddG", v, sep = "")
    if (v <= model@numViews) {
      data[,ddGlabel] = data[,ddGlabel] + Intercept@I.values[1,v,paste("Round.", data$Round, sep = "")]
    } else {
      data[,ddGlabel] = data[,ddGlabel] + Intercept@I.values[2,(2*model@numViews-v+1),paste("Round.", data$Round, sep = "")]
    }
  }
  return (data)
}

#' Converts view-specific ddG values to affinity scores.
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}.
ddg2affinity = function(data, model) {
  maxW = (1+1*(model@rcSymmetric == FALSE))*model@numViews
  ddGcols = paste("ddG", 1:maxW, sep = "")
  max.ddG = by(data[,ddGcols], data$Round, max)
  for (l in unique(data$Round)) {
    data[data$Round == l,ddGcols] = data[data$Round == l,ddGcols]-max.ddG[as.character(l)]
  }
  for (i in 1:maxW) {
    ddG.lab = paste("ddG", i, sep = "")
    affinity.lab = paste("Affinity", i, sep = "")
    data[,ddG.lab] = exp(data[,ddG.lab])
    colnames(data)[colnames(data)==ddG.lab] = affinity.lab
  }
  return(data)
}


#' Removes view-specific affinity scores from data frame. To be used after selection of probe-specific top binding window.
#' 
#' @param data Count table for DNA probes used for filtering. 
#' @param model Object of class \linkS4class{model}. 
cleanAffinityInfo = function(data, model) {
  maxW = 2*model@numViews
  if (model@rcSymmetric == TRUE) {maxW = model@numViews}
  for (i in 1:maxW) {
    affinity.lab = paste("Affinity", i, sep = "")
    data[,affinity.lab] = NULL
  }
  return(data)
}


