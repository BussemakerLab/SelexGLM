#' @include topModelMatch.R
NULL

#' Score K-mer Sequences
#' 
#' Scores Kmers using model beta values for all feature parameters estimated with regression.
#' 
#' @param data Table of k-mers.
#' @param model Object of class \linkS4class{model}.
#' @param l.index Left-most position of PSAM to be used for scoring,
#' @param PSAMcol Name for PSAM affinity column. 
#' @param seqCol Name for k-mer variable to be scored.
#' @export
scoreKmers = function(data, model, l.index, PSAMcol = "PredictedAffinity", seqCol = "Kmer") {
  PSAM = getPSAM(model)
  nc = ncol(PSAM)
  k = nchar(data[1,seqCol])
  if (nc == k) {
    if (!(is.null(l.index))) {
      if (l.index != 1) {
        stop('l.index must be equal to 1 for PSAM length = table k-mer length.')
      }
    }
  } else if (nc < k) {
    stop('PSAM length must be >= table k-mer length.')
  } else {
    if (!is.null(l.index)) {
      if (l.index+k > nc) {
        stop('l.index+table k-mer > PSAM length. Choose a smaller value of l.index.')
      } 
    } else {
      if (nc %% 2 == k %% 2) {
        l.index = (nc-k)/2+1
      } else {
        stop('l.index must be defined explicitly for PSAM and data table kmer lengths of different parity.')
      }
    }
  }
  PSAM = PSAM[,l.index:(l.index+k-1)]
  for (i in 1:k) {
    Alab = paste("N.A", i, sep = "")
    Clab = paste("N.C", i, sep = "")
    Glab = paste("N.G", i, sep = "")
    Tlab = paste("N.T", i, sep = "")
    data[,Alab] = 0
    data[(stringi::stri_sub(data[,seqCol], i, i) == "A"), Alab] = 1
    data[,Clab] = 0
    data[(stringi::stri_sub(data[,seqCol], i, i) == "C"), Clab] = 1
    data[,Glab] = 0
    data[(stringi::stri_sub(data[,seqCol], i, i) == "G"), Glab] = 1
    data[,Tlab] = 0
    data[(stringi::stri_sub(data[,seqCol], i, i) == "T"), Tlab] = 1
  }
  N.mat = log(PSAM)
  rownames(N.mat) = NULL
  colnames(N.mat) = NULL
  flatCM = as.vector(N.mat)
  names(flatCM) = NULL
  ddGlabel = "ddG"
  aMatchVar = "N.A1"
  cIndex.start = match(aMatchVar, colnames(data))
  cIndex.end = cIndex.start+4*k-1
  rownames(data) = NULL
  ddGhold = as.matrix(data[,cIndex.start:cIndex.end])%*%cbind(flatCM[1:length(flatCM)])
  data[,ddGlabel] = ddGhold
  for (i in 1:k) {
    Alab = paste("N.A", i, sep = "")
    Clab = paste("N.C", i, sep = "")
    Glab = paste("N.G", i, sep = "")
    Tlab = paste("N.T", i, sep = "")
    data[,Alab] = NULL
    data[,Clab] = NULL
    data[,Glab] = NULL
    data[,Tlab] = NULL
  } 
  
  
  if (model@includeShape) {
    if (model@exUpBases > 0) {
      for (w in 1:model@exUpBases) {
        for (x in rownames(model@features@Shape@Shape.values)) {
          sLabel = paste(x, w, sep = "")
          data[,sLabel] = 0
        }
      }
    }
    for (j in (1+model@exUpBases):(k-model@exDownBases)) {
      shapeLengthAdd = model@exUpBases+model@exDownBases
      for (s in rownames(model@features@Shape@Shape.values)) {
        sLabel = paste(s,j, sep = "")
        sP = gsub("Shape.", "", s)
        data[,sLabel] = model@shapeTable[stringi::stri_sub(data$Kmer, j-model@exUpBases, j+model@exUpBases),sP]
      }
    }
    if (model@exDownBases > 0) {
      for (y in (k-model@exDownBases+1):k) {
        for (z in rownames(model@features@Shape@Shape.values)) {
          sLabel = paste(z, y, sep = "")
          data[,sLabel] = 0
        }
      }
    }
    Shape.mat = model@features@Shape@Shape.values
    rownames(Shape.mat) = NULL
    colnames(Shape.mat) = NULL
    flatCM = as.vector(Shape.mat)
    matchNum = 1
    matchShape = rownames(model@features@Shape@Shape.values)[1]
    matchVar = paste(matchShape, matchNum, sep = "")
    cIndex.start = match(matchVar, colnames(data))
    cIndex.end = cIndex.start+nrow(model@features@Shape@Shape.values)*(model@fpLen)-1
    ddGhold = as.matrix(data[,cIndex.start:cIndex.end])%*%cbind(flatCM[1:length(flatCM)])
    data[,ddGlabel] = data[,ddGlabel]+ddGhold
    for (l in 1:k) {
      for (x in rownames(model@features@Shape@Shape.values)) {
        sLabel = paste(x, l, sep = "")
        data[,sLabel] = NULL
      }
    }
  }
  
  data[,PSAMcol] = exp(data$ddG)
  if (model@includeShape == TRUE) {
    data[,PSAMcol] = data[,PSAMcol]/max(data[,PSAMcol])
  }

  data$ddG = NULL
  data = data[order(-data[,PSAMcol]),]
  rownames(data) = NULL
  return(data)
}





