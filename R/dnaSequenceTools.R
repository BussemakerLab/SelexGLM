
#' Convert DNA iupac sequence to regular expression
#' 
#' Converts an iupac sequence to R regular expression. If the sequence is already in regular expression form (i.e. it only contains the letters A, C, G, T with brackets and numbers), it will not be changed. 
#' 
#' @param refSeq Sequence in iupac or regular expression format.
#' @return Sequence in regular expression format recognized by R.
#' @export
convertSeq2RegEx = function(refSeq) {
  refSeqCheck = gsub("A", "", refSeq)
  refSeqCheck = gsub("C", "", refSeqCheck)
  refSeqCheck = gsub("G", "", refSeqCheck)
  refSeqCheck = gsub("T", "", refSeqCheck)
  if (nchar(refSeqCheck) != 0) {
    refSeq = gsub("N", "[ACGT]", refSeq)
    refSeq = gsub("B", "[CGT]", refSeq)
    refSeq = gsub("D", "[AGT]", refSeq)
    refSeq = gsub("H", "[ACT]", refSeq)
    refSeq = gsub("V", "[ACG]", refSeq)
    refSeq = gsub("R", "[AG]", refSeq)
    refSeq = gsub("Y", "[CT]", refSeq)
    refSeq = gsub("K", "[GT]", refSeq)
    refSeq = gsub("M", "[AC]", refSeq)
    refSeq = gsub("S", "[CG]", refSeq)
    refSeq = gsub("W", "[AT]", refSeq)
  }
  return (refSeq)
}

#' Gets reverse complement of DNA sequence in iupac or regular expression format.
#' 
#' Returns reverse complement of DNA sequence . Extends 'Biostrings' functionality to include regular expression characters and iupac symbols not equal to A, C, G, and T.
#' 
#' @param refSeq Sequence in iupac or regular expression format.
#' @return Reverse complement of input in the same format as the input (i.e. regular expression or iupac).
#' @export
revComp = function(refSeq) {
  refSeq = gsub("A", "F", refSeq)
  refSeq = gsub("T", "A", refSeq)
  refSeq = gsub("F", "T", refSeq)
  
  refSeq = gsub("C", "F", refSeq)
  refSeq = gsub("G", "C", refSeq)
  refSeq = gsub("F", "G", refSeq)
  
  refSeq = gsub("R", "F", refSeq)
  refSeq = gsub("Y", "R", refSeq)
  refSeq = gsub("F", "Y", refSeq)
  
  refSeq = gsub("K", "F", refSeq)
  refSeq = gsub("M", "K", refSeq)
  refSeq = gsub("F", "M", refSeq)
  
  refSeq = gsub("B", "F", refSeq)
  refSeq = gsub("V", "B", refSeq)
  refSeq = gsub("F", "V", refSeq)
  
  
  refSeq = gsub("D", "F", refSeq)
  refSeq = gsub("H", "D", refSeq)
  refSeq = gsub("F", "H", refSeq)
  
  refSeq = paste(rev(strsplit(refSeq, NULL)[[1]]), collapse="")
  refSeqList = strsplit(refSeq, split = character(0))[[1]]
  for (i in 1:length(refSeqList)) {
    if (refSeqList[i] == "]") {
      refSeqList[i] = "X"
    } 
  }
  for (i in 1:length(refSeqList)) {
    if (refSeqList[i] == "[") {
      refSeqList[i] = "]"
    } 
  }
  for (i in 1:length(refSeqList)) {
    if (refSeqList[i] == "X") {
      refSeqList[i] = "["
    } 
  }
  refSeq = paste(refSeqList, collapse= "") 
  
  return (refSeq)
}


#' Count sequence length of a DNA regular expression.
#' 
#' Count number of characters implied by a regular expression for DNA sequences.
#' 
#' @param reSeq Sequence in regular expression format (i.e. containing A,C,G,T,'[',']','{','}').
#' @return Length of DNA sequence implied by reSeq input.
#' @export
countDNAseqLength = function(reSeq) {
  splitRE = strsplit(reSeq, NULL)[[1]]
  l.splitRE = length(splitRE)
  i = 1
  x = 0
  for (loop in 1:l.splitRE) {
    addNum = 1
    if (splitRE[i] == "[") {
      for (j in (i+1):l.splitRE) {
        addNum = addNum + 1
        if (splitRE[j] == "]") {
          if (j != l.splitRE) {
            if (splitRE[j+1] == "{") {
              addNum = addNum+1
              for (k in (j+2):l.splitRE) {
                addNum = addNum+1
                if (splitRE[k] == "}") {
                  x = x+as.numeric(paste(splitRE[(j+2):(k-1)], collapse = ""))
                  break
                  
                }
              }
              break
            } else {
              x = x+1
              break
            }
          } else {
            x = x+1
            break
          }
        }
      }
    } else {
      x = x+1
    }
    
    i = i+addNum
    if (i > l.splitRE) {break}
  }
  return (x)
}





