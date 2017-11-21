#' Build Probe Count Table
#' 
#' Create probe count table suitable as starting point for GLM analysis of SELEX data.
#' 
#' @param SampleA SELEX sample handle used to define the probe 'universe'. 
#' If no \code{SampleB} argument is provided, this sample is first randomly split into two 
#' equal-sized halves,which are then used to define the counts and the probe universe, 
#' respectively.
#' @param SampleB SELEX Sample handle used to define the probe counts (optional).
#' @param markovModel Markov model for rounds 0.
#' @return A data frame containing counts for each probe in the universe.
#' @export
getProbeCounts = function(SampleA, markovModel, SampleB) {
  k.A = as.numeric(as.character(SELEX::selex.getAttributes(SampleA)$VariableRegionLength))
  if (missing(SampleB)) {
    r.split = SELEX::selex.split(sample = SampleA)
    rA= r.split$train
    rB = r.split$test
    CountTableA = SELEX::selex.counts(sample = rA, k =k.A, minCount = 1, markovModel = markovModel)
    CountTableB = SELEX::selex.counts(sample = rB, k = k.A, minCount = 1, markovModel = markovModel)
  } else {
    CountTableA = SELEX::selex.counts(sample = SampleA, k =k.A, minCount = 1, markovModel = markovModel)
    CountTableB = SELEX::selex.counts(sample = SampleB, k = k.A, minCount = 1, markovModel = markovModel)
  }
  MergeB = data.frame(Kmer=CountTableB$Kmer, OCtest=CountTableB$ObservedCount)
  tA = merge(CountTableA, MergeB, by="Kmer", all.x=TRUE, sort=FALSE)
  tA[is.na(tA$OCtest), "ObservedCount"] = 0
  tA$OCtest = NULL
  tA = tA[order(tA$ObservedCount, decreasing=TRUE), ]
  tA$Round = SELEX::selex.getAttributes(SampleA)$Round
  tA$ExpectedCount = NULL
  colnames(tA)[colnames(tA) == "Kmer"] = "Probe"
  return(tA)
}
