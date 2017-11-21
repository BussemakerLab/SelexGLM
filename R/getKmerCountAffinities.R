#' Create k-mer affinity table for seeding GLM analysis of SELEX data.
#'
#' @param Sample Sample handle to the dataset on which k-mer counting should be perfomed. 
#' @param k K-mer length to be counted.
#' @param markovModel Markov model handle to use to predict previous round probabilities and expected counts.
#' @param minCount The minimum number of counts for a k-mer to be output.
#' @param symmetrize logical: indicates whether count cutoff and sorting should be imposed on the symmetrized counts/affinities.
#' @return A data frame containing k-mer affinities (standard and reverse complement symmetric versions).
#' @export
getKmerCountAffinities = function(sample, k, markovModel, minCount = 100, symmetrize = TRUE) {
  KmerTable = SELEX::selex.affinities(sample = sample, k = k, minCount = 1, markovModel = markovModel)
  RevCompTable = KmerTable
  KmerTable$Kmer = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(KmerTable$Kmer)))
  MergedTable = merge(KmerTable, RevCompTable, by = "Kmer", all = TRUE, sort = FALSE)
  MergedTable$Kmer = as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(MergedTable$Kmer)))
  names(MergedTable) = c("Kmer", "ObservedCount", "Probability", "ExpectedCount", "Affinity", "SE",
                         "OC.RC", "Prob.RC", "ExpectedCount.RC", "Affinity.RC", "SE.RC")
  MergedTable$ObservedCount[is.na(MergedTable$ObservedCount)] = 0
  MergedTable$Probability[is.na(MergedTable$Probability)] = 0
  MergedTable$ExpectedCount[is.na(MergedTable$ExpectedCount)] = 0
  MergedTable$Affinity[is.na(MergedTable$Affinity)] = 0
  MergedTable$SE[is.na(MergedTable$SE)] = 0
  MergedTable$OC.RC[is.na(MergedTable$OC.RC)] = 0
  MergedTable$Prob.RC[is.na(MergedTable$Prob.RC)] = 0
  MergedTable$ExpectedCount.RC[is.na(MergedTable$ExpectedCount.RC)] = 0
  MergedTable$Affinity.RC[is.na(MergedTable$Affinity.RC)] = 0
  MergedTable$SE.RC[is.na(MergedTable$SE.RC)] = 0
  MergedTable$ObservedCountSym = MergedTable$ObservedCount+MergedTable$OC.RC
  MergedTable$ProbabilitySym = MergedTable$Probability+MergedTable$Prob.RC
  MergedTable$ExpectedCountSym = MergedTable$ExpectedCount+MergedTable$ExpectedCount.RC
  if (symmetrize == TRUE) {
    MergedTable = MergedTable[MergedTable$ObservedCountSym >= minCount,]
    MergedTable$AffinitySym = (MergedTable$ObservedCountSym/MergedTable$ExpectedCountSym)^(1/SELEX::selex.getAttributes(sample)$Round)
    MergedTable$AffinitySym = MergedTable$AffinitySym/max(MergedTable$AffinitySym)
    MergedTable$SESym = 2/sqrt(MergedTable$ObservedCountSym)
    MergedTable$OC.RC = NULL
    MergedTable$Prob.RC = NULL
    MergedTable$ExpectedCount.RC = NULL
    MergedTable$Affinity.RC = NULL
    MergedTable$SE.RC = NULL
    MergedTable = MergedTable[order(-MergedTable$AffinitySym, -MergedTable$Affinity),]
  } else {
    MergedTable = MergedTable[(MergedTable$ObservedCount >= minCount),]
    MergedTable$Affinity = MergedTable$Affinity/max(MergedTable$Affinity)
    MergedTable$AffinitySym = (MergedTable$ObservedCountSym/MergedTable$ExpectedCountSym)^(1/SELEX::selex.getAttributes(sample)$Round)
    MergedTable$AffinitySym = MergedTable$AffinitySym/max(MergedTable$AffinitySym)
    MergedTable$SESym = 2/sqrt(MergedTable$ObservedCountSym)
    MergedTable$OC.RC = NULL
    MergedTable$Prob.RC = NULL
    MergedTable$ExpectedCount.RC = NULL
    MergedTable$Affinity.RC = NULL
    MergedTable$SE.RC = NULL
    MergedTable = MergedTable[order(-MergedTable$Affinity, -MergedTable$AffinitySym),]
  }
  rownames(MergedTable) = NULL
  return (MergedTable)
}



