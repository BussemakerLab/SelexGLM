# Updates values slot in featureSet and its feature-type specific slots (e.g. N, Shape, etc.)

setGeneric("updateValues<-", 
           function(object,value){standardGeneric("updateValues<-")})
setGeneric("updateErrors<-", 
           function(object,value){standardGeneric("updateErrors<-")})
setGeneric("updateZ<-", 
           function(object,value){standardGeneric("updateZ<-")})
setGeneric("updateSig<-", 
           function(object,value){standardGeneric("updateSig<-")})

# Get Generic Functions
setGeneric("getValues", 
           function(object,value){standardGeneric("getValues")})
setGeneric("getErrors", 
           function(object,value){standardGeneric("getErrors")})
setGeneric("getZ", 
           function(object,value){standardGeneric("getZ")})
setGeneric("getSig", 
           function(object,value){standardGeneric("getSig")})
setGeneric("getOldValues", 
           function(object,value){standardGeneric("getOldValues")})
setGeneric("getOldErrors", 
           function(object,value){standardGeneric("getOldErrors")})
setGeneric("getOldZ", 
           function(object,value){standardGeneric("getOldZ")})
setGeneric("getOldSig", 
           function(object,value){standardGeneric("getOldSig")})
setGeneric("getEquivMat", 
           function(object,value){standardGeneric("getEquivMat")})
setGeneric("getSeedLen", 
           function(object,value){standardGeneric("getSeedLen")})
setGeneric("getUpFootprintExtend", 
           function(object,value){standardGeneric("getUpFootprintExtend")})
setGeneric("getFsUpFootprintExtend", 
           function(object,value){standardGeneric("getFsUpFootprintExtend")})
setGeneric("getDownFootprintExtend", 
           function(object,value){standardGeneric("getDownFootprintExtend")})
setGeneric("getFsDownFootprintExtend", 
           function(object,value){standardGeneric("getFsDownFootprintExtend")})
setGeneric("getFpLen", 
           function(object,value){standardGeneric("getFpLen")})
setGeneric("getPositions", 
           function(object,value){standardGeneric("getPositions")})
setGeneric("getNumViews", 
           function(object,value){standardGeneric("getNumViews")})
setGeneric("getRcSymmetric", 
           function(object,value){standardGeneric("getRcSymmetric")})
setGeneric("getRounds", 
           function(object,value){standardGeneric("getRounds")})
setGeneric("getShapeParamsUsed", 
           function(object,value){standardGeneric("getShapeParamsUsed")})
setGeneric("getName", 
           function(object,value){standardGeneric("getName")})
setGeneric("getVarRegLen", 
           function(object,value){standardGeneric("getVarRegLen")})
setGeneric("getLeftFixedSeq", 
           function(object,value){standardGeneric("getLeftFixedSeq")})
setGeneric("getRightFixedSeq", 
           function(object,value){standardGeneric("getRightFixedSeq")})
setGeneric("getLeftFixedSeqOverlap", 
           function(object,value){standardGeneric("getLeftFixedSeqOverlap")})
setGeneric("getRightFixedSeqOverlap", 
           function(object,value){standardGeneric("getRightFixedSeqOverlap")})
setGeneric("getConsensusSeq", 
           function(object,value){standardGeneric("getConsensusSeq")})
setGeneric("getAffinityType", 
           function(object,value){standardGeneric("getAffinityType")})
setGeneric("getConfidenceLevel", 
           function(object,value){standardGeneric("getConfidenceLevel")})
setGeneric("getMinAffinity", 
           function(object,value){standardGeneric("getMinAffinity")})
setGeneric("getMissingValueSuppression", 
           function(object,value){standardGeneric("getMissingValueSuppression")})
setGeneric("getMinSeedValue", 
           function(object,value){standardGeneric("getMinSeedValue")})
setGeneric("getRegressionFormula", 
           function(object,value){standardGeneric("getRegressionFormula")})
setGeneric("getIteration", 
           function(object,value){standardGeneric("getIteration")})
setGeneric("getIncludeDNAstrand", 
           function(object,value){standardGeneric("getIncludeDNAstrand")})
setGeneric("getIncludeView", 
           function(object,value){standardGeneric("getIncludeView")})
setGeneric("getIncludeShape", 
           function(object,value){standardGeneric("getIncludeShape")})
setGeneric("getUseFixedValuesOffset.N", 
           function(object,value){standardGeneric("getUseFixedValuesOffset.N")})
setGeneric("getUseFixedValuesOffset.Shape", 
           function(object,value){standardGeneric("getUseFixedValuesOffset.Shape")})
setGeneric("getExUpBases", 
           function(object,value){standardGeneric("getExUpBases")})
setGeneric("getExDownBases", 
           function(object,value){standardGeneric("getExDownBases")})
setGeneric("getShapeTable", 
           function(object,value){standardGeneric("getShapeTable")})
setGeneric("getFeatures", 
           function(object,value){standardGeneric("getFeatures")})
setGeneric("getVerbose", 
           function(object,value){standardGeneric("getVerbose")})
setGeneric("getDateInitialized", 
           function(object,value){standardGeneric("getDateInitialized")})
setGeneric("getDateLastModified", 
           function(object,value){standardGeneric("getDateLastModified")})
setGeneric("getDateBetasLastAdded", 
           function(object,value){standardGeneric("getDateBetasLastAdded")})
setGeneric("getGlmFits", 
           function(object,value){standardGeneric("getGlmFits")})
setGeneric("getDesignMatrixSummary", 
           function(object,value){standardGeneric("getDesignMatrixSummary")})
setGeneric("getN",
           function(object,value){standardGeneric("getN")})
setGeneric("getIntercept",
           function(object,value){standardGeneric("getIntercept")})
setGeneric("getShape",
           function(object,value){standardGeneric("getShape")})
setGeneric("getDesignMatrix",
           function(object,design){standardGeneric("getDesignMatrix")})
setGeneric("getFeatureDesign",
           function(object,...){standardGeneric("getFeatureDesign")})
setGeneric("getPlotValues",
           function(object,...){standardGeneric("getPlotValues")})
setGeneric("getPSAM",
           function(object, ...){standardGeneric("getPSAM")})
setGeneric("getTopSeqLabels",
           function(object, ...){standardGeneric("getTopSeqLabels")})



# Set Generic Functions
setGeneric("setValues<-", 
           function(object,value){standardGeneric("setValues<-")})
setGeneric("setErrors<-", 
           function(object,value){standardGeneric("setErrors<-")})
setGeneric("setZ<-", 
           function(object,value){standardGeneric("setZ<-")})
setGeneric("setSig<-", 
           function(object,value){standardGeneric("setSig<-")})
setGeneric("setOldValues<-", 
           function(object,value){standardGeneric("setOldValues<-")})
setGeneric("setOldErrors<-", 
           function(object,value){standardGeneric("setOldErrors<-")})
setGeneric("setOldZ<-", 
           function(object,value){standardGeneric("setOldZ<-")})
setGeneric("setOldSig<-", 
           function(object,value){standardGeneric("setOldSig<-")})
setGeneric("setEquivMat<-", 
           function(object,value){standardGeneric("setEquivMat<-")})
setGeneric("setSeedLen<-", 
           function(object,value){standardGeneric("setSeedLen<-")})
setGeneric("setUpFootprintExtend<-", 
           function(object,value){standardGeneric("setUpFootprintExtend<-")})
setGeneric("setFsUpFootprintExtend<-", 
           function(object,value){standardGeneric("setFsUpFootprintExtend<-")})
setGeneric("setDownFootprintExtend<-", 
           function(object,value){standardGeneric("setDownFootprintExtend<-")})
setGeneric("setFsDownFootprintExtend<-", 
           function(object,value){standardGeneric("setFsDownFootprintExtend<-")})
setGeneric("setFpLen<-", 
           function(object,value){standardGeneric("setFpLen<-")})
setGeneric("setPositions<-", 
           function(object,value){standardGeneric("setPositions<-")})
setGeneric("setNumViews<-", 
           function(object,value){standardGeneric("setNumViews<-")})
setGeneric("setRcSymmetric<-", 
           function(object,value){standardGeneric("setRcSymmetric<-")})
setGeneric("setRounds<-", 
           function(object,value){standardGeneric("setRounds<-")})
setGeneric("setShapeParamsUsed<-", 
           function(object,value){standardGeneric("setShapeParamsUsed<-")})
setGeneric("setName<-", 
           function(object,value){standardGeneric("setName<-")})
setGeneric("setVarRegLen<-", 
           function(object,value){standardGeneric("setVarRegLen<-")})
setGeneric("setLeftFixedSeq<-", 
           function(object,value){standardGeneric("setLeftFixedSeq<-")})
setGeneric("setRightFixedSeq<-", 
           function(object,value){standardGeneric("setRightFixedSeq<-")})
setGeneric("setLeftFixedSeqOverlap<-", 
           function(object,value){standardGeneric("setLeftFixedSeqOverlap<-")})
setGeneric("setRightFixedSeqOverlap<-", 
           function(object,value){standardGeneric("setRightFixedSeqOverlap<-")})
setGeneric("setConsensusSeq<-", 
           function(object,value){standardGeneric("setConsensusSeq<-")})
setGeneric("setAffinityType<-", 
           function(object,value){standardGeneric("setAffinityType<-")})
setGeneric("setConfidenceLevel<-", 
           function(object,value){standardGeneric("setConfidenceLevel<-")})
setGeneric("setMinAffinity<-", 
           function(object,value){standardGeneric("setMinAffinity<-")})
setGeneric("setMissingValueSuppression<-", 
           function(object,value){standardGeneric("setMissingValueSuppression<-")})
setGeneric("setMinSeedValue<-", 
           function(object,value){standardGeneric("setMinSeedValue<-")})
setGeneric("setRegressionFormula<-", 
           function(object,value){standardGeneric("setRegressionFormula<-")})
setGeneric("setIteration<-", 
           function(object,value){standardGeneric("setIteration<-")})
setGeneric("setIncludeDNAstrand<-", 
           function(object,value){standardGeneric("setIncludeDNAstrand<-")})
setGeneric("setIncludeView<-", 
           function(object,value){standardGeneric("setIncludeView<-")})
setGeneric("setIncludeShape<-", 
           function(object,value){standardGeneric("setIncludeShape<-")})
setGeneric("setUseFixedValuesOffset.N<-", 
           function(object,value){standardGeneric("setUseFixedValuesOffset.N<-")})
setGeneric("setUseFixedValuesOffset.Shape<-", 
           function(object,value){standardGeneric("setUseFixedValuesOffset.Shape<-")})
setGeneric("setExUpBases<-", 
           function(object,value){standardGeneric("setExUpBases<-")})
setGeneric("setExDownBases<-", 
           function(object,value){standardGeneric("setExDownBases<-")})
setGeneric("setShapeTable<-", 
           function(object,value){standardGeneric("setShapeTable<-")})
setGeneric("setFeatures<-", 
           function(object,value){standardGeneric("setFeatures<-")})
setGeneric("setVerbose<-", 
           function(object,value){standardGeneric("setVerbose<-")})
setGeneric("setDateInitialized<-", 
           function(object,value){standardGeneric("setDateInitialized<-")})
setGeneric("setDateLastModified<-", 
           function(object,value){standardGeneric("setDateLastModified<-")})
setGeneric("setDateBetasLastAdded<-", 
           function(object,value){standardGeneric("setDateBetasLastAdded<-")})
setGeneric("setGlmFits<-", 
           function(object,value){standardGeneric("setGlmFits<-")})
setGeneric("setDesignMatrixSummary<-", 
           function(object,value){standardGeneric("setDesignMatrixSummary<-")})
setGeneric("finalizeFeatureBetas",
            function(object){standardGeneric("finalizeFeatureBetas")})
setGeneric("setN<-",
           function(object,value){standardGeneric("setN<-")})
setGeneric("setIntercept<-",
           function(object,value){standardGeneric("setIntercept<-")})
setGeneric("setShape<-",
           function(object,value){standardGeneric("setShape<-")})


setGeneric("buildSymmetricEquivalenceMatrix",
           function(object){standardGeneric("buildSymmetricEquivalenceMatrix")})
setGeneric("buildSymSquaredEquivalenceMatrix",
           function(object){standardGeneric("buildSymSquaredEquivalenceMatrix")})
setGeneric("buildNullEquivalenceMatrix",
           function(object){standardGeneric("buildNullEquivalenceMatrix")})
setGeneric("formatEquivalenceMatrix",
           function(object){standardGeneric("formatEquivalenceMatrix")})




setGeneric("addSeedPsam<-",
           function(object, value){standardGeneric("addSeedPsam<-")})

setGeneric (
  name = "seedTable2psam",
  def = function(object, seedTable){standardGeneric("seedTable2psam")}
)
setGeneric("addNewBetas",
           function(object, design, value,...){standardGeneric("addNewBetas")})


setGeneric("plotZs",
           function(object, ...){standardGeneric("plotZs")})

