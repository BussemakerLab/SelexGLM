## ------------------------------------------------------------------------
library(SELEX)
library(stringi)
library(Biostrings)
library(SelexGLM)
library(devtools)
library(reshape2)
library(ggplot2)
library(Rmisc)

## ---- message=FALSE, warning=FALSE, results='hide'-----------------------
options(java.parameters = "-Xmx4000M")
workDir = tempdir()
selex.config(workingDir=workDir, maxThreadNumber=4)

## ------------------------------------------------------------------------
selex.loadAnnotation(system.file("extdata", "config.xml", package="SELEX"))
selex.sampleSummary()

## ------------------------------------------------------------------------
r0.train = selex.sample(seqName = 'R0.libraries', sampleName='R0.barcodeGC', round = 0)
r0.test = selex.sample(seqName = 'R0.libraries', sampleName='R0.barcodeCG', round = 0)
dataSample = selex.sample(seqName = 'R2.libraries', sampleName = 'ExdHox.R2', round = 2)

## ------------------------------------------------------------------------
# MARKOV MODEL BUILT
kmax = selex.kmax(sample = r0.test)
# Train Markov model on Hm 16bp library Round 0 data
mm = selex.mm(sample = r0.train, order = NA, crossValidationSample =r0.test, Kmax = kmax, mmMethod = "TRANSITION")
mmscores = selex.mmSummary(sample = r0.train)
ido = which(mmscores$R==max(mmscores$R))
mm.order = mmscores$Order[ido]

## ------------------------------------------------------------------------
# INFOGAIN USED TO CALCULATE KLEN
libLen = as.numeric(as.character(selex.getAttributes(dataSample)$VariableRegionLength))
selex.infogain(sample = dataSample, k = c((mm.order+1):libLen), markovModel = mm)
infoscores = selex.infogainSummary(sample = dataSample)

#information gain barplot
idx = which(infoscores$InformationGain==max(infoscores$InformationGain))
colstring = rep('BLUE', nrow(infoscores))
colstring[idx] = 'RED'
barplot(height=infoscores$InformationGain, names.arg=infoscores$K, col=colstring,
        xlab="Oligonucleotide Length (bp)", ylab="Information Gain (bits)")
kLen = infoscores$K[idx]

## ------------------------------------------------------------------------
# For the sake of previous analysis on the Hox data used in this example, I will use kLen.f = 12 as my k-mer length, even though kLen identified through the information gain analysis has kLen = 13
data.kmerTable = selex.affinities(sample=dataSample, k=kLen, markovModel=mm)
data.kmerTable = data.kmerTable[order(-data.kmerTable$Affinity), ]
rownames(data.kmerTable) = NULL

data.probeCounts = getProbeCounts(dataSample, markovModel = mm)
summary(data.probeCounts)
print(data.probeCounts[1:10,])

## ------------------------------------------------------------------------
# Inputs about library are data specific 
model = new("model",
             varRegLen = libLen,
             leftFixedSeq =  "GTTCAGAGTTCTACAGTCCGACGATCTGG", 
             rightFixedSeq ="CCAGCTGTCGTATGCCGTCTTCTGCTTG", 
             seedLen = kLen, 
             leftFixedSeqOverlap = 4,
             initialAffinityCutoff = 0.00,
             missingValueSuppression = 1,
             minSeedValue = .001, 
             upFootprintExtend = 2,
             includeWindowFactor = FALSE,
             confidenceLevel = .95, 
             verbose = FALSE, 
             useFixedValuesOffset.N = FALSE,
             rounds = list(c(2)),
             rcSymmetric = FALSE,
             minAffinity = 0.01
          )

## ------------------------------------------------------------------------
model@features@N

## ------------------------------------------------------------------------
# Model nucleotide Betas before seed PSAM is added
addSeedPsam(model) = seedTable2psam(model, data.kmerTable)
# Model nucleotide Betas after seed PSAM is added
model@features@N

## ------------------------------------------------------------------------
#Use this definition of data for complete analysis
data = data.probeCounts

data = topModelMatch(data, model)
# Uses aligned probes to build design matrix
data = addDesignMatrix(data, model)
# Constructs regression expression with independent features using design matrix
regressionFormula = updatedRegressionFormula(data, model)
fit = glm(regressionFormula, 
          data=data, 
          family = poisson(link="log"))

model = addNewBetas(model, data, fit)
# Nucleotide Features after first round of fitting

# GABRIELLA: this plotting commmand is not working, can you fix it?
#plot(model, Nplot.ddG = TRUE, verticalPlots = TRUE)

data = data.probeCounts

data.nrow = nrow(data)
for (i in 2:3) {
  data = topModelMatch(data, model)
  data = addDesignMatrix(data, model)
  if (data.nrow == nrow(data)) {
    print ("Stability Reached")
    break
  } else {
    data.nrow = nrow(data)
  }

          
  regressionFormula = updatedRegressionFormula(data, model)
  fit = glm(regressionFormula, 
            data=data, 
            family = poisson(link="log"))

  model = addNewBetas(model,data,fit)
  # Nucleotide Features after i'th round of fitting
}
model@features@N@N.values

## ------------------------------------------------------------------------
save(model, file = "HowToFitMononucleotideModel.Result.RData")

