require(parallel)
require(e1071)
options(stringsAsFactors = FALSE)

load("data/gbm.validated.RData")
genes <- as.character(unique(c(gbm.validated[,1], gbm.validated[,2])))
load("data/gbm.gene.expr.RData")
load("data/gbm.svms.RData")
svms <- gbm.svms

##########################################################
##  CERNIA main function do the following:
##
##  1. Launch DT-Hybrid
##  2. Compute scores for the selected pairs
##  3. Perform the classification and collect the results
##########################################################

  
cat("Loaded", length(genes), "genes, continue...\n")

if(!exists("cerna.recommendations")) source("classification/RunDTHybrid.R")

# Collect the recommendations
cerna.recommendations <- launch.DTHybrid("BRCA")
cerna.recommendations <- cerna.recommendations[genes, genes]

source("classification/ScoringFunction.R")
source("classification/TissueClassification.R")

new.data  <- final.scores$scores.matrix.noNA
classification.results <- cerna.classification.scores(classify(new.data, svms))
write.table(classification.results, file = "classification.results.txt", quote = FALSE, row.names = FALSE, sep = "\t")

