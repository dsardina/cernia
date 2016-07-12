
source("classification/DTHybrid/Recommendation.R")

load("data/interactions.RData")

present <- genes %in% interactions$Target
if(any(!present)) stop(sum(!present), "genes are not in the miRNA-target dataset. Execution stopped!")

################################################
##  Load the interactions and create the matrix
################################################

launch.DTHybrid <- function(tissue.type){
  cat("Loading miRNA-target interactions...\n")
  prj.file <- NULL
  
  if(tissue.type == "BRCA") prj.file <- "data/dty/projectionMtx.brca.cernia.expr.RData"
  else if(tissue.type == "PRAD") prj.file <- "data/dty/projectionMtx.prad.cernia.expr.RData"
  else if(tissue.type == "GBM") prj.file <- "data/dty/projectionMtx.gbm.cernia.expr.RData"
  else{
    cat("Tissue type not valid. Set default type: BRCA")
    prj.file <- "data/dty/projectionMtx.brca.cernia.expr.RData"
  }
  
  if(file.exists(prj.file)){
    cat("Found recommendations in the folder, loading recommendations for", tissue.type,"...")
    load(prj.file)
  }
  else{
    if(!exists("interactions")) stop("Cannot load miRNA-target interactions from data directory.")
    
    tar <- NULL
    if(tissue.type == "BRCA"){
      load("data/dty/brca.genes.RData")
      tar <- brca.genes
    }
    else if(tissue.type == "PRAD"){
      load("data/dty/prad.genes.RData")
      tar <- prad.genes
    }
    else if(tissue.type == "GBM"){
      load("data/dty/gbm.genes.RData")
      tar <- gbm.genes
    }
    else {
      cat("Tissue type not valid. Set default type: BRCA")
      load("data/dty/brca.genes.RData")
      tar <- brca.genes
    }
    
    # Extract miRs and their targets
    mir <- unique(interactions[,1])
    
    # Create the matrix of the interactions
    cat("Creating matrix from interactions...\n")
    A <- matrix(nrow = length(tar), ncol = length(mir), data = 0)
    colnames(A) <- mir
    rownames(A) <- tar
    
    for(i in 1:nrow(interactions)){
      A[as.character(interactions[i,2]), as.character(interactions[i,1])] <- 1
    }
    
    #####################
    ##  Launch DT-Hybrid
    #####################
    
    cat("Launching DT-Hybrid...\n")
    
    # Make projection from bipartite network using DT-hybrid sources
    cl <- makeCluster(detectCores()-1)
    M <- computeRecommendation(A, cl = cl)
    W <- graphWeights(nrow(M), ncol(M), M, cl = cl)
    
    stopCluster(cl)
    cat("Done!\n")
    
    save(W, file = "DTY_recommendations.RData")
  }
  
  return(W)
}