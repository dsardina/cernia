
##################################################
##  Utility function for correlation computation
##################################################

# Generate all pairs from a gene's vector
generate.pairs <- function(genes){
  return(unlist(lapply(1:length(genes), function(i, l){
    if(i != l) return(paste0(genes[i], "-", genes[(i+1):l]))
  }, length(genes)))
  )
}

# Compute the correlation between a pair of gene expressions
computePairCorrelation <- function(x, y, gene.expr){
  
  if(!is.character(x) || !is.character(y)) stop("Error: x and y must be character")
  
  # Get the expression for current target and its miRNA
  xexpr <- gene.expr[x,]
  yexpr <- gene.expr[y,]
  
  return(cor(as.numeric(xexpr), as.numeric(yexpr), method = "pearson", use = "pairwise"))
}

# Compute pairwise correlations
compute.pairwise.correlations <- function(pairs, gene.expr, cl){
  parLapply(cl = cl, pairs, function(x){
    pair <- unlist(strsplit(x, split = "-", fixed = TRUE))
    return(computePairCorrelation(pair[1], pair[2], gene.expr))
  })
}

# Create a matrix from correlation results
correlations.matrix <- function(genes, corrs){
  m <- matrix(data = 0, nrow = length(genes), ncol = length(genes))
  rownames(m) <- genes
  colnames(m) <- genes
  
  lapply(1:length(corrs), function(i){
    pair <- unlist(strsplit(names(corrs)[i], split = "-", fixed = TRUE))
    m[pair[1], pair[2]] <<- corrs[i]
  })
  
  return(m)
}


################################################
##  Compute correlations for all pairs of genes
################################################

compute.correlations <- function(genes, gene.expr){
  if(class(genes) != "character") stop("genes must be a character vector")
  
  cat("Computing pairwise correlations for all genes using tissue gene expression...\n")
  tissue.pairs <- generate.pairs(genes)
  
  cl <- makeCluster(detectCores()-1)
  clusterExport(cl, c("computePairCorrelation", "gene.expr"))
  
  system.time(
    tissue.correlations <- unlist(compute.pairwise.correlations(tissue.pairs, gene.expr, cl))
  )
  # names(hermes.correlations) <- hermes.pairs
  names(tissue.correlations) <- tissue.pairs
  
  stopCluster(cl)
  
  return(correlations.matrix(genes, tissue.correlations))
}

