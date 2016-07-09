
prepare <- function(genes, correlations, recommendations, scores){
  
  cat("Filtering and preparing...\n")
  scores.1 <- lapply(genes, function (g) {
    
    flt <- intersect(setdiff(names(scores[[g]]), g), genes)
    sco <- scores[[g]][flt]
    
    if (length(sco) == 0) {
      return (matrix(NA, nrow = 0, ncol = 8))
    }
    
    sco <- t(sapply(sco, function(x) { x }))
    
    if (nrow(sco) == 0) {
      return (matrix(NA, nrow = 0, ncol = 8))
    }
    rownames(sco) <- paste0(g,"-",rownames(sco))
    
    sco[is.infinite(sco)] <- NA
    sco[,1]               <- log(sco[,1])
    sco[,7]               <- log(recommendations[g, flt])
    corr                  <- correlations[g, flt]
    
    if (!all(!is.na(corr))) {
      corr.erro <<- corr
      stop("ERRORE")
    }
    
    sco <- cbind(sco, corr)
    return (sco)
  })
  
  cat("Merging...\n")
  
  scores.matrix <- do.call(rbind, scores.1)
  colnames(scores.matrix) <- c("P", "MUT1", "MUT2", "MUT3", "MUT4", "MUT5", "PRJ", "CORR")
  
  scores.matrix.noNA <- na.omit(scores.matrix)
  
  return(list(scores.matrix=scores.matrix, scores.matrix.noNA=scores.matrix.noNA))
}