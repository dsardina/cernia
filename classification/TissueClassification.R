

# library(party)

cerna.classification.scores <- function(resz){
  tmp <- do.call(rbind, resz)
  c.pairs <- colMeans(tmp)
  
  u.cerna <- as.data.frame(do.call(rbind, strsplit(names(c.pairs), split = "-")))
  colnames(u.cerna) <- c("ceRNA1", "ceRNA2")
  u.cerna$Score <- c.pairs
  
  return(u.cerna)
}

classify <- function(new.data, svms){
  # Create a vector of result for the total amount of svm's batches
  resz <- list("vector", length(svms))
  
  for (n.r in 1:length(svms)) {
    cat("Batch", n.r ,"...\n")
    conv.factor <- function (f) {
      lvls <- as.numeric(levels(f))
      f1 <- as.numeric(f)
      for (i in 1:length(lvls)) {
        f1[f1 == i] <- lvls[i]
      }
      names(f1) <- names(f)
      return (f1)
    }
    
    svm.batch <- svms[[n.r]]
    l     <- length(svm.batch)
    
    results <- conv.factor(predict(svm.batch[[1]], new.data))
    for (i in 2:l) {
      cat("TS",i,"...")
      results <- results + conv.factor(predict(svm.batch[[i]], new.data))
    }
    
    resz[[n.r]]  <- results / l
  }
  
  return(resz)
}