library(ROCR)

# compute the optimal point
opt.cut <- function(perf, pred){
  cut.ind <- mapply(FUN=function(x, y, p){
    d <- (x - 0)^2 + (y-1)^2
    ind <- which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

# Function that plot the curves
plot.rocs <- function(type, save.plot = TRUE){
  # Compute predictions
  preds <- prediction(resz, dts)
  
  # Compute performances
  roc.perf <- performance(preds, measure = "tpr", x.measure = "fpr")
  
  # Plots
  if(save.plot){
    size <- 8*600
    jpeg(file = paste0(type,".roc.curves.jpg"), width = size, height = size, res=600)
    
    plot(roc.perf, col = "blue", avg = "vertical", spread.estimate = "stderror", lwd = 4, main=paste(toupper(type),"- ROC Curves"), cex.lab=2, cex.main = 2.7)
    plot(roc.perf, col = "red", lty = 2, lwd = 1.8, add = TRUE)
    
    dev.off()
  }
  
  # Compute the AUC
  auc.perf <- performance(preds, measure = "auc")
  cat(type, " AUC values: ", unlist(auc.perf@y.values))
  
  # Compute the optimal cut
  cat(type, " Optimal cuts:")
  print(opt.cut(roc.perf, preds))

  return(list(performance = unlist(auc.perf@y.values), opt.cut = opt.cut(roc.perf, preds)))
  #   save(list = c("preds", "roc.perf", "pr.perf", "acc.perf", "auc.perf"), file = paste0(type,".rocr.results.RData"))
}
