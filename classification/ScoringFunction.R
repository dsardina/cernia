cat("Initialize data for vector of scores computation...\n")


#############
# Load MREs 
#############

cat("Loading MREs...\n")
load("data/mres.RData")

if(!exists("mres") || ncol(mres) != 5){
  stop("Invalid data.frame for mres, check the example below:\n
    Mirna\tTarget\tenergy\tgap_l\tgap_r
    MIMAT0000062\t23208\t-26.01\t202\t222
    MIMAT0000063\t29128\t-26.62\t494\t515")
}


######################
##  Utility functions 
######################

# Return the miRs interacting with given target
getInteractingMirs <- function(target) return(interactions[interactions[,2] == target, 1])

# Search for the MREs between target and miR
getMres    <- function(target, miR) mres[mres[,2] == target & mres[,1] == miR, ]
getAllMres <- function(t1, t2, cm) mres[mres[,2] %in% c(t1,t2) & mres[,1] %in% cm, ]

# Calculate the maximum distance from leftmost and rightmost MRE
getMaxDist <- function(positions) (abs(max(positions[,1]) - min(positions[,2])))

# Compute score 1
getScore1 <- function(commonMiRs, miRs) {
  log(length(commonMiRs)/length(miRs))
}

# Compute score 2
getScore2 <- function(cm, shMre){
  # compute the MREs' weighted density for every miR
  return (sum(sapply(cm, function (miR) {
    mres <- shMre[shMre[,1] == miR,]
    if (nrow(mres) <= 0) return (1)
    d <- getMaxDist(mres[,c("gap_l", "gap_r")])
    return (
      log(nrow(mres) / d)
    )
  })))
}

# Compute score 2 bis
getScore2bis <- function(cm, shMre){
  # compute the MREs' weighted density for every miR
  return (sum(sapply(cm, function (miR) {
    mres <- shMre[shMre[,1] == miR,]
    if (nrow(mres) <= 0) return (1)
    d <- getMaxDist(mres[,c("gap_l", "gap_r")])
    return (
      log(sum(abs(mres[,3])) / d)
    )
  })))
}

# Compute score 3
getScore3 <- function(cm, shMre){
  return (sum(log(sapply(cm, function (miR) {
    # mmres <- shMre[,1] == miR # get only the mre for the current miR
    positions <- shMre[shMre[,1] == miR, c("gap_l", "gap_r")] # get positions of the MREs
    if (nrow(positions) <= 0) return (1)
    return ((getMaxDist(positions))^2/sum((positions[,2]-positions[,1])^2))
  }))))
}

# Compute score 4
getScore4 <- function(shMre){
  b <- nrow(shMre)
  if (b == 0) return (1)
  # mjmre <- length(unique(shMre[,1])) # number of different miRs referred to MREs found
  return (log((b - length(unique(shMre[,1])) + 1)/b))
}

# Compute the scores for the current pair
inner <- function (j, i, mi, zz){
  mj    <- getInteractingMirs(j)
  cat(paste("Working on S[",i,",",j,"]\n"),file=zz)
  if (length(mi) == 0 || length(mj) == 0) return (-Inf)
  cm    <- intersect(mi, mj)
  if (length(cm) <= 0) { 
    return (c(1, -Inf, -Inf, -Inf, -Inf, -Inf, cerna.recommendations[i,j]))
  }
  ft <- fisher.test(matrix(c(length(mi),length(cm),length(cm),length(mj)),2,2))
  shMre <- getAllMres(i, j, cm)
  return (c(ft$p.value, getScore1(cm,mi), getScore2(cm, shMre),
            getScore2bis(cm, shMre), getScore3(cm, shMre), 
            getScore4(shMre), cerna.recommendations[i,j]))
}

# Call the inner function for the current target
outer <- function(i, cl){
  l   <- length(targets)
  idx <- which(targets == i)
  
  mi  <- getInteractingMirs(i)
  zz <- file("all.Rout.txt", open = "at")
  
  r <- lapply(targets[(idx+1):l], inner, i, mi, zz)
  names(r) <- targets[(idx+1):l]
  close(zz)
  
  return (r)
}

# Execute the scoring functions for each pair
execute <- function (cl) {
  tars <- targets[-length(targets)]
  r <- parLapply(cl=cl,X=tars,fun=outer, NA)
  names(r) <- tars
  return (r)
}


#############################################################
##  Load targets from the dataset and compute scoring vector
#############################################################

targets <- genes

cat("Init parallel execution...\n")

# Init parallel execution
cl <- makeCluster(detectCores()-1)

clusterExport(cl,c(
  "mres", "targets", "interactions", "cerna.recommendations",
  "getAllMres", "getMaxDist", "getInteractingMirs", "inner",
  "getScore1", "getScore2", "getScore2bis", "getScore3", "getScore4"
))

cat("Computing vector of scores for all pairs...\n")
cernia.scores <- execute(cl)

stopCluster(cl)
rm(cl)

source("classification/ComputeCorrelations.R")
source("classification/Prepare.R")

correlations <- compute.correlations(targets, gene.expr)
final.scores <- prepare(targets[-length(targets)], correlations, cerna.recommendations, cernia.scores)

cat("Done!\n")