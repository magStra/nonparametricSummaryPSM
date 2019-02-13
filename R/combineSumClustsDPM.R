#' Compute summary PSMs (posterior similarity matrices) from a set of multiple PSMs obtained for instance
#' by means of subsampling. This implements the Dirichlet process and Pitman-Yor process based methods
#' for combining PSMs proposed in Strauss et al. Unravelling shared pseudo-trajectories at single-cell resolution.


library(PReMiuM)
library(mcclust)
library(ClusterR)


#' Internal function
#' @title computeSumClustPear
#' @param PSM posterior similarity matrix
#' @return Summary clustering computed using the PEAR criterion (Fritsch and Ickstadt, 2009, using the mcclust package
#'  (Fritsch, 2012))
#' @author Magdalena Strauss
#' @export
computeSumClustPEAR = function(PSM,maxCl=10)
{
  return(norm.label(maxpear(PSM,method="comp",max.k=maxCl)$cl))
}

#' Internal function
#' @title computeWeightsSumClust
#' @param allocs
#' @return PSM and summary clustering obtained from Pitman-Yor process with allocs as input,
#' weights used to compute the summary PSM from the individual PSMs
#' @author Magdalena Strauss
#' @export

computeWeightsSumClust = function(allocs)
{
  colN <- c()
  for (j in 1:dim(allocs)[2])
  {
    colN <- c(colN,sprintf("var%d",j))
  }
  colnames(allocs) <- colN
  allocs <- as.data.frame(allocs)
  output <- c()
  clust <- profRegr(covNames=colN, data=allocs,nSweeps=1000, nBurn=1000, nProgress=500, nFilter=1,
           nClusInit=20,  xModel="Discrete",  sampler="Truncated",  excludeY=T,
           varSelectType="BinaryCluster",dPitmanYor =0.8,alpha=5)

   output$clusObj <- clust
   ng <- dim(allocs)[1]
   dissimObj <- calcDissimilarityMatrix(clust)
   PSM <- matrix(0,ng,ng)
   PSM[lower.tri(PSM, diag=F)] <- 1-dissimObj$disSimMat
   diag(PSM) <- 0.5
   PSM <- t(PSM) + PSM
   output$PSM <- PSM
   clusObj <- calcOptimalClustering(dissimObj)
   output$sumCl <- clusObj
   #get the weights from the output file
   output$rho <- read.table("output_rho.txt",header=F)
   output
}


#' Internal function
#' @title computeWeightsSumClustDPM
#' @param allocs
#' @return PSM and summary clustering obtained from Dirichlet process with allocs as input,
#' weights used to compute the summary PSM from the individual PSMs
#' @author Magdalena Strauss
#' @export
computeWeightsSumClustDPM = function(allocs)
{
  ## This function computes weights for the weighted average of the PSM, using the Dirichlet process method.
  colN <- c()
  for (j in 1:dim(allocs)[2])
  {
    colN <- c(colN,sprintf("var%d",j))
  }
  colnames(allocs) <- colN
  allocs <- as.data.frame(allocs)
  output <- c()
  clust <- profRegr(covNames=colN, data=allocs,nSweeps=1000, nBurn=1000, nProgress=500, nFilter=1,
                    nClusInit=20,  xModel="Discrete",   excludeY=T,
                    varSelectType="BinaryCluster")

  output$clusObj <- clust
  ng <- dim(allocs)[1]
  dissimObj <- calcDissimilarityMatrix(clust)
  PSM <- matrix(0,ng,ng)
  PSM[lower.tri(PSM, diag=F)] <- 1-dissimObj$disSimMat
  diag(PSM) <- 0.5
  PSM <- t(PSM) + PSM
  output$PSM <- PSM
  clusObj <- calcOptimalClustering(dissimObj)
  output$sumCl <- clusObj
  #get the weights from the output file
  output$rho <- read.table("output_rho.txt",header=F)
  output
}

#' @title processPSMs
#' @param PSMs 3-dimensional array of PSMs, for each j PSMs[,,j] is the PSM of subsampled chain j
#' @return weightedPSM: weighted summary PSM obtained using a Pitman-Yor process mixture model
#' with variable selection
#' @return sumClustPEAR: final summary clustering obtained from weightedPSM using the PEAR criterion
#' @return weightedPSM_DP weighted summary PSM obtained using a Dirichlet process mixture model
#' with variable selection
#' @return sumClustPEAR_DP: final summary clustering obtained from weightedPSM_DP using the PEAR criterion
#' @return weights: weights which were used for the computation of the summary PSM (Pitman-Yor based model)
#' @return weights_DP: weights which were used for the computation of the summary PSM (Dirichlet based model)
#' @author Magdalena Strauss
#' @example load("examplePSMs.rda") sumPSM <- processPSMs(PSMS)
#' @export

processPSMs <- function(PSMs)
{
  PSMs[PSMs>1]=1
  PSMs[PSMs<0]=0
  d = dim(PSMs)
  for (k in 1:d[3])
  {diag(PSMs[,,k]) <- 1
  PSMs[,,k] = 0.5*(PSMs[,,k]+t(PSMs[,,k]))}

  sumClusts <- apply(PSMs,3,computeSumClustPEAR)
  clustOb <- computeWeightsSumClust(sumClusts)
  clustObDPM <- computeWeightsSumClustDPM(sumClusts)
  meanWeights <- apply(clustOb$rho,2,mean)
  meanWeightsDPM <- apply(clustObDPM$rho,2,mean)
  #summary PSM obtained from the PY process clustering of the summary clusterings

  weights <-meanWeights/sum(meanWeights)
  weightsDPM <-meanWeightsDPM/sum(meanWeightsDPM)
  PSM_PY <- clustOb$PSM
  PSM_DP <- clustObDPM$PSM
  sumClustPSM <- clustOb$sumCl$clustering
  sumClustPSMDP <- clustObDPM$sumCl$clustering
  weightedPSM <- matrix(0,d[1],d[1])
  for (j in 1:d[3])
  {
    weightedPSM <- weightedPSM + PSMs[,,j]*weights[j]
  }
  weightedPSM_DP <- matrix(0,d[1],d[1])
  for (j in 1:d[3])
  {
    weightedPSM_DP<- weightedPSM_DP + PSMs[,,j]*weightsDPM[j]
  }
  weightedPSM[weightedPSM>1] <- 1
  weightedPSM[weightedPSM<0] <- 0
  diag(weightedPSM) <- 1
  weightedPSM = 0.5*(weightedPSM+t(weightedPSM))
  weightedPSM_DP[weightedPSM_DP>1] <- 1
  weightedPSM_DP[weightedPSM_DP<0] <- 0
  diag(weightedPSM_DP) <- 1
  weightedPSM_DP = 0.5*(weightedPSM_DP+t(weightedPSM_DP))
  #summary clustering using PEAR from the weighted PSM
  sumClustPEAR <- computeSumClustPEAR(weightedPSM)
  sumClustPEAR_DP <- computeSumClustPEAR(weightedPSM_DP)
  tempFiles <- dir(path=getwd(), pattern="output")
  file.remove(tempFiles)
  return(list(weightedPSM=weightedPSM,sumClustPEAR=sumClustPEAR,
                    sumClustPEAR_DP = sumClustPEAR_DP,
                    weights=weights, weightedPSM_DP=weightedPSM_DP,weights_DP=weightsDPM))

}

