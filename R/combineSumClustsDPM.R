#' Compute summary PSMs (posterior similarity matrices) from a set of multiple PSMs obtained for instance
#' by means of subsampling. This implements the Dirichlet process and Pitman-Yor process based methods
#' for combining PSMs proposed in Strauss et al. Unravelling shared pseudo-trajectories at single-cell resolution.


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


#' @title checkConvergence
#' @param PSMs 3-dimensional array of PSMs, for each j PSMs[,,j] is the PSM of subsampled chain j, j = 1, ..., m
#' @return results: list of length m-1 of results of processPSM function applied to the following subsets of PSMs: 1,2; 1,2,3; 1,2,3,4; ... 1,2,...,m
#' @return coph: vector of length m-1 of cophenetic correlation coefficients measuring how well the summary PSM obtained using the PY+PEAR method on the first
#' k submatrices (k=1, ..., l) is described by a hierarchcial clustering tree for a summary PSM based on the first l subsampled chains
#' @return coph_DP: as coph, but for summary PSMs obtained by the DPM+PEAR method
#' @return distPY: vector of length m-2 of Frobenius (Euclidean) norm of distances between consecutive summary PSMs, that is summary PSMs obtained from chains 1,2,..,k,k+1 and 1,2,...k (PY+PEAR)
#' @return distDP: as distPY, but for DPM+PEAR method
#' @author Magdalena Strauss
#' @export
checkConvergence <- function(PSMs)
{
  #use processPMs function to compute summary PSMs for chains 1:M, where M = 2, ..., number of chains
  #then compute cophenetic correlation for all of these, to see how it changes across the number of PSMs
  #also includes test for testing if PSMs are the same
  L = dim(PSMs)[3]
  results <- list()
  coph <- list()
  coph_DP <- list()
  distPY <- c()
  distDP <- c()
  for (l in 2:L)
  {
    PSMs1 = PSMs[,,1:l]
    results[[l-1]] <- processPSMs(PSMs1)
    MPY <- results[[l-1]]$weightedPSM
    MDP <- results[[l-1]]$weightedPSM_DP
  }
  for (l in 2:(L-1))
  {
    distPY <- c(distPY,norm(results[[l]]$weightedPSM-results[[l-1]]$weightedPSM,"f"))
    distDP <- c(distDP,norm(results[[l]]$weightedPSM_DP-results[[l-1]]$weightedPSM_DP,"f"))
    hClust <- hclust(dist(1-results[[l]]$weightedPSM),method="average")
    cophHClust <- cophenetic(hClust)
    hClustDP <- hclust(dist(1-results[[l]]$weightedPSM_DP),method="average")
    cophHClustDP <- cophenetic(hClustDP)
    aa <- rep(0,l)
    bb <- rep(0,l)
    for (k in 1:l)
    {
      aa[k] <- cor(dist(1-results[[k]]$weightedPSM),cophHClust)
      bb[k] <- cor(dist(1-results[[k]]$weightedPSM_DP),cophHClustDP)
    }
    coph[[l-1]] <- aa
    coph_DP[[l-1]] <- bb
  }
  return(list(results=results,coph=coph,coph_DP=coph_DP,distPY=distPY,distDP=distDP))
}

#' @title RhatConcentration
#' @param concentrationSamples: nIterations x nChains matrix of nIterations samples of the concentration
#' parameter alpha for each of nChains subsampled chains
#' @param nChainsTest: number of subsamples for which we want to test if the number of sufficient for convergence
#' @return GR-statistics across subsampled chains (see Strauss et al. 2019) for groups of
#' @author Magdalena Strauss
#' @export
RhatConcentration <- function(concentrationSamples,nChainsTest)
{
  nChains <- ncol(concentrationSamples)
  nn <- floor(nrow(concentrationSamples)/2)
  concentrationSamples <- concentrationSamples[-(1:nn),]
  a <- floor(nChains/nChainsTest)
  output <- c()
  for (k in 1:10){
    yy <- 1:nChains
    data <- c()
    for (l in 1:a){
      xx <- sample(yy,nChainsTest,replace = F)
      yy <- setdiff(yy,xx)
      x <- as.vector(concentrationSamples[,xx])
      names(x) <- NULL
      data <- rbind(data,x)
    }
    x1MC = as.mcmc.list(lapply(as.data.frame(t(data)), mcmc))
    XX1 = gelman.diag(x1MC,transform = F,autoburnin = F)
    output <- append(output,XX1[[1]][[1]])
  }
  return(output)
}
