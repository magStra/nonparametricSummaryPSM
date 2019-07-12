## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(nonparametricSummaryPSM)

## ----echo=T, results='hide'----------------------------------------------
data("examplePSMs")
npSummaryPSM <- processPSMs(PSMs)

## ----echo=T, results='hide'----------------------------------------------
convExample <- checkConvergence(PSMs)

## ----echo=T--------------------------------------------------------------
library(ggplot2)
dist_Frobenius <- c(convExample$distPY,convExample$distDP)/nrow(convExample$results[[1]]$weightedPSM)^2
method <- c(rep("DPM+PEAR",length(convExample$distPY)),rep("PY+PEAR",length(convExample$distPY)))
number_Chains <- seq(3,24,1)
dist_Frobenius <- data.frame(dist_Frobenius=dist_Frobenius,number_Chains=number_Chains,method=method)
p <- ggplot(dist_Frobenius,aes(x=number_Chains,y=dist_Frobenius)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "Frobenius distance/# matrix elements")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains")+theme(legend.position="top")
p

## ----echo=T--------------------------------------------------------------
cophPY <- convExample$coph[[length(convExample$coph)]]
cophDPM <- convExample$coph_DP[[length(convExample$coph)]]
method <- c(rep("DPM+PEAR",length(cophPY)),rep("PY+PEAR",length(cophPY)))
coph <- c(cophPY,cophDPM)
number_Chains <- seq(2,24,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains")+theme(legend.position="top")
p

## ----echo=T--------------------------------------------------------------
cophPYRel <- cophPY/cophPY[length(cophPY)]
cophDPMRel <- cophDPM/cophDPM[length(cophDPM)]
coph <- c(cophPYRel,cophDPMRel)
number_Chains <- seq(2,24,1)
cophDF <- data.frame(coph=coph,method=method,number_Chains=number_Chains)
p <- ggplot(cophDF,aes(x=number_Chains,y=coph)) + geom_line(aes(color=method)) + geom_point(aes(color=method))+
  scale_y_continuous(name = "ratio of cophenetic correlation")+theme(text=element_text(size=12),axis.text=element_text(size=12))+
  scale_x_continuous(name="number of chains")+theme(legend.position="top")
p

