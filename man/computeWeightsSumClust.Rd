% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combineSumClustsDPM.R
\name{computeWeightsSumClust}
\alias{computeWeightsSumClust}
\title{computeWeightsSumClust}
\usage{
computeWeightsSumClust(allocs)
}
\arguments{
\item{allocs}{matrix where each row is a different full set of cluster allocations of all the samples}
}
\value{
PSM and summary clustering obtained from Pitman-Yor process with allocs as input,
weights used to compute the summary PSM from the individual PSMs
}
\description{
Internal function
}
\author{
Magdalena Strauss
}
