% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifyCB.R
\name{classifyCB}
\alias{classifyCB}
\title{Function for classifying cross breed combination based on genomic data (GeneSeek 50K (Illumina) SNP chip).}
\usage{
classifyCB(
  DNA = NULL,
  Predictions = NULL,
  Train = "TrainP+",
  alpha0 = 73.58105,
  Unkn_LogLike = NULL,
  PriorDist = NULL,
  pred_min = 10^(-10)
)
}
\arguments{
\item{DNA:}{matrix (n x m) with 0-1-2 coding of SNP's to be evaluated, typically data returned from ped2num()}

\item{alpha0:}{Deafault value is 73.58105, i.e. the free parameter \eqn{\alpha_0} used
in calculation of \eqn{V(\theta)}.}

\item{Unkn_LogLike:}{the log likelihood for the uniform distribution for unknown breed in PLS-QDA.
Defaults value is NULL, for which the value is calculated based on the m dimensional
space spanned bu the scores in the PLS-model for PB's.}

\item{PriorDist:}{"NULL" (default) for which a totally flat prior is applied,
'informative' for which the informative prior in Vinje et.al.#' is applied,
or a scalar >0 and <1 defining the prior probability for class "Unknown"
and a flat prior is applied to
the rest of the classes,}

\item{pred_min:}{Parameter to be passed to "KL_dist_CBpred_func()". Default to 10^(-10).}
}
\value{
a list bla bla
}
\description{
Function for classifying cross breed combination based on genomic data (GeneSeek 50K (Illumina) SNP chip).
}
\examples{
To come
}
\references{
Vinje H, Brustad HK, Heggli A, Sevillano CA, Van Son M,
Gangsei LE. Classification of breed combinations for slaughter pigs based
on genotypes-modeling DNA samples of crossbreeds as fuzzy sets from purebred
founders. Front Genet. 2023 Dec 4;14:1289130. doi: 10.3389/fgene.2023.1289130.
PMID: 38116292; PMCID: PMC10729766.
}
\author{
Lars Erik Gangsei & Hilde Vinje
}
