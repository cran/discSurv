#' Discrete Survival Analysis
#' 
#' Includes functions for data transformations, estimation, evaluation and
#' simulation of discrete survival analysis. The most important functions are listed below: \itemize{
#' \item\code{\link{contToDisc}}: Discretizes continuous time variable into a
#' specified grid of censored data for discrete survival analysis.
#' \item\code{\link{dataLong}}: Transform data from short format into long
#' format for discrete survival analysis and right censoring.
#' \item\code{\link{dataLongCompRisks}}: Transforms short data format to long
#' format for discrete survival modelling in the case of competing risks with
#' right censoring. \item\code{\link{dataLongTimeDep}}: Transforms short data
#' format to long format for discrete survival modelling of single event
#' analysis with right censoring. \item\code{\link{cIndex}}: Calculates
#' the concordance index for discrete survival models (independent measure of
#' time). \item\code{\link{dataLongSubDist}}: Converts
#' the data to long format suitable for applying discrete subdistribution
#' hazard modelling (competing risks). } 
#' 
#' "DataShort" format is defined as data without repeated measurements. 
#' "DataSemiLong" format consists of repeated measurements, but there are gaps between 
#' the discrete time intervals. "DataLong" format is expanded to include all time 
#' intervals up to the last observation per individual.
#' 
#' \tabular{ll}{ Package: \tab discSurv \cr Type: \tab Package \cr Version:
#' \tab 2.0.0 \cr Date: \tab 2022-03-02 \cr License: \tab GPL-3 \cr }
#' 
#' @name discSurv-package
#' @aliases discSurv-package
#' @docType package
#' @author Thomas Welchowski \email{welchow@@imbie.meb.uni-bonn.de}
#' 
#' Moritz Berger \email{moritz.berger@@imbie.uni-bonn.de} 
#' 
#' David Koehler \email{koehler@@imbie.uni-bonn.de}
#' 
#' Matthias Schmid \email{matthias.schmid@@imbie.uni-bonn.de}
#' 
#' @references
#' \insertRef{bergerTutorial}{discSurv} \cr\cr
#' \insertRef{bergerClassTree}{discSurv} \cr\cr
#' \insertRef{bergerSubdist}{discSurv} \cr\cr
#' \insertRef{schmidDiscMeasure}{discSurv} \cr\cr
#' \insertRef{tutzModelDisc}{discSurv}
#'
#' @keywords package
NULL

# Namespace code
#' @importFrom graphics abline lines plot
#' @importFrom stats aggregate as.formula formula binomial glm glm.control loess model.matrix na.omit pnorm predict qqline qqnorm terms update quantile .checkMFClasses delete.response model.frame na.pass
#' @importFrom utils tail flush.console
#' @importFrom mgcv gam
#' @importFrom functional Curry
#' @importFrom mvtnorm rmvnorm
#' @importFrom data.table rbindlist
#' @importFrom Rdpack reprompt
#' @importFrom VGAM vglm multinomial predictvglm
#' @importFrom geepack geeglm
#' @importFrom treeClust rpart.predict.leaves
#' @importFrom rpart rpart
#' @importFrom ranger ranger
#' @importFrom utils getFromNamespace
#' @importFrom mvnfast rmvn
NULL
