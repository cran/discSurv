\name{simCompRisk}
\alias{simCompRisk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of Competing Risks data
}
\description{
Simulates responses and covariates of discrete competing risk models. Responses and covariates are modelled separately by gaussian copulas given kendalls tau. This data can further be analysed with discrete competing risk models. 
}
\usage{
simCompRisk(responseCorr, covariateCorr, sampleSize, covariateSeed = NULL, 
covariateQuantFunc, covariateArgs, intercept = TRUE, trueCoef, responseSeed, 
responseQuantFunc, responseFixArgs, responseLinPredArgs, censorRN, censorArgs, censorSeed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{responseCorr}{
Regulates correlation between event survival times. Numeric Matrix with kendalls tau version b rank correlations. Each cell is restricted to be between -1 and 1. Diagonal elements are always 1.
}
  \item{covariateCorr}{
Regulates correlation between event covariates. Numeric Matrix with kendalls tau version b rank correlations (each cell is restricted to be between -1 and 1). Diagonal elements are always 1. Uses singular, value decomposition for invertation of covariance matrices.
}
  \item{sampleSize}{
Integer scalar specifying the number of rows of the data frame.
}
  \item{covariateSeed}{
Integer scalar, specifying seed of covariate random number generation.
}
  \item{covariateQuantFunc}{
Character vector, giving the marginal quantile functions of all covariates
}
  \item{covariateArgs}{
List of Arguments for each quantile function of the marginal distributions. Each list element should be a named numeric vector (names giving the arguments)
}
  \item{intercept}{
Logical vector, giving if intercept is given in true coefficients for each Event (Default all TRUE)
}
  \item{trueCoef}{
List of numeric vectors containing the beta of true coefficients, e. g. linear predictor eta = X %*% beta; with X as simulated design matrix. Each list element has all beta coefficients of one event
}
  \item{responseSeed}{
Integer scalar, specifying seed of response random number generation
}
  \item{responseQuantFunc}{
Character vector, giving the marginal quantile functions of all survival
}
  \item{responseFixArgs}{
List of Arguments for each quantile function of the marginal distributions. Each list element should be a named numeric vector.
}
  \item{responseLinPredArgs}{
List of lists, specifying the relationship of linear predictor and parameters of the marginal distributions. Each list element is a list of all functional relationships between linear predictor and parameters of one marginal distribution. Each list element is a function giving the functional relationship between linear predictor and one parameter.
}
  \item{censorRN}{
Integer scalar, specifying seed of censoring random number generation
}
  \item{censorArgs}{
Named numeric vector, giving all arguments of the marginal censoring distribution
}
  \item{censorSeed}{
Integer scalar, specifying seed of censoring random number generation
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
List with following components
\itemize{
	\item {Data: } {Simulated data of class "data.frame"}
	\item {covariateCorr: } {Original input rank correlation structure of covariates}
	\item {samleSize: } {Sample size}
  \item {covariateSeed: } {Covariate seed}
	\item {...} {(all arguments specified in Input are saved other the corresponding names)}
}
}
\references{
Wiji Narendranathan and Mark B. Stewart, (1993),
\emph{Modelling the probability of leaving unemployment: 
competing risks models with flexible base-line hazards},
Applied Statistics, pages 63-83

Roger B. Nelsen, (2006), 
\emph{An introduction to Copulas},
Springer Science+Business Media, Inc.

Steele Fiona and Goldstein Harvey and Browne William, (2004),
\emph{A general multilevel multistate competing risks model 
for event history data}
Statistical Modelling, volume 4, pages 145-159
%@article{narendranathan1993modelling,
%  title={Modelling the probability of leaving unemployment: 
%  competing risks models with flexible base-line hazards},
%  author={Narendranathan, Wiji and Stewart, Mark B},
%  journal={Applied Statistics},
%  pages={63--83},
%  year={1993}
%}
%@article{steele2004general,
%  title={A general multilevel multistate competing risks model for event history data, 
%  with an application to a study of contraceptive use dynamics},
%  author={Steele, Fiona and Goldstein, Harvey and Browne, William},
%  journal={Statistical Modelling},
%  volume={4},
%  number={2},
%  pages={145--159},
%  year={2004}
%}
%@book{IntroCop,
%title={An introduction to Copulas},
%author={Roger B. Nelsen},
%publisher={Springer Science+Business Media, Inc.},
%year={2006}
%}
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dataLongCompRisks}}, \code{\link[VGAM:vgam]{vgam}}
}
\examples{
library(Matrix)
library(matrixcalc)

########################
# Design of Simulation

# Monte carlo samples for each design = 10
# SampleSize = 100 in each replicate

# Responses: 
# R1 i.n.d. ~ Weibull (shape=0.5, scale=exp(eta)) 
# E(R1) = gamma (3)*exp(eta)
# R2 i.n.d. ~ Weibull (shape=1, scale=exp(eta))
# E(R2) = gamma (2)*exp(eta)
# R3 i.n.d. ~ Weibull (shape=2, scale=exp(eta))
# E(R3) = gamma (1.5)*exp(eta)

# Random Censoring
# Z = 2/3 * (gamma (3) + gamma (2) + gamma (1.5)) = 2.590818
# Independent of survival times C_i i.i.d ~ Exp (lambda=Z)

# Correlation structure of survival times:
# Specify with kendalls tau -> spearmans rho
# Kendalls tau matrix:
# R1   R2   R3
# R1  1   0.3  0.4
# R2  0.3   1  0.5
# R3  0.4 0.5    1

# Covariates: 
# V1: Binary variable ~ Bin (n=2, pi=0.5) -> E (V1) = 1
# V2: Continuous positive variable ~ Gamma (shape=1, scale=1) -> E (V2) = 1
# V3: Continuous variable ~ Normal (mu=0, sigma^2=1) -> E (V5) = 0

# True linear predictor:
# eta = X %*% beta
# beta_1 = c(-0.5, 1, 0.5)
# beta_2 = c(1, -1, 1)
# beta_3 = c(-0.5, -0.5, 2)

# Correlation Structure in Covariates:
# Specify with kendalls tau -> pearsons rho
# V1    V2     V3
# V1       1 -0.25      0
# V2   -0.25     1   0.25
# V3       0  0.25      1

# Design Correlation Matrix of covariates
DesignCovariateCor <- diag(3)
DesignCovariateCor [lower.tri(DesignCovariateCor)] <- c(-0.25, 0, 0.25)
DesignCovariateCor [upper.tri(DesignCovariateCor)] <- c(-0.25, 0, 0.25)

# Check if correlation matrix is symmetric
sum(DesignCovariateCor-t(DesignCovariateCor))==0 # TRUE -> ok
# Check if correlation matrix is positive definite
is.positive.definite (DesignCovariateCor) # TRUE -> ok
# Check if correlation matrix is transformed pearson matrix is positive, definite
is.positive.definite (apply(DesignCovariateCor, c(1,2), tauToPearson)) # TRUE -> ok

# Design Correlation Matrix positive definite after transformation
DesignResponseCor <- diag(3)
DesignResponseCor [lower.tri(DesignResponseCor)] <- c(0.3, 0.4, 0.5)
DesignResponseCor [upper.tri(DesignResponseCor)] <- c(0.3, 0.4, 0.5) 

# Check if response correlation matrix is symmetric
sum(DesignResponseCor-t(DesignResponseCor))==0
# Check if response correlation matrix is positive definite
is.positive.definite (DesignResponseCor)
# Check if response correlation matrix is transformed pearson matrix is positive, definite
is.positive.definite (apply(DesignResponseCor, c(1,2), tauToPearson))

# Simulate raw data
DesignSampleSize <- 100
MonteCarloSamples <- 10
SimOutput <- vector("list", MonteCarloSamples)
for(j in 1:MonteCarloSamples) {
  SimOutput [[j]] <- simCompRisk (responseCorr=DesignResponseCor, covariateCorr=DesignCovariateCor, 
	  covariateSeed=NULL, sampleSize=100, covariateQuantFunc=c("qbinom", "qgamma", "qnorm"), 
	  covariateArgs=list(c(size=2, prob=0.5), c(shape=1, scale=1), c(mean=0, sd=1)), 
	  intercept=c(FALSE, FALSE, FALSE), 
	  trueCoef=list(c(-0.5, 1, 0.5), c(1, -1, 1), c(-0.5, -0.5, 2)), 
	  responseSeed=NULL, responseQuantFunc=c("qweibull", "qweibull", "qweibull"), 
	  responseFixArgs=list(c(shape=0.5), c(shape=1), c(shape=2)), 
	  responseLinPredArgs=list(list(scale=function (x) {exp(x)}), 
	  list(scale=function (x) {exp(x)}), 
	  list(scale=function (x) {exp(x)})), censorRN="rexp", 
	  censorArgs=c(rate=2/3 * (gamma (3) + gamma (2) + gamma (1.5))), censorSeed=NULL)
}
SimOutput [[1]]
SimOutput [[5]]

OverviewData <- data.frame(Min=rep(NA, MonteCarloSamples), 
Max=rep(NA, MonteCarloSamples), Censrate=rep(NA, MonteCarloSamples))
for(j in 1:MonteCarloSamples) {
  
  # Calculate Range of observed time
  OverviewData [j, 1:2] <- range (SimOutput [[j]]$Data [, "Time"])
  
  # Calculate censoring rate
  OverviewData [j, "Censrate"] <- mean(SimOutput [[j]]$Data [, "Censor"])
  
}
OverviewData
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line