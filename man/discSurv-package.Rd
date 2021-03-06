\name{discSurv-package}
\alias{discSurv-package}
\alias{discSurv}
\docType{package}
\title{
Discrete Survival Analysis
}
\description{
Includes functions for data transformations, estimation, evaluation and simulation of discrete survival analysis. Also discrete life table estimates are available. The most important functions are listed below:
\itemize{
\item{\code{\link{contToDisc}}: }{Discretizes continuous time variable into a specified grid of censored data for discrete survival analysis.}
\item{\code{\link{dataLong}}: }{Transform data from short format into long format for discrete survival analysis and right censoring.}
\item{\code{\link{dataLongCompRisks}}: }{Transforms short data format to long format for discrete survival modelling in the case of competing risks with right censoring.}
\item{\code{\link{dataLongTimeDep}}: }{Transforms short data format to long format for discrete survival modelling of single event analysis with right censoring.}
\item{\code{\link{concorIndex}}: }{Calculates the concordance index for discrete survival models (independent measure of time).}
\item{\code{\link{simCompRisk}}: }{Simulates responses and covariates of discrete competing risk models.}
\item{\code{\link{dataLongSubDist}}: }{Converts the data to long format suitable for applying discrete subdistribution hazard modelling (competing risks).}
}
}
\details{
\tabular{ll}{
Package: \tab discSurv \cr
Type: \tab Package \cr
Version: \tab 1.4.1 \cr
Date: \tab 2019-12-10 \cr
License: \tab GPL-3 \cr
}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de} (Maintainer)

Matthias Schmid \email{matthias.schmid@imbie.uni-bonn.de}
}
\references{
\strong{Main references}:

Moritz Berger, Thomas Welchowski, Steffen Schmitz-Valckenberg and Matthias Schmid, (2019), \emph{A classification tree approach for the modeling of competing risks in discrete time}, Advances in Data Analysis and Classification, vol. 13, issue 4, pages 965-990

Moritz Berger, Matthias Schmid, Thomas Welchowski, Steffen Schmitz-Valckenberg and Jan Beyersmann, (2018), \emph{Subdistribution Hazard Models for Competing Risks in Discrete Time}, Biostatistics, Doi: 10.1093/biostatistics/kxy069

Matthias Schmid, Gerhard Tutz and Thomas Welchowski, (2017), 
\emph{Discrimination Measures for Discrete Time-to-Event Predictions}, Econometrics and Statistics, Elsevier, Doi: 10.1016/j.ecosta.2017.03.008

Gerhard Tutz and Matthias Schmid, (2016), \emph{Modeling discrete time-to-event data}, Springer series in statistics, Doi: 10.1007/978-3-319-28158-2

\strong{Further references}:

Gerhard Tutz, (2012), \emph{Regression for Categorical Data}, Cambridge University Press

Hajime Uno and Tianxi Cai and Lu Tian and L. J. Wei, (2007), 
\emph{Evaluating Prediction Rules for t-Year Survivors With Censored Regression Models},
Journal of the American Statistical Association

Gerds T. A. and M. Schumacher, (2006),
\emph{Consistent estimation of the expected Brier score in general survival models with right-censored event times},
Biometrical Journal, Vol. 48, No. 6, pages 1029-1040

Roger B. Nelsen, (2006), \emph{An introduction to Copulas},
Springer Science+Business Media, Inc.

Patrick J. Heagerty and Yingye Zheng, (2005), 
\emph{Survival Model Predictive Accuracy and ROC Curves}, 
Biometrics 61, 92-105

Steele Fiona and Goldstein Harvey and Browne William, (2004),
\emph{A general multilevel multistate competing risks model 
for event history data}
Statistical Modelling, volume 4, pages 145-159

Jerald F. Lawless, (2000), 
\emph{Statistical Models and Methods for Lifetime Data, 2nd edition},
Wiley series in probability and statistics

Jason P. Fine and Robert J. Gray, (1999), \emph{A proportional hazards model for the subdistribution of a competing risk}, Journal of the American statistical association, Vol. 94, No. 446, pages 496-509

Ludwig Fahrmeir, (1997), \emph{Discrete failure time models},
LMU Sonderforschungsbereich 386, Paper 91, \url{http://epub.ub.uni-muenchen.de/}

Tilmann Gneiting and Adrian E. Raftery, (2007),
\emph{Strictly proper scoring rules, prediction, and
estimation},
Journal of the American Statistical Association 102 (477), 359-376

M. J. van der Laan and J. M. Robins, (2003),
\emph{Unified Methods for Censored Longitudinal Data and Causality},
Springer, New York

Wiji Narendranathan and Mark B. Stewart, (1993),
\emph{Modelling the probability of leaving unemployment: 
competing risks models with flexible base-line hazards},
Applied Statistics, pages 63-83

W. A. Thompson Jr., (1977), 
\emph{On the Treatment of Grouped Observations in Life Studies},
Biometrics, Vol. 33, No. 3

William H. Kruskal, (1958),
\emph{Ordinal Measures of Association},
Journal of the American Statistical Association, Vol. 53, No. 284, pages 814-861

}
\keyword{ package }
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
%%~~ documentation directory ~~
%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%%\examples{
%%~~ simple examples of the most important functions ~~
%%}