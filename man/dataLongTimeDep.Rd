\name{dataLongTimeDep}
\alias{dataLongTimeDep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Long Time Dependent Covariates
}
\description{
Transforms short data format to long format for discrete survival modelling of single event analysis with right censoring. Covariates may vary over time. 
}
\usage{
dataLongTimeDep(dataSet, timeColumn, censColumn, idColumn, timeAsFactor=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in short format. Must be of class "data.frame".
}
  \item{timeColumn}{
Character giving the column name of the observed times. It is required that the observed times are discrete (integer).
}
  \item{censColumn}{
Character giving the column name of the event indicator. It is required that this is a binary variable with 1=="event" and 0=="censored".
}
  \item{idColumn}{
Name of column of identification number of persons as character.
}
  \item{timeAsFactor}{
Should the time intervals be coded as factor? Default is to use factor. If the argument is false, the column is coded as numeric. 
}
}
\details{
There may be some intervals, where no additional information on the covariates is observed (e. g. observed values in interval one and three but two is missing). In this case it is assumed, that the values from the last observation stay constant over time until a new measurement was done. 
}
\value{
Original data.frame with three additional columns:
\itemize{
	\item {obj: } {Index of persons as integer vector}
	\item {timeInt: } {Index of time intervals (factor)}
	\item {y: } {Response in long format as binary vector. 1=="event happens in period timeInt" and 0 otherwise}
}
}
\references{
Ludwig Fahrmeir, (1997), \emph{Discrete failure time models},
LMU Sonderforschungsbereich 386, Paper 91, \url{http://epub.ub.uni-muenchen.de/}

W. A. Thompson Jr., (1977), 
\emph{On the Treatment of Grouped Observations in Life Studies},
Biometrics, Vol. 33, No. 3
%@article{DiscSurvFahrmeir,
%author={Ludwig Fahrmeir},
%title={Discrete failure time models}, 
%journal={LMU Sonderforschungsbereich 386, Paper 91, http://epub.ub.uni-muenchen.de/}, 
%year={1997}
%}
%@article{DiscLogitModel,
%author={W. A. Thompson Jr.},
%title={On the Treatment of Grouped Observations in Life Studies}, 
%journal={Biometrics, Vol. 33, No. 3}, 
%year={1977}
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
\code{\link{contToDisc}}, \code{\link{dataLong}}, \code{\link{dataLongCompRisks}}
}
\examples{
# Example Primary Biliary Cirrhosis data
library(survival)
dataSet1 <- pbcseq

# Only event death is of interest
dataSet1$status [dataSet1$status==1] <- 0
dataSet1$status [dataSet1$status==2] <- 1
table(dataSet1$status)

# Convert to months
dataSet1$day <- ceiling(dataSet1$day/30)+1
names(dataSet1) [7] <- "month"

# Convert to long format for time varying effects
pbcseqLong <- dataLongTimeDep (dataSet=dataSet1, timeColumn="month", 
censColumn="status", idColumn="id")
pbcseqLong [pbcseqLong$obj==1, ]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line