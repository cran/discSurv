\name{dataCensoringShort}
\alias{dataCensoringShort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Censoring Transformation for short formats
}
\description{
Function for transformation of discrete survival times in censoring encoding. In contrast to \code{\link{dataCensoring}} this function needs the original data in short format as argument. With the new generated variable "yCens", the discrete censoring process can be analyzed instead of the discrete survival process. In discrete survival analysis this information is used to constructs weights for predictive evaluation measures. It is applicable in single event survival analysis.
}
\usage{
dataCensoringShort(dataSet, eventColumns, timeColumn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in transformed long format.
}
  \item{eventColumns}{
Name of columns of event columns (character vector). The event columns have to be in binary format. If the sum of all events equals zero in a row, then this observation is interpreted as censored. 
}
  \item{timeColumn}{
Name of column with discrete time intervals (character scalar).
}
}
%\details{
%The standard procedure is to use functions such as \code{\link{dataLong}}, \code{\link{dataLongTimeDep}}, \code{\link{dataLongCompRisks}} to augment the data set from short format to long format before using \code{\link{dataCensoringShort}}. 
%}
\value{
Original data set as argument *dataSet*, but with added censoring process as first variable in column "yCens".
}
\references{
Gerhard Tutz and Matthias Schmid, (2016), \emph{Modeling discrete time-to-event data}, Springer series in statistics, Doi: 10.1007/978-3-319-28158-2

Ludwig Fahrmeir, (1997), \emph{Discrete failure time models},
LMU Sonderforschungsbereich 386, Paper 91, \url{https://epub.ub.uni-muenchen.de/}

W. A. Thompson Jr., (1977), 
\emph{On the Treatment of Grouped Observations in Life Studies},
Biometrics, Vol. 33, No. 3
%@article {DiscSurvFahrmeir,
%author={Ludwig Fahrmeir},
%title={Discrete failure time models}, 
%journal={LMU Sonderforschungsbereich 386, Paper 91, http://epub.ub.uni-muenchen.de/}, 
%year={1997}
%}
%@article {DiscLogitModel,
%author={W. A. Thompson Jr.},
%title={On the Treatment of Grouped Observations in Life Studies}, 
%journal={Biometrics, Vol. 33, No. 3}, 
%year={1977}
%}
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dataCensoring}}, \code{\link{contToDisc}}, \code{\link{dataLong}}, \code{\link{dataLongTimeDep}}, \code{\link{dataLongCompRisks}}
}
\examples{
library(pec)
data(cost)
head(cost)
IntBorders <- 1:ceiling(max(cost$time)/30)*30
subCost <- cost [1:100, ]

# Convert from days to months
CostMonths <- contToDisc (dataSet=subCost, timeColumn="time", intervalLimits=IntBorders)
head(CostMonths)

# Generate censoring process variable in short format
CostMonthsCensorShort <- dataCensoringShort (dataSet=CostMonths, 
eventColumns="status", timeColumn="time")
head(CostMonthsCensorShort)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line