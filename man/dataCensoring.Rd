\name{dataCensoring}
\alias{dataCensoring}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Censoring Transformation
}
\description{
Function for transformation of discrete survival times in censoring encoding. Prior this function the data has to be already transformed to long format. With this new generated variable, the discrete censoring process can be analysed instead of the discrete survival process. In discrete survival analysis this information is used to constructs weights for predictive evaluation measures. It is applicable in single event survival analysis.
}
\usage{
dataCensoring(dataSetLong, respColumn, idColumn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSetLong}{
Original data in transformed long format .
}
  \item{respColumn}{
Name of column of discrete survival response as character.
}
  \item{idColumn}{
Name of column of identification number of persons as character.
}
}
\details{
The standard procedure is to use functions such as \code{\link{dataLong}}, \code{\link{dataLongTimeDep}}, \code{\link{dataLongCompRisks}} to augment the data set from short format to long format before using \code{\link{dataCensoring}}. 
}
\value{
Original data set as argument *dataSetLong*, but with added censoring process as first variable in column "yCens"
}
\references{
Ludwig Fahrmeir, (1997), \emph{Discrete failure time models},
LMU Sonderforschungsbereich 386, Paper 91, \url{http://epub.ub.uni-muenchen.de/}

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
\code{\link{contToDisc}}, \code{\link{dataLong}}, \code{\link{dataLongTimeDep}}, \code{\link{dataLongCompRisks}}
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

# Convert to long format based on months
CostMonthsLong <- dataLong (dataSet=CostMonths, timeColumn="timeDisc", censColumn="status")
head(CostMonthsLong, 20)

# Generate censoring process variable
CostMonthsCensor <- dataCensoring (dataSetLong=CostMonthsLong, respColumn="y", idColumn="obj")
head(CostMonthsCensor)
tail(CostMonthsCensor [CostMonthsCensor$obj==1, ], 10)
tail(CostMonthsCensor [CostMonthsCensor$obj==3, ], 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line