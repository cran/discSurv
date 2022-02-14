\name{dataLongCompRisksTimeDep}
\alias{dataLongCompRisksTimeDep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Long Competing Risks Time Dependent Covariates Transformation
}
\description{
Transforms short data format to long format for discrete survival modelling in the case of competing risks with right censoring. Covariates may vary over time. 
}
\usage{
dataLongCompRisksTimeDep(dataSet, timeColumn, eventColumns, 
eventColumnsAsFactor=FALSE, idColumn, timeAsFactor=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in short format. Must be of class "data.frame".
}
  \item{timeColumn}{
Character giving the column name of the observed times. It is required that the observed times are discrete (integer).
}
  \item{eventColumns}{
Character vector giving the column names of the event indicators (excluding censoring column). It is required that all events are binary encoded. If the sum of all event indicators is zero, then this is interpreted as a censored observation. Alternatively a column name of a factor representing competing events can be given. In this case the argument "eventColumnsAsFactor" has to be set TRUE and the first level is assumed to represent censoring. 
}
  \item{eventColumnsAsFactor}{
Should the argument eventColumns be intepreted as column name of a factor variable(logical scalar)? Default is FALSE. 
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
Original data set in long format with additional columns
\itemize{
	\item {obj: } {Gives identification number of objects (row index in short format) (integer)}
	\item {timeInt: } {Gives number of discrete time intervals (factor)}
	\item {responses: } {Columns with dimension count of events + 1 (censoring)
		\itemize{ 
			\item {e0: } {No event (observation censored in specific interval)}
			\item {e1: } {Indicator of first event, 1 if event takes place and 0 otherwise}
			\item {...} {...}
			\item {ek: } {Indicator of last k-th event, 1 if event takes place and 0 otherwise }
		}
	}
}
}
\references{
Ludwig Fahrmeir, (1997), \emph{Discrete failure time models},
LMU Sonderforschungsbereich 386, Paper 91, \url{https://epub.ub.uni-muenchen.de/}

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

# Convert to months
pbcseq$day <- ceiling(pbcseq$day/30)+1
names(pbcseq) [7] <- "month"
pbcseq$status <- factor(pbcseq$status)

# Convert to long format for time varying effects
pbcseqLong <- dataLongCompRisksTimeDep(dataSet=pbcseq, timeColumn="month", 
eventColumns="status", eventColumnsAsFactor=TRUE, idColumn="id", 
timeAsFactor=TRUE)
head(pbcseqLong)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line