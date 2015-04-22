\name{lifeTable}
\alias{lifeTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Life Table Construction and Estimates
}
\description{
Constructs a life table and estimates discrete hazard rates, survival functions, cumulative hazard rates and their standard errors without covariates.
}
\usage{
lifeTable(dataSet, timeColumn, censColumn, intervalBorders = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in short format. Must be of class "data.frame".
}
  \item{timeColumn}{
Name of the column with discrete survival times. Must be a scalar character vector.
}
  \item{censColumn}{
Gives the column name of the event indicator (1=observed, 0=censored). Must be of type "character".
}
  \item{intervalBorders}{
Optional names of the intervals for each row, e. g. [a_0, a_1), [a_1, a_2), ..., [a_{q-1}, a_q)
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
List containing an object of class "data.frame" with following columns
\itemize{
		\item{n: } {Number of individuals at risk in a given time interval (integer)}
		\item{events: } {Observed number of events in a given time interval (integer)}
		\item{dropouts: } {Observed number of dropouts in a given time interval (integer)}
		\item{atRisk: } {Estimated number of individuals at risk, corrected by dropouts (numeric)}
		\item{hazard: } {Estimated risk of death (without covariates) in a given time interval}
		\item{seHazard: } {Estimated standard deviation of estimated hazard}
		\item{S: } {Estimated survival curve}
		\item{seS: } {Estimated standard deviation of estimated survival function}
		\item{cumHazard: } {Estimated cumulative hazard function}
		\item{seCumHazard:} {Estimated standard deviation of the estimated cumulative hazard function}	
	}
}
\references{
Jerald F. Lawless, (2000), 
\emph{Statistical Models and Methods for Lifetime Data, 2nd edition},
Wiley series in probability and statistics
%@book{StatModLifeTime,
%author = {Jerald F. Lawless},
%title = {Statistical Models and Methods for Lifetime Data, 2nd edition},
%publisher = {Wiley series in probability and statistics},
%year = {2000}
%}
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}

Matthias Schmid \email{matthias.schmid@imbie.uni-bonn.de}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
# Example with unemployment data
library(Ecdat)
data(UnempDur)

# Extract subset of all persons smaller or equal the median of age
UnempDurSubset <- subset(UnempDur, age<=median(UnempDur$age))
LifeTabUnempDur <- lifeTable (dataSet=UnempDurSubset, timeColumn="spell", censColumn="censor1")
LifeTabUnempDur

# Example with monoclonal gammapothy data
library(survival)
head(mgus)
# Extract subset of mgus
subMgus <- mgus [mgus$futime<=median(mgus$futime), ]
# Transform time in days to intervals [0, 1), [1, 2), [2, 3), ... , [12460, 12461)
mgusInt <- subMgus
mgusInt$futime <- mgusInt$futime + 1
LifeTabGamma <- lifeTable (dataSet=mgusInt, timeColumn="futime", censColumn="death")
head(LifeTabGamma$Output, 25)
plot(x=1:dim(LifeTabGamma$Output)[1], y=LifeTabGamma$Output$hazard, type="l", 
xlab="Time interval", ylab="Hazard", las=1, 
main="Life table estimated marginal hazard rates")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line