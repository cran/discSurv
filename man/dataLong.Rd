\name{dataLong}
\alias{dataLong}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Long Transformation
}
\description{
Transform data from short format into long format for discrete survival analysis and right censoring. Data is assumed to include no time varying covariates, e. g. no follow up visits are allowed. It is assumed that the covariates stay constant over time, in which no information is available.
}
\usage{
dataLong(dataSet, timeColumn, censColumn, timeAsFactor=TRUE, 
remLastInt=FALSE, aggTimeFormat=FALSE, lastTheoInt=NULL)
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
  \item{timeAsFactor}{
Should the time intervals be coded as factor? Default is to use factor. If the argument is false, the column is coded as numeric. 
}
 \item{remLastInt}{
Should the last theoretical interval be removed in long format? Default is no deletion. This is only important, if the short format data includes the last theoretic interval [a_q, Inf). There are only events in the last theoretic interval, so the hazard is always one and these observations have to be excluded for estimation. 
}
 \item{aggTimeFormat}{
Instead of the usual long format, should every obseration have all time intervals? (logical scalar) Default is standard long format. In the case of nonlinear risk score models, the time effect has to be integrated out before these can be applied to the C-index \code{\link{concorIndex}}. 
}
\item{lastTheoInt}{
Gives the number of the last theoretic interval (integer scalar). Only used, if aggTimeFormat==TRUE.
}
}
\details{
If the data has continuous survival times, the response may be transformed to discrete intervals using function \code{\link{contToDisc}}. If the data set has time varying covariates the function \code{\link{dataLongTimeDep}} should be used instead. In the case of competing risks and no time varying covariates see function \code{\link{dataLongCompRisks}}.
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

Matthias Schmid \email{matthias.schmid@imbie.uni-bonn.de}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{contToDisc}}, \code{\link{dataLongTimeDep}}, \code{\link{dataLongCompRisks}}
}
\examples{
# Example unemployment data
library(Ecdat)
data(UnempDur)

# Select subsample
subUnempDur <- UnempDur [1:100, ]
head(subUnempDur)

# Convert to long format
UnempLong <- dataLong (dataSet=subUnempDur, timeColumn="spell", censColumn="censor1")
head(UnempLong, 20)

# Is there exactly one observed event of y for each person?
splitUnempLong <- split(UnempLong, UnempLong$obj)
all(sapply(splitUnempLong, function (x) sum(x$y))==subUnempDur$censor1) # TRUE

# Second example: Acute Myelogenous Leukemia survival data
library(survival)
head(leukemia)
leukLong <- dataLong (dataSet=leukemia, timeColumn="time", censColumn="status")
head(leukLong, 30)

# Estimate discrete survival model
estGlm <- glm(formula=y ~ timeInt + x, data=leukLong, family=binomial())
summary(estGlm)

# Estimate survival curves for non-maintained chemotherapy
newDataNonMaintained <- data.frame(timeInt=factor(1:161), x=rep("Nonmaintained"))
predHazNonMain <- predict(estGlm, newdata=newDataNonMaintained, type="response")
predSurvNonMain <- cumprod(1-predHazNonMain)

# Estimate survival curves for maintained chemotherapy
newDataMaintained <- data.frame(timeInt=factor(1:161), x=rep("Maintained"))
predHazMain <- predict(estGlm, newdata=newDataMaintained, type="response")
predSurvMain <- cumprod(1-predHazMain)

# Compare survival curves
plot(x=1:50, y=predSurvMain [1:50], xlab="Time", ylab="S(t)", las=1, 
type="l", main="Effect of maintained chemotherapy on survival of leukemia patients")
lines(x=1:161, y=predSurvNonMain, col="red")
legend("topright", legend=c("Maintained chemotherapy", "Non-maintained chemotherapy"), 
col=c("black", "red"), lty=rep(1, 2))
# The maintained therapy has clearly a positive effect on survival over the time range
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line