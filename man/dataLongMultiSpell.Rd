\name{dataLongMultiSpell}
\alias{dataLongMultiSpell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Long Transformation for multi spell analysis
}
\description{
Transform data from short format into long format for discrete multi spell survival analysis and right censoring. 
}
\usage{
dataLongMultiSpell(dataSet, timeColumn, censColumn, 
idColumn, timeAsFactor=FALSE)
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
Character giving the column name of the event status. The event can take multiple values on a discrete scale (0, 1, 2, ...) and repetion of events is allowed. It is assumed that the number zero corresponds to censoring and all number > 0 represent the observed states between transitions. 
}
  \item{idColumn}{
Name of column of identification number of persons as character.
}  
  \item{timeAsFactor}{
Should the time intervals be coded as factor? Default is to use factor. If the argument is false, the column is coded as numeric. 
}
}
\details{
If the data has continuous survival times, the response may be transformed to discrete intervals using function \code{\link{contToDisc}}. The discrete time variable needs to be strictly increasing for each person, because otherwise the order of the events is not distinguishable. Here is an example data structure in short format prior augmentation with three possible states: \\
idColumn=1, 1, ... , 1, 2, 2, ... , n \\
timeColumn= t_ID1_1 < t_ID1_1 < ... < t_ID1_k, t_ID2_1 < t_ID2_2 < ... < t_ID2_k, ... \\
censColumn=0, 1, ... , 2, 1, 0, ... , 0
}
\value{
Original data.frame with three additional columns:
\itemize{
\item {obj: } {Index of persons as integer vector}
\item {timeInt: } {Index of time intervals (formated as factor or integer)}
\item {e0: } {Response in long format as binary vector. Event "e0" is assumed to correspond to censoring. If "e0" is coded one in the in the last observed time interval "timeInt" of a person, then this observation was censored.}
\item {e1: } {Response in long format as binary vector. The event "e1" is the first of the set of possible states "1, 2, 3, ..., X".}
\item {... } {Response in long format as binary vectors. These events correspond to the following states "e2, e3, ...".}
\item {eX } {Response in long format as binary vector. The event "eX" is the last state out of the set of possible states "1, 2, 3, ..., X".}
}
}
\references{
Gerhard Tutz and Matthias Schmid, (2016), \emph{Modeling discrete time-to-event data}, Springer series in statistics, Doi: 10.1007/978-3-319-28158-2

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
\code{\link{contToDisc}}, \code{\link{dataLongTimeDep}}, \code{\link{dataLongCompRisks}}, \code{\link{dataLongCompRisks}}
}
\examples{
############################
# Example of artificial data

# Seed specification
set.seed(-2578)

# Three possible states (0, 1, 2) including censoring
# Discrete time intervals (1, 2, ... , 10)

datFrame <- data.frame(
ID=c(rep(1, 5), rep(2, 3), rep(3, 2), rep(4, 1), rep(5, 3)), 
time=c(c(2, 5, 6, 8, 10), c(1, 6, 7), c(9, 10), c(6), c(2, 3, 4)), 
state=c(c(0, 0, 2, 1, 0), c(1, 2, 2), c(0, 1), c(2), c(0, 2, 1)), 
x=rnorm(n=5+3+2+1+3) )

# Transformation to long format
datFrameLong <- dataLongMultiSpell(dataSet=datFrame, timeColumn="time", 
censColumn="state", idColumn="ID")
head(datFrameLong, 25)

# Fit multi state model without autoregressive terms
library(VGAM)
cRm <- vglm(cbind(e0, e1, e2) ~ timeInt + x, data=datFrameLong, 
family="multinomial")
summary(cRm)
# -> There is no significant effect of x (as expected).

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line