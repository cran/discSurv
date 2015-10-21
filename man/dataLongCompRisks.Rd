\name{dataLongCompRisks}
\alias{dataLongCompRisks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Long Competing Risks Transformation
}
\description{
Transforms short data format to long format for discrete survival modelling in the case of competing risks with right censoring. It is assumed that the covariates are not time varying. 
}
\usage{
dataLongCompRisks(dataSet, timeColumn, eventColumns)
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
Character vector giving the column names of the event indicators. It is required that all events are binary and the sum of all event indicators at a given time point, equals one.
}
}
\details{
It is assumed, that only one event happens at a specific time point (competing risks). Either the observation is censored or one of the possible events takes place.
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
Steele Fiona and Goldstein Harvey and Browne William, (2004),
\emph{A general multilevel multistate competing risks model 
for event history data}
Statistical Modelling, volume 4, pages 145-159

Wiji Narendranathan and Mark B. Stewart, (1993),
\emph{Modelling the probability of leaving unemployment: 
competing risks models with flexible base-line hazards},
Applied Statistics, pages 63-83
%@article{narendranathan1993modelling,
%  title={Modelling the probability of leaving unemployment: competing risks models with flexible base-line hazards},
%  author={Narendranathan, Wiji and Stewart, Mark B},
%  journal={Applied Statistics},
%  pages={63--83},
%  year={1993}
%}
%@article{steele2004general,
%  title={A general multilevel multistate competing risks model for event history data, with an application to a study of contraceptive use dynamics},
%  author={Steele, Fiona and Goldstein, Harvey and Browne, William},
%  journal={Statistical Modelling},
%  volume={4},
%  number={2},
%  pages={145--159},
%  year={2004}
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
\code{\link{contToDisc}}, \code{\link{dataLongTimeDep}}
}
\examples{
# Example with unemployment data
library(Ecdat)
data(UnempDur)

# Select subsample
SubUnempDur <- UnempDur [1:100, ]

# Convert competing risk data to long format
SubUnempDurLong <- dataLongCompRisks (dataSet=SubUnempDur, timeColumn="spell", 
eventColumns=c("censor1", "censor2", "censor3", "censor4"))
head(SubUnempDurLong, 20)

# Fit multinomial logit model with VGAM package
# with one coefficient per response
library(VGAM)
multLogitVGM <- vgam(cbind(e0, e1, e2, e3, e4) ~ timeInt + ui + age + logwage,
                    family=multinomial(refLevel=1), 
                    data = SubUnempDurLong)
coef(multLogitVGM)

# Alternative: Use nnet
# Convert response to factor
rawResponseMat <- SubUnempDurLong[, c("e0", "e1", "e2", "e3", "e4")]
NewFactor <- factor(unname(apply(rawResponseMat, 1, function(x) which(x == 1))), 
                    labels = colnames(rawResponseMat))

# Include recoded response in data
SubUnempDurLong <- cbind(SubUnempDurLong, NewResp=NewFactor)

# Construct formula of mlogit model
mlogitFormula <- formula(NewResp ~ timeInt + ui + age + logwage)

# Fit multinomial logit model
# with one coefficient per response
library(nnet)
multLogitNNET <- multinom(formula=mlogitFormula, data=SubUnempDurLong)
coef(multLogitNNET)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line