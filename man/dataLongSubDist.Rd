\name{dataLongSubDist}
\alias{dataLongSubDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data Matrix and Weights for Discrete Subdistribution Hazard Models
}
\description{
Generates the augmented data matrix and the weights required for discrete subdistribution hazard modeling with right censoring.
}
\usage{
dataLongSubDist(dataSet, timeColumn, eventColumns, eventFocus,
timeAsFactor=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataSet}{
Original data in short format. Must be of class "data.frame".
}
  \item{timeColumn}{
Character specifying the column name of the observed event times. It is required that the observed times are discrete (integer).
}
  \item{eventColumns}{
Character vector specifying the column names of the event indicators (excluding censoring events). It is required that a 0-1 coding is used for all events. If the sum of the event indicators is zero this will be interpreted as a censoring event. 
}
  \item{eventFocus}{
Column name of the event of interest (type 1 event).
}
  \item{timeAsFactor}{
Logical indicating whether time should be coded as a factor in the augmented data matrix. If FALSE, a numeric coding will be used. 
}
}
\details{
This function sets up the augmented data matrix and the weights that are needed for weighted maximum likelihood (ML) estimation of the discrete subdistribution model proposed by Berger et al. (2018). The model is a discrete-time extension of the original subdistribution model proposed by Fine and Gray (1999).
}
\value{
Data frame with additional column "subDistWeights". The latter column contains the weights that are needed for fitting a weighted binary regression model, as described in Berger et al. (2018). The weights are calculated by a life table estimator for the censoring event.
}
\references{
Moritz Berger, Matthias Schmid, Thomas Welchowski, Steffen Schmitz-Valckenberg and Jan Beyersmann (2018). \emph{Subdistribution Hazard Models for Competing Risks in Discrete Time}. Submitted manuscript.

Jason P. Fine and Robert J. Gray (1999). \emph{A proportional hazards model for the subdistribution of a competing risk}. Journal of the American Statistical Association 94, pp. 496-509.
}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dataLong}}
}
\examples{
# Example with unemployment data
library(Ecdat)
data(UnempDur)

# Generate subsample, reduce number of intervals to k = 5
SubUnempDur <- UnempDur [1:500, ]
SubUnempDur$time <- as.numeric(cut(SubUnempDur$spell, c(0,4,8,16,28)))

# Convert competing risks data to long format
# The event of interest is re-employment at full job
SubUnempDurLong <- dataLongSubDist (dataSet=SubUnempDur, timeColumn="time", 
eventColumns=c("censor1", "censor2", "censor3"), eventFocus="censor1")
head(SubUnempDurLong)

# Fit discrete subdistribution hazard model with logistic link function
logisticSubDistr <- glm(y ~ timeInt + ui + age + logwage,
                    family=binomial(), data = SubUnempDurLong, 
                    weights=SubUnempDurLong$subDistWeights)
summary(logisticSubDistr)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line