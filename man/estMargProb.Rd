\name{estMargProb}
\alias{estMargProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimated Marginal Probabilities
}
\description{
Estimates the marginal probability P(T=t|x) based on estimated hazard rates. The hazard rates may or may not depend on covariates. The covariates have to be equal across all estimated hazard rates. Therefore the given hazard rates should only vary over time.
}
\usage{
estMargProb(haz)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{haz}{
Numeric vector of estimated hazard rates.
}
}
\details{
The argument *haz* must be given for the first q intervals [a_0, a_1), [a_1, a_2), ..., [a_{q-1}, a_q). The last interval is calculated automatically. 
}
\value{
Named vector of estimated marginal probabilities.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{estSurv}}
}
\examples{
# Example unemployment data
library(Ecdat)
data(UnempDur)

# Select subsample
subUnempDur <- UnempDur [1:100, ]

# Convert to long format
UnempLong <- dataLong (dataSet=subUnempDur, timeColumn="spell", censColumn="censor1")
head(UnempLong)

# Estimate binomial model with logit link
Fit <- glm(formula=y ~ timeInt + age + logwage, data=UnempLong, family=binomial())

# Estimate discrete survival function given age, logwage of first person
hazard <- predict(Fit, newdata=subset(UnempLong, obj==1), type="response")

# Estimate marginal probabilities given age, logwage of first person
MarginalProbCondX <- estMargProb (hazard)
MarginalProbCondX
sum(MarginalProbCondX)==1 # TRUE: Marginal probabilities must sum to 1!
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line