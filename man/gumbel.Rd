\name{gumbel}
\alias{gumbel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gumbel Link Function
}
\description{
Constructs the link function with gumbel distribution in approriate format for use in generalized, linear models
}
\usage{
gumbel()
}
%- maybe also 'usage' for other objects documented here.
\details{
Insert this function into a binary regression model
}
%%\value{
%%}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Thomas Welchowski \email{welchow@imbie.meb.uni-bonn.de}

Matthias Schmid \email{matthias.schmid@imbie.uni-bonn.de}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{glm}}
}
\examples{
# Example with copenhagen stroke study
library(pec)
data(cost)
head(cost)

# Take subsample and convert time to months
costSub <- cost [1:50, ]
costSub$time <- ceiling(costSub$time/30)
costLong <- dataLong(dataSet=costSub, timeColumn="time", censColumn="status")
gumbelModel <- glm(formula=y ~ timeInt + diabetes, data=costLong, family=binomial(link=gumbel()))

# Estimate hazard given prevStroke nad no prevStroke
hazPrevStroke <- predict(gumbelModel, newdata=data.frame(timeInt=factor(1:143), 
diabetes=factor(rep("yes", 143), levels=c("no", "yes"))), type="response")
hazWoPrevStroke <- predict(gumbelModel, newdata=data.frame(timeInt=factor(1:143), 
diabetes=factor(rep("no", 143), levels=c("no", "yes"))), type="response")

# Estimate survival function
SurvPrevStroke <- cumprod(1-hazPrevStroke)
SurvWoPrevStroke <- cumprod(1-hazWoPrevStroke)

# Example graphics of survival curves with and without diabetes
plot(x=1:143, y=SurvWoPrevStroke, type="l", xlab="Months", 
ylab="S (t|x)", las=1, lwd=2, ylim=c(0,1))
lines(x=1:143, y=SurvPrevStroke, col="red", lwd=2)
legend("topright", legend=c("Without diabetes", "Diabetes"), 
lty=1, lwd=2, col=c("black", "red"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line