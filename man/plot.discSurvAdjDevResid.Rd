\name{plot.discSurvAdjDevResid}
\alias{plot.discSurvAdjDevResid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method of class "discSurvAdjDevResid"
}
\description{
Is called implicitly by using standard plot function on an object of class "discSurvAdjDevResid". It plots a qqplot against the normal distribution. If the model fits the data well, it should be approximately normal distributed.
}
\usage{
\method{plot}{discSurvAdjDevResid}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class "discSurvAdjDevResid"
}
\item{...}{
Additional arguments to plot function
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}
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
\code{\link{adjDevResid}}
}
%%\examples{
%%
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line