\name{plot.discSurvMartingaleResid}
\alias{plot.discSurvMartingaleResid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot method of class "discSurvMartingaleResid"
}
\description{
Gives a different plot of each marginal covariate against the martingale residuals. Additionally a nonparametric \code{\link{loess}}  estimation is done. 
}
\usage{
\method{plot}{discSurvMartingaleResid}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class "discSurvMartingaleResid"
}
  \item{\dots}{
Additional arguments to the plot function
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
Thomas Welchowski
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{martingaleResid}}
}
%%\examples{
%%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
