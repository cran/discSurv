\name{tauToPearson}
\alias{tauToPearson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transformation of Tau to Pearson correlation
}
\description{
If the two variables are bivariate normal distributed, this formula maps kendalls tau to bravais pearson correlation coefficient rho.
}
\usage{
tauToPearson(Tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tau}{
Numeric vector in the interval [-1, 1]
}
}
\details{
This relationship holds only in the case if the two random variables are jointly normal distributed.
}
\value{
Vector of bravais pearson correlation coefficient rho
}
\references{
William H. Kruskal, (1958),
\emph{Ordinal Measures of Association},
Journal of the American Statistical Association, Vol. 53, No. 284, pp. 814-861
%@article{OrdMeas,
%  title={Ordinal Measures of Association},
%  author={William H. Kruskal},
%  journal={Journal of the American Statistical Association, Vol. 53, No. 284, pp. 814-861},
%  year={1958}
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
\code{\link{simCompRisk}}
}
\examples{
# Plot of relation between kendalls tau and pearson correlation
MaxDisc1 <- nlminb(start=0.5, objective=function (x) -abs(tauToPearson (x)-x))
MaxDisc2 <- nlminb(start=-0.5, objective=function (x) -abs(tauToPearson (x)-x))
plot(x=seq(-1,1,length.out=500), y=tauToPearson(seq(-1,1,length.out=500)), 
xlab=expression(tau), ylab=expression(rho), type="l", las=1, 
main="Relationship between tau and pearson correlation (bivariate normal)", lwd=2)
lines(x=seq(-1,1,length.out=500), y=seq(-1,1,length.out=500), lty=2)
segments(x0=0, y0=-1.25, x1=0, y1=0, lty=2)
segments(x0=-1.25, y0=0, x1=0, y1=0, lty=2)
segments(x0=MaxDisc1$par, y0=-1.25, x1=MaxDisc1$par, y1=tauToPearson (MaxDisc1$par), lty=2)
segments(x0=MaxDisc2$par, y0=-1.25, x1=MaxDisc2$par, y1=tauToPearson (MaxDisc2$par), lty=2)

# The maximum discrepancy between pearson and spearman is at 
# a kendalls tau value about 0.56 and -0.56
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line