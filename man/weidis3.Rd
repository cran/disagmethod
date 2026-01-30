\name{weidis3}
\alias{weidis3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wei Stram Disaggregation Method Function
}
\description{
This uses the Wei Stram method for temporal disaggregation of time series
}
\usage{
weidis3(x, m = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Aggregate Series; must be a ts, xts, or zoo object
}
  \item{m}{
order of disaggregation; 3, 4, 12
}
}
\details{
This function uses the method found in Wei and Stram (1990, <doi:10.1111/j.2517-6161.1990.tb01799.x>).  We fit an aggregate (p,d,q) model, and produce a disaggregate model of (p,d,(p+d+1)).  We generate the disaggregate series based on the disaggregate model.
}
\value{
%%  If it is a LIST, use
\item{bigy }{order of the disaggregate model}
\item{fin1 }{final disaggregate series}
}
\references{
Wei and Stram (1990, <doi:10.1111/j.2517-6161.1990.tb01799.x>)
}
\author{
Erin Hodgess
}

\examples{
library(tswge)
data(tx.unemp.adj)
#Monthly seasonally adjusted Texas unemployment data
#Create a quarterly sum
my.un.q <- aggregate(tx.unemp.adj,nfreq=4)
e.wei <- weidis3(my.un.q,3)
sum(e.wei$fin1[1:3])
my.un.q[1]
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ts}
\keyword{models}
% Use only one keyword per line.

