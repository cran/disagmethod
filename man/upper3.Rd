\name{upper3}
\alias{upper3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Upper Bound Disaggregation Method Function
}
\description{
This uses the Upper Bound method for temporal disaggregation of time series
}
\usage{
upper3(x, m = 1)
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
This function uses the upper bound method found in  Hodgess and Wei (1996, "Temporal Disaggregation of Time Series"). We fit an aggregate (p,d,q) model, and produce a disaggregate model of (p,d,(p+d)).  We generate the disaggregate series based on the disaggregate model.
}
\value{
%%  If it is a LIST, use
\item{bigy }{order of the disaggregate model}
\item{fin1 }{final disaggregate series}
}
\references{
Hodgess and Wei (1996, "Temporal Disaggregation of Time Series"), in M. Ahsanullah and D. Bhoj (Eds), "Applied Statistical Science I".
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
e.upp <- upper3(my.un.q,3)
sum(e.upp$fin1[1:3])
my.un.q[1]
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ts}
\keyword{models}
% Use only one keyword per line.

