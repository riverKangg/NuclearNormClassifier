\name{coef_img}
\alias{linear combination}
\title{coefficient image}

\description{
  This function performs linear combination of train images.
}
\usage{
  coef_img(A, x)
}

\arguments{
  \item{A}{
  "A" is a train data set in the form of a three-dimensional array.
  This data set contains several face images used for algorithm learning.
  }
  \item{x}{
  "x" is a coefficient vector for the train images.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
\link{https://ieeexplore.ieee.org/document/7420697}
}
\author{Jisun Kang}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (A, x)
{
    p <- dim(A)[1]
    q <- dim(A)[2]
    n <- dim(A)[3]
    A_x <- matrix(rep(0, p * q), dim = c(p, q))
    for (i in 1:n) {
        A_x <- A_x + x[i] * A[, , i]
    }
    return(A_x)
  }
}

\keyword{Linear combination}
