\name{svt}
\alias{singular value thresholding}
\title{This function performs singular value thresholding. }
\description{This function performs singular value thresholding. Singular value thresholding is used to update y that cannot be represented in closed form. Entering the parameters mu, A (x), B, and Z returns the y after performing singular value thresholding.}
\usage{svt(mu, A_x, B, Z)}
\arguments{
  \item{mu}{
Parameter to be set when defining a function.
}
  \item{A_x}{
A_x is a result obtained by correlation coefficient x with A, a test image data set. This can be made using the coef_img function.
}
  \item{B}{
B represents an image of a face with occlusion.
}
  \item{Z}{
Z can be found using A, B, and y. It can be calculated naturally during algorithm implementation.
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
%% ~put references to the literature/web site here ~
}
\author{Jisun Kang}
\note{
%%  ~~further notes~~
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mu, A_x, B, Z)
{
    svt_in <- A_x - B + 1/mu * Z
    result_svd <- svd(svt_in)
    d <- pmax(result_svd$d - 1/mu, 0)
    svt_out <- result_svd$u \%*\% d \%*\% result_svd$v
    return(svt_out)
  }
}
\keyword{Singular value thresdholding}
