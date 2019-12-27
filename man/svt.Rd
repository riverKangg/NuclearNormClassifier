\name{svt}
\alias{svt}
\title{Singular value thresholding}

\description{
  This function performs singular value thresholding.
  Singular value thresholding is used to update y that cannot be represented in closed form.
  Entering the parameters mu, A (x), B, and Z returns the y after performing singular value thresholding.
}

\usage{svt(mu, A_x, B, Z)}

\arguments{
  \item{mu}{
  "mu" is a parameter that must be set when defining a function.
  }
  \item{A_x}{
  "A_x" is a result obtained by correlation coefficient x with A, a test image data set.
  This can be made using the coef_img function.
  }
  \item{B}{
  "B" is an image of a face with occlusion.
  }
  \item{Z}{
  "Z" can be found using A, B, and y. It can be calculated naturally during algorithm implementation.
  }
}

\details{}

\value{
  \item{svt_out}{
  "svt_out" is the y after performing singular value thresholding.
  }
}

\references{
  \link{https://ieeexplore.ieee.org/document/7420697}
}

\author{Jisun Kang}

\seealso{
  \code{\link{https://arxiv.org/abs/0810.3286}}
}

\examples{

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
