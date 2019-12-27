\name{admm_nmr_clf_top}
\alias{admm_nmr_clf_top}
\title{}

\description{

}
\usage{
  admm_nmr_clf_top(A, B, label, top_num)
}

\arguments{
  \item{A}{
  "A" is a train data set in the form of a three-dimensional array.
  This data set contains several face images used for algorithm learning.
  }
  \item{B}{
  "B" is an image of a face with occlusion.
  }
  \item{label}{
  "label" is a vector representing a label for train dataset A.
  }
  \item{top_num}{
  "top_num" represents the number of classification candidates.
  If "top_num" is 1, this function is the same as admm_nmr_clf.
  }
}
\details{}
\value{
  \item{val_min[0]}{

  }
  \item{val_top}{

  }
}
\references{
  \link{https://ieeexplore.ieee.org/document/7420697}
}
\author{Jisun Kang}
\note{}
\seealso{}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (A, B, label, top_num)
{
    x <- ADMM_NMR_fit(A, B)
    A_x <- coef_img(A, x)
    label_unique <- unique(label)
    error_lst <- c()
    for (i in target_unique) {
        label_tf <- (label == label_unique)
        x_i <- x[label_tf]
        A_i <- A[, , label_tf]
        A_x_i <- coef_img(A_i, x_i)
        diff <- A_x - A_x_i
        error <- nuclear(diff)
        error_lst <- append(error_lst, error)
    }
    error_top <- sort(error_lst)[1:top_num]
    error_top_tf <- is.element(error_lst, error_top)
    error_min_tf <- (error_lst == min(error))
    val_min <- label[error_top_tf]
    return(val_min[0], val_top)
  }
}

\keyword{Image classification}
\keyword{top-n error}
