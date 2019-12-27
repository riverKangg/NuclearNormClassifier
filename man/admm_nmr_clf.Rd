\name{admm_nmr_clf}
\alias{admm_nmr_clf}
\title{Classification using nuclear norm regression}

\description{
  Entering the training data set and the test data returns the optimal label.
}

\usage{
  admm_nmr_clf(A, B, label)
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
}

\value{
  \item{val}{
  The result of classification. It represents the optimal label.
  }
}

\references{
  \link{https://ieeexplore.ieee.org/document/7420697}
}
\author{Jisun Kang}

\examples{

function (A, B, label)
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
    error_tf <- error == min(error)
    val <- label_unique[error_tf]
    return(val)
  }
}

\keyword{Image data}
\keyword{Classification}
