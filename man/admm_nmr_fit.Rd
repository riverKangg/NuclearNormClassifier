\name{admm_nmr_fit}
\alias{admm_nmr_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
This function fit the ADMM algorithm for the nuclear norm regression.
It returns the optimal regression coefficient vector.
You can use the output of this function for classification.
}
\usage{
admm_nmr_fit(max_iter, ep_abs, ep_rel, A, B)
}

\arguments{
  \item{max_iter}{
  Maximum number of iterations
  }
  \item{ep_abs}{
  Parameter used for termination conditions
  }
  \item{ep_rel}{
  Parameter used for termination conditions
  }
  \item{A}{
  "A" is a train data set in the form of a three-dimensional array.
  This data set contains several face images used for algorithm learning.
  }
  \item{B}{
  "B" is an image of a face with occlusion.
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

\seealso{
\code{\link{https://web.stanford.edu/~boyd/papers/pdf/admm_slides.pdf}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (max_iter, ep_abs, ep_rel, A, B)
{
    p <- dim(A)[1]
    q <- dim(A)[2]
    n <- dim(A)[3]
    h <- array(rep(0, p * q), dim = c(p * q, n))
    for (n_iter in 1:n) {
        h[, n_iter] <- as.vector(A[, , n_iter])
        M <- solve(t(H) \%*\% H + lambd/mu) \%*\% t(H)
    }
    Y <- -B
    Z <- 0
    k <- 0
    x <- rep(0, n)
    for (i in 1:max_iter) {
        g_in <- B + Y - 1/mu * Z
        g <- as.vector(g_in)
        x_new <- M \%*\% g
        A_x_new <- coef_img(A, x_new)
        Y_new <- svt(mu, A_x_new, B, Z)
        Z_new <- Z + mu * (A_x_new - Y_new - B)
        r_pri <- A_x_new - Y_new - B
        ep_pri <- ep_abs * sqrt(p * q) + ep_rel
        s_dual <- 1/mu * h \%*\% as.vector(Y - Y_new)
        ep_dual <- ep_abs * sqrt(n) + ep_rel
    }
    x <- x_new
    Y <- Y_new
    Z <- Z_new
    A_x <- A_x_new
    return(X_new)
  }
}

\keyword{The alternating direction method of multipliers}
\keyword{The nuclear norm regression}
