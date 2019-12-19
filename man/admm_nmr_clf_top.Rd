\name{admm_nmr_clf_top}
\alias{admm_nmr_clf_top}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
admm_nmr_clf_top(A, B, label, top_num)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{label}{
%%     ~~Describe \code{label} here~~
}
  \item{top_num}{
%%     ~~Describe \code{top_num} here~~
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
\author{
%%  ~~who you are~~
}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
