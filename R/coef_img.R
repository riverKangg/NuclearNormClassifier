coef_img <- function(A,x) {
  p<-dim(A)[1]; q<-dim(A)[2]; n<-dim(A)[3];
  A_x <- matrix(rep(0, p*q), dim=c(p,q))
  for (i in 1:n){
    A_x <- A_x + x[i]*A[,,i]
  }
  return(A_x)
}
