
coef_img <- function(A,x) {
  p<-dim(A)[1]; q<-dim(A)[2]; n<-dim(A)[3];
  A_x <- matrix(rep(0, p*q), dim=c(p,q))
  for (i in 1:n){
    A_x <- A_x + x[i]*A[,,i]
  }
  return(A_x)
}

# vec 모듈 생략, as.vector 사용

svt <- function(mu,A_x,B,Z) {
  svt_in <- A_x-B+1/mu*Z
  result_svd <- svd(svt_in)
  # thresholding
  d <- pmax(result_svd$d-1/mu, 0)
  # return
  svt_out <- result_svd$u %*% d %*% result_svd$v
  return(svt_out)
}


###############################################################
########################## FAST ADMM ##########################
###############################################################

fast_admm_nmr_fit <- function(){
  p<-dim(A)[1]; q<-dim(A)[2]; n<-dim(A)[3];
  # step 1
  h <- array(rep(0,p*q),dim=c(p*q,n))
  theta <- lambd*(1+gam)
  M <- solve(t(H)%*%H+lambd/mu)%*%t(H)

  # step 2
  x<-rep(0,n); xhat<-rep(0,n)
  Z<-0; Zhat<-0; Y<-0; alpha<-1; k<-0;

  # step 3~5
  for (i in 1:max_iter){
    # updating Y
    A_xhat <- coef_img(A,xhat)
    Y_new <- svt(A_xhat,B,Z)

    # updating x
    g_in <- B+Y_new-1/mu*Zhat
    g <- as.vector(g_in)
    x_new <- M%*%g

    # updating Z
    A_x_new <- coef_img(A,x_new)
    Z_new <- Zhat + mu*(A_x_new-Y_new-B)

    # updating alpha
    alpha_new <- (1+sqrt(1+4*alpha*alpha))/2

    # updating xhat
    xhat_new <- x+(alpha-1)/alpha_new*(x_new-x)

    # updating Zhat
    Zhat_new <- Z+(alpha-1)/alpha_new*(Z_new-Z)

    # termination
    r_pri <- A_x_new-Y_new-B
    ep_pri <- ep_abs*sqrt(p*q)+ep_rel*max(c(nuclear(A_x_new),nuclear(Y_new),nuclear(B)))
    s_dual <- 1/mu*t(H) %*% as.vector(Y-Y_new)
    ep_dual <- ep_abs*sqrt(n) + ep_rel


  }

  x<-x_new; xhat<-xhat_new; Y<-Y_new; Z<-Z_new; Zhat<-Zhat_new;

}



########################## clf ##########################

fast_admm_nmr_clf <- function(A,B,label){
  x <- FAST_ADMM_NMR_fit(A,B)
  A_x <- coef_img(A,x)

  # A_x_i 구하기 + error 구하기
  label_unique <- unique(label)
  error_lst <- c()
  for (i in target_unique){
    label_tf <- (label==label_unique)
    x_i <- x[label_tf]
    A_i <- A[,,label_tf]
    A_x_i <- coef_img(A_i,x_i)

    diff <- A_x-A_x_i
    error <- nuclear(diff)
    error_lst <- append(error_lst,error)
  }
  # min error
  error_tf <- error==min(error)
  val <- label_unique[error_tf]

  # return
  return(val)
}


####################### clf top #######################

fast_admm_nmr_clf_top <- function(A,B,label,top_num){
  x <- FAST_ADMM_NMR_fit(A,B)
  A_x <- coef_img(A,x)
  # A_x_i 구하기 + error 구하기
  label_unique <- unique(label)
  error_lst <- c()
  for (i in target_unique){
    label_tf <- (label==label_unique)
    x_i <- x[label_tf]
    A_i <- A[,,label_tf]
    A_x_i <- coef_img(A_i,x_i)

    diff <- A_x-A_x_i
    error <- nuclear(diff)
    error_lst <- append(error_lst,error)
  }

  # top_num
  error_top <- sort(error_lst)[1:top_num]
  error_top_tf <- is.element(error_lst,error_top)

  # min error
  error_min_tf <- (error_lst==min(error))
  val_min <- label[error_top_tf]

  # return
  return(val_min[0],val_top)
}
