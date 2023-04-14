predict.mars <- function(object,newdata) {
  if(missing(newdata) || is.null(newdata)) {
    B <- as.matrix(object$B)
  }
  else {
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  drop(B %*% beta)
}


make_B <- function(B, Bfuncs){
  B.output <- matrix(0,nrow=nrow(B),ncol=length(Bfuncs))
  for(i in 1:length(Bfuncs)) {
    B.output[,i] <- out.make_B(B,Bfuncs[[i]])
    B.output
  }
}


out.make_B <- function(B,s) {
  B.output <- rep(1,nrow(B))
  if(nrow(s)>1) {
    for(i in 2:nrow(s)){
      B.output <- B.output * h(B[,s[i,"v"]],s[i,"s"],s[i,"t"])
    }
  }
  B.output
}

