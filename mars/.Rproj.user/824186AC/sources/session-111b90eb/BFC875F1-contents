#' Multivariate Adaptive Regression Splines (MARS)
#' @description Mars is a form of stepwise linear regression to improve
#'  performance of the model in the regression setting.
#'
#' @param formula an R formula
#' @param data a data frame containing the data for the model
#' @param control an object of class 'mars.control'
#' @details The MARS function is applied for constructing Multivariate
#' Adaptive Regression Splines models. It refers to a stepwise linear
#' regression format that produces the required model. As it proceeds
#' it is ensured that for example, it will first look for the single
#' point across the range of X values where two different linear relationships
#' between Y and X achieve the smallest error via GCV computation. The "mars"
#' function takes as input a formula representing the relationship between the
#' response and explanatory variables, the data set to fit the model on and a
#' control object. It returns an object of class "mars". This returned object
#' is later used for further output results on different methods.
#'
#' The mars function comprises several methods, including anova.mars,
#' plot.mars, predict.mars, print.mars, and summary.mars. These methods
#' generate outputs through the application of the corresponding
#' techniques - anova, plot, predict, print, and summary.
#'
#' @return an object of class 'mars' that include the final regression
#' and a description of the basis functions. There are plot, predict, anova,
#' summary and print methods for mars objects.
#' @export
#' @author Jiawei Lin, Yanning Chen, Isabah Tabenda Hasan
#'
#' @references
#' Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
#' The Annals of Statistics, 19(1), 1–67. http://www.jstor.org/stable/2241837
#'
#' Boehmke, B. (2020, February). Chapter 7 Multivariate Adaptive Regression
#' Splines | Hands-On Machine Learning with R.
#' Github.io. https://bradleyboehmke.github.io/HOML/mars.html
#'
#' @seealso anova, predict, print, plot, summary
#'
#' @examples
#' library(ISLR)
#' data("Auto")
#' data("College")
#' data("Hitters")
#' # example1:
#' Auto_mars<- mars(Auto$mpg~Auto$displacement+Auto$weight, data=Auto)
#' # example2:
#' College_mars=mars(College$Grad.Rate~College$F.Undergrad+College$P.Undergrad,data=College)
#' # example3:
#' Hitters_mars=mars(Hitters$Salary~Hitters$HmRun+Hitters$Hits, data=Hitters)
#'
mars <- function(formula,data,control=mars.control()) {
  cc <- match.call() # save the call
  mf <- model.frame(formula,data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)[,-1,drop=FALSE]
  x_names <- colnames(x)
  control <- validate_mars.control(control)
  fwd <- fwd_stepwise(y,x,control)
  bwd <- bwd_stepwise(fwd,control)
  fit <- lm(y~.-1,data=data.frame(y=y,bwd$B)) # notice -1 added
  out <- c(list(call=cc,formula=formula,y=y,B=bwd$B,Bfuncs=bwd$Bfuncs,
                x_names=x_names),fit)
  class(out) <- c("mars",class(fit))
  out
}


fwd_stepwise <- function(y,x,control=mars.control()){
  N <- length(y) # sample size
  n <- ncol(x) # number of predictors
  B <- init_B(N,control$Mmax)
  Bfuncs<-vector(mode="list",length=control$Mmax+1)
for(i in 1:(control$Mmax/2)) {
    M<-2*i-1
    if(control$trace) cat("M",M,"\n")
    lof_best <- Inf
    for(m in 1:M) {
      vset <- setdiff(1:n,Bfuncs[[m]][,"v"])
      if(control$trace) cat("M, m, vset", M,m,vset,"\n")
      for(v in vset){ # select a variable to split on
        tt <- split_points(x[,v],B[,m])
        for(t in tt) {
          Bnew <- data.frame(B[,1:M],
                             Btem1=B[,m]*h(x[,v],-1,t),
                             Btem2=B[,m]*h(x[,v],1,t))
                    gdat <- data.frame(y=y,Bnew)
                    lof <- LOF(y~.,gdat,control)
                    if(lof < lof_best) {
                      lof_best <- lof
                      best_split <- c(m=m,v=v,t=t)
          }
        }
      }
    }

    mstar<-best_split["m"]
    vstar<-best_split["v"]
    tstar<-best_split["t"]

    Bfuncs[[M+1]]<-rbind(Bfuncs[[mstar]],c(s=-1,vstar,tstar))
    Bfuncs[[M+2]]<-rbind(Bfuncs[[mstar]],c(s=+1,vstar,tstar))

    B[,M+1:2] <- cbind(B[,mstar]*h(x[,vstar],-1,tstar),B[,mstar]*h(x[,vstar],+1,tstar))
  }

  colnames(B) <- paste0("B",(0:(ncol(B)-1))) # optional
  return(list(y=y, B=B, Bfuncs=Bfuncs))
}



init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}


bwd_stepwise<- function(fwd,control){
  #fwd is a list with elements y, B and Bfuncs
  Mmax <- ncol(fwd$B)-1
  Jstar <- 2:(Mmax+1) #Note fwd$B[,Jstar] doesn't include 1st column
  Kstar <- Jstar
  dat <-data.frame(y=fwd$y,fwd$B)
  lofstar <- LOF(y~.-1,dat,control)
  for(M in (Mmax+1):2){
    b <- Inf
    L <- Kstar
    if(control$trace) cat("L:",L,"/n")
    for(m in L){
      K <- setdiff(L,m)
      #L <- # Rcode
      dat <- data.frame(y=fwd$y,fwd$B[,K])

      lof <- LOF(y~.,dat,control)
      # Rcode: to update Kstar and Jstar
      # Jstar lists the set of remaining columns (predictors)
      if(control$trace) cat("M:K:lof",M,":",K,":",lof,"\n")
      if(lof<b) {
        b <- lof
        Kstar <- K
      }
      if(lof<lofstar){
        lofstar<- lof
        Jstar <- K
      }
    }
    if(control$trace) cat("M:Jstar:lofstar",M,":",Jstar,":",lofstar,"\n")
  }

  Jstar <- c(1,Jstar)
  return(list(y=fwd$y,B=fwd$B[,Jstar],Bfuncs=fwd$Bfuncs[Jstar]))
}

LOF <- function(form,data,control) {
   # update this LOF to GCV
  mod<-lm(form,data)
  RSS<-sum((mod$res)^2)
  N=nrow(data)
  M=length(coef(mod))-1
  Ctilde<-sum(diag(hatvalues(mod)))+control$d*M
  return(RSS*N/(N-Ctilde)^2)
}

h <- function(x,s,t) {
    return(pmax(0,s*(x-t)))
  }
  # if x>t, s=+1, this return max(0,x-t)
  # if x<t, s=-1, this return max(0,t-x)



split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
  ...
}


#------------------------------------------------------------------------
# constructor, validator and helper for class mars.control
#------------------------------------------------------------------------
#
new_mars.control <- function(control) {
  structure(control, class="mars.control")
}

validate_mars.control <- function(control) {
  stopifnot(is.integer(control$Mmax), is.numeric(control$d),is.logical(control$trace))
  if(control$Mmax<2){
    warning("Mmax must be >=2; setting to 2")
    control$Mmx<-2
  }
  if(control$Mmax %% 2 >0){
    control$Mmax<-2*ceiling(control$Mmax/2)
    warning("Mmax should be an even integer. Setting to ", control$Mmax)
  }
  control
}



mars.control <- function(Mmax=2,d=3,trace=FALSE) {
  Mmax <- as.integer(Mmax)
  control <- list(Mmax=Mmax,d=d,trace=trace)
  control <- validate_mars.control(control)
  new_mars.control(control)
}

