print.mars <- function(object) {
  print("mars object")
  print(object$call)
  cat("coefficients \n",object$coefficients,"\n")
  print("knots")
  print(object$Bfuncs)
  print("Model(first 6)")
  colnames(object$model)<-c("y","intercept",object$x_names)
  print(head(object$model))

}





