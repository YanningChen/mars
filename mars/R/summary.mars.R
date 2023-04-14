## summary()

summary.mars <- function(object){
  xn <- object$x_names
  bn <- names(object$coefficients)
  knots <- object$Bfuncs
  cat("Basis functions:\nB0:\n  Intercept\n")
  for(i in 2:length(knots)) {
    cat(paste0(bn[i],":\n"))
    for(j in 1:nrow(knots[[i]])) {
      cat(paste0("  Component ",j,"; variable: ",xn[knots[[i]][j,"v"]],";"))
      cat(paste0(" sign: ",knots[[i]][j,"s"],";"))
      cat(paste0(" knots at: ",knots[[i]][j,"t"],"\n"))
    }
  }
  summary=summary.lm(object)
  return(summary)
}

