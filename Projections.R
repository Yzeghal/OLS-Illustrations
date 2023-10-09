library(matrixcalc)

project<-function(X, threshold = 1e-13){
  M = t(X)%*%X
  if (abs(det(M))<threshold){
    stop("Non invertible matrix X'X")
  }
  mat = matrix.inverse(t(X)%*%X)%*%t(X)
  return(mat)
}

coefs<-function(Y,X, threshold = 1e-13){
  mat = project(X, threshold)
  return(mat%*%Y)
}

#Coded a R2 that is faster to call than $r.squared
R2 = function(m){
  #takes a lm object and yields the Rsquared
  if (class(m)!="lm"){
    stop("Input i not a linear model")
  }
  x = m$model[,1]
  SCT = sd(x)^2
  SCR = sd(m$residuals)^2
  R2 = 1-SCR/SCT
  return(R2)
}