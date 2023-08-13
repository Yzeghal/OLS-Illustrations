library(MASS)  #mvnorm
library(stats) #summarise regressions with asymptotic sd estimators (non heteroscedastic robust, use sandwich library for that)
library(matrixcalc)
library(AER)


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

V = diag(3) # variance matrix 
E = mvrnorm(n=50000, mu = c(0,0,0), Sigma = V)

e1 = E[,1]
e2 = E[,2]
e3 = E[,3]

# e1,e2,e3 are observations of 3  random variables that form an orthonormal base of a 3 dimension subspace of RV in |R

# 1st step : generating y

y = e3+e2
reg1 = lm(y~e2-1) # -1 specifies not to add a constant to the model.
reg1
x1 = 1/sqrt(2)*e1+1/sqrt(2)*e2 # X1 is a random variable of L2 norm 1.
reg2 = lm(y~e2+x1-1)
reg2

instrument = cbind(e2,e1)
Pinstru = project(instrument)
P1 = Pinstru %*% cbind(e2,x1)
res1 = instrument %*% P1

sum(abs(res1[,1]-e2)>1e-13)
sum(abs(res1[,2]-sqrt(2)/2*e3)>1e-13)
coefs(y,res1)
ivreg(y~x1 + e2-1 | e2 + e1-1)$coefficients

lm(y~e2+x1-1)



