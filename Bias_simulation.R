#####
# Script containing regressions of simulated variables to illustrate included and excluded variable bias.
# In the first part, we construct included variable bias
# In the second part, we construct ommitted variable bias
# ...
#####
library(MASS)  #mvnorm
library(stats) #summarise regressions with asymptotic sd estimators (non heteroscedastic robust, use sandwich library for that)
library(matrixcalc)  # Used for matrix inverse
library (matlib)
library(rgl)#usedfor 3d representation
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

###################################################
# Regression of centered variables without constant.
set.seed(18)
V = diag(3) #variance matrix : base will be orthogonal and normed vectors.
X = mvrnorm(n=50000, mu = c(0,0,0), Sigma = V)

x1 = X[,1]
x2 = X[,2]
d = rep(1,50000)

y = 1/sqrt(2)*x1 + 1/sqrt(2)*x2
mean(y)     # means of variables are "almost 0"
mean(x1)
mean(x2)
reg_centered = lm(y~x1+x2-1)  #adding -1 in the formula specifies that
#no constant should be added to the regressors
reg_centered

# Regression of non centred variables without constant
yo = y+400
x1o = x1+100
x2o = x2+300
mean(yo)    # variables a no longer centered
mean(x1o)
mean(x2o)
reg_offset = lm(yo~x1o+x2o-1)

reg_offset
sum(reg_offset$coefficients * c(100,300))

# Regression of non centred  variables with a constant

reg_constant = lm(yo~x1o+x2o+d-1) # Remember that x3 contains ones
reg_constant






###########################
# Illustrating the use of a constant in 3D

base <- diag(3)
rownames(base) <- c("x1", "x2", "1")

vec = matrix(1/sqrt(2)*base[1,]+sqrt(1-1/sqrt(2)^2)*base[2,], nrow = 1)
rownames(vec) = "Y"
open3d()
planes3d(0, 0, 1, 0, col="gray", alpha=0.2)

vectors3d(base, col = "black")
vectors3d(vec, col = "red")
segments3d(rbind(c(1/sqrt(2),0,0),c(1/sqrt(2),1/sqrt(2),0)))
segments3d(rbind(c(0,1/sqrt(2),0),c(1/sqrt(2),1/sqrt(2),0)))
# close3d()
open3d()

names(base) = NULL
vecX1O = matrix(base[1,]+1*base[3,], nrow = 1)
vecX2O = matrix(base[2,]+3*base[3,], nrow = 1)
vecY = matrix(1/sqrt(2)*base[1,]+1/sqrt(2)*base[2,]+4*base[3,], nrow = 1)
vec = rbind(vecX1O, vecX2O, vecY)
rownames(vec) = c("x1o", "x2o", "yo")

planes3d(0, 0, 1, 0, col="gray", alpha=0.2)
vectors3d(base, col = "black")
vectors3d(vec, col = "blue")
segments3d(rbind(c(1/sqrt(2),0,0),c(1/sqrt(2),1/sqrt(2),0)))
segments3d(rbind(c(0,1/sqrt(2),0),c(1/sqrt(2),1/sqrt(2),0)))
segments3d(rbind(c(0,1,0),c(0,1,3)))
segments3d(rbind(c(1,0,0),c(1,0,1)))
segments3d(rbind(c(1/sqrt(2),1/sqrt(2),0),c(1/sqrt(2),1/sqrt(2),4)))

open3d()

names(base) = NULL
vecX1O = matrix(base[1,]+2*base[3,], nrow = 1)
vecX2O = matrix(base[2,]+6*base[3,], nrow = 1)
vecY = matrix(1/sqrt(2)*base[1,]+1/sqrt(2)*base[2,]+8*base[3,], nrow = 1)
vec = rbind(vecX1O, vecX2O, vecY)
rownames(vec) = c("x1o", "x2o", "yo")

planes3d(0, 0, 1, 0, col="gray", alpha=0.2)
vectors3d(base, col = "black")
vectors3d(vec, col = "green")
segments3d(rbind(c(1/sqrt(2),0,0),c(1/sqrt(2),1/sqrt(2),0)))
segments3d(rbind(c(0,1/sqrt(2),0),c(1/sqrt(2),1/sqrt(2),0)))
segments3d(rbind(c(0,1,0),c(0,1,6)))
segments3d(rbind(c(1,0,0),c(1,0,2)))
segments3d(rbind(c(1/sqrt(2),1/sqrt(2),0),c(1/sqrt(2),1/sqrt(2),8)))

close3d()
close3d()
close3d()


#Settings 
#graph_1

#####
vec <- diag(3)
rownames(vec) <- c("e1", "e2", "e3")
open3d()
vectors3d(vec, color=rep("black",3), lwd=2, cex.lab =2)
# draw the XZ plane, whose equation is Y=0
planes3d(0, 0, 1, 0, col="gray", alpha=0.2)


#####
# First regression 
set.seed(18)
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
R2(reg1)          # ~0.5 as expected : ||e2||^2/||e2+e3||^2

# 3D representation of reg1


open3d()
planes3d(0, 0, 1, 0, col="lightgrey")

vectors3d(base[c(1,3),], col = "black",  lwd = 2)
vectors3d(c(0,1,1), col ="red",  lwd = 2, labels = "Y")
vectors3d(base[2,], col = "blue4",  lwd = 2, labels ="e2")
segments3d(rbind(c(0,1,0),c(0,1,1)), col = "red")
corner(c(0,0,0),c(0,1,0),c(0,1,1), col = "red")
light3d(0, 0)
light3d(0, 0)
close3d()



# 2nd step : adding a correlated variable does not necessarily generate bias if it is orthogonal to epsilon

x1 = 1/sqrt(2)*e1+1/sqrt(2)*e2 # X1 is a random variable of L2 norm 1.
mean(x1^2) # estimates the norm of X1
cor(x1,e2)
reg2 = lm(y~e2+x1-1)
reg2
R2(reg2)
base <- diag(3)
rownames(base) <- c("e1", "e2", "e3")

# 3D representation of reg2

vec = matrix(1/sqrt(2)*base[1,]+1/sqrt(2)*base[2,], nrow = 1)

rownames(vec) = "X1"
open3d()
planes3d(0, 0, 1, 0, col="cornflowerblue")

vectors3d(base[c(1,3),], col = "black",  lwd = 2)
vectors3d(c(0,1,1), col ="red",  lwd = 2, labels = "Y")
vectors3d(base[2,], col = "blue4",  lwd = 2, labels ="e2")
vectors3d(vec, col = "blue4", lwd = 2)
segments3d(rbind(c(0,1,0),c(0,1,1)), col = "red")
corner(c(0,0,0),c(0,1,0),c(0,1,1), col = "red")

light3d(0, 0)
light3d(0, 0)
light3d(0, 0)
close3d()

# 3rd step : A control orthogonal to e2 cannnot change its coefficient.

xo = 1/sqrt(2)*e1+1/sqrt(2)*e3

cor(xo,e2)  # "Almost" 0

reg3 = lm(y~e2+xo-1)
reg3

R2(reg3)
vec =rbind( 1/sqrt(2)*c(1,0,1), c(0,1,0)) %*% base
rownames(vec) = c("Xo", "e2")
pointfall = t(vec)%*%project(t(vec))%*% matrix(c(0,1,1), nrow = 3)
open3d()
planes3d(0, 0, 1, 0, col="lightgrey")
planes3d(1, 0,-1, 0, col = "cornflowerblue")
vectors3d(base[c(1,3),], col = "black",  lwd = 2)
vectors3d(vec, col = "blue4", lwd = 2)
vectors3d(c(0,1,1),label = "Y", col = "red", lwd = 2)
segments3d(rbind(c(0,0,0),t(pointfall)), col = "red")
segments3d(rbind(t(pointfall),c(0,1,1)), col = "red")
segments3d(rbind(t(pointfall),c(0,1,0)), col = "black")
corner(c(0,0,0),as.vector(pointfall),c(0,1,1), col = "red")
corner(c(0,0,0),c(0,1,0),as.vector(pointfall), col = "black")

light3d(0,0)
light3d(0,0)

open3d()
planes3d(0, 0, 1, 0, col="lightgrey")
vectors3d(base[c(1,3),], col = "black",  lwd = 2)
vectors3d(vec, col = "blue4", lwd = 2)
vectors3d(c(0,1,1),label = "Y", col = "red", lwd = 2)
segments3d(rbind(c(0,1,0),c(0,1,1)), col = "red")
corner(c(0,0,0),c(0,1,0),c(0,1,1), col = "red")
light3d(0, 0)
light3d(0, 0)

close3d()
close3d()

# 3rd step : adding more noise to x1, that would reduce actual noise.
x1e = x1+e3 
reg3 = lm(y~e2+x1e)
# Here, giving a non 0 coeff to x1 enables us to reduce the mean sq error through e3, even though it requires to increase e1 and e2 coefficients.
# Nothing can be done for e1, but e2 then compensate for the excess due to this minimisation choice.

summary(reg3)
R2(reg3)

#Theoretical counterpart :
A = matrix( c(sqrt(2)/2,sqrt(2)/2,1,0,1,0), nrow = 3) 
M = matrix.inverse(t(A)%*%A)%*%t(A)  # Vector-> coordinates in (x1e,e2) base.
Beta_0 = M%*%c(0,1,1) #theoretical coefficients
Beta_0

#step 4.1 : controlling for e3 in step 3
reg4.1 = lm(y~e2+x1e+e3)
summary(reg4.1)
#here, allowing a non null coeff on x1e has a cost in the direction of e1.
#step 4.2 : see what happens without this cost
x3e = 1/sqrt(2) * e2 +e3
reg4.2 = lm(y~e2+x3e)
summary(reg4.2)

#step 4.3 : controlling for e1 instead of e3 is even worse than not controlling
reg4.3 = lm(y~e2+x1e+e1)
summary(reg4.3)
#here e1 reinforces the bias by making possible to kill the cost of putting weight on the biased variables

#4.1 shows how controlling can kill included variable bias
#4.2 compared to 3 shows that the existence of a cost of adding a variable may be helpful
#4.3 shows that the wrong control can bias the results even more by killing this cost

#step 5 : more understandable illustration
x2 = e2+e3/3
reg5= lm(y~x2+e2)
summary(reg5)

#step 6 : cases 3 to 5 work  if y is tainted by an eps orthogonal to e1,e2,e3.
# It just lowers R2 but we cannot reduce MSE in eps direction, so no bias is created on the theoretical counterpart
ye = y+eps

reg4.2 = lm(ye~e2+x1e+e3)
summary(reg4.2)

reg5.2= lm(y~x2+e2)
summary(reg5.2)

# More realistic example :

V2 = matrix(rep(1,9), nrow = 3) #variance matrix : base will be orthogonal and normed vectors.
V2 = V2+diag(3)
V2[c(3,7)] = 0
V2=V2/2
V2
X = mvrnorm(50000,c(0,0,0), V2)
eps = mvrnorm(50000, c(0,0,0), diag(3)) #Indept from X

x1 = X[,1]
x2 = X[,2]
x3 = X[,3]
eps1 = eps[,1]
eps2 = eps[,2]
eps3 = eps[,3]


y = x2 + x3 + eps1 + eps2
x3e = x3 + eps1/4

reg = lm(y~x2) # we want the effect of x2 => omitted vriable bias
summary(reg)

reg2 = lm(y~x2+x3)
summary(reg2) # we have cleaned this omitted variable bias by adding x3

#Now suppose we only have an x3 that adds a bias
reg3 = lm(y~x2+x3e) #here, x3 is tainted with eps1. We can reduce MSE by adding weight on x3e
summary(reg3)


# Extra : common measurement errors on variable 

# suppose y is tainted by a measurement error. If this error is independant from explanatory variables, all is fine.

ye =y+eps3
reg_mes_err_orth = lm(ye~x2+x3)
summary(reg_mes_err_orth)

# if for some reason, the measurement error is shared with an explanory variable :
x2e = x2+eps3/2
reg_mes_err_non_orth = lm(ye~x2e+x3)
summary(reg_mes_err_non_orth)
#here x2 coefficient has been boosted to reduce MSE due to eps3, x3 has been reduced to compensate x2 boost since they are positively correlated.

#Difference with omited variable bias
#####

#####
#Controlling for a constant




sum((lm(y1~x1+x2-1)$coefficients)*c(1,2))





#Categorical variable
normR2 = function(x){
  return(c(sqrt(t(x)%*%x)))
}

cat = as.integer(y>0.5)
ryc = lm(y~cat)
mean(y) - ryc$coefficients[2] * mean(cat)

pred = c(y%*%(cat-mean(cat))/normR2(cat-mean(cat))^2) * (cat-mean(cat))+ mean(y)
ryc$fitted.values-pred
normR2(ryc$fitted.values - pred)/50000
