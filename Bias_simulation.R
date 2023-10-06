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
source("C:/Users/tayoy/Documents/GitHub/OLS-Illustrations/Projections.R")

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
V = diag(2) #variance matrix : basis will be orthogonal and normed vectors.
E = mvrnorm(n=50000, mu = c(0,0), Sigma = V)

e1 = E[,1]
e2 = E[,2]
One = rep(1,50000)

y = 1/sqrt(2)*e1 + 1/sqrt(2)*e2
mean(y)     # means of variables are "almost 0"
mean(e1)
mean(e2)
reg_centered = lm(y~e1+e2-1)  #adding -1 in the formula specifies that
                              #no constant should be added to the regressors
reg_centered

# Regression of non centred variables without constant
yo = y+400
x1o = e1+100
x2o = e2+300
mean(yo)    # variables a no longer centered
mean(x1o)
mean(x2o)
reg_offset = lm(yo~x1o+x2o-1)

reg_offset
sum(reg_offset$coefficients * c(100,300))

# Regression of non centred  variables with a constant

reg_constant = lm(yo~x1o+x2o+One-1) # Remember that ones contains ones
reg_constant






###########################
# Illustrating the use of a constant in 3D

basis <- diag(3)
rownames(basis) <- c("E1", "E2", "1")

vec = matrix(1/sqrt(2)*basis[1,]+sqrt(1-1/sqrt(2)^2)*basis[2,], nrow = 1)
rownames(vec) = "Y"
open3d()
planes3d(0, 0, 1, 0, col="cornflowerblue", alpha=0.4)

vectors3d(basis, col = "black")
vectors3d(vec, col = "red")
segments3d(rbind(c(1/sqrt(2),0,0),c(1/sqrt(2),1/sqrt(2),0)))
segments3d(rbind(c(0,1/sqrt(2),0),c(1/sqrt(2),1/sqrt(2),0)))
light3d(0,0)
close3d()
open3d()

names(basis) = NULL
vecX1O = matrix(basis[1,]+1*basis[3,], nrow = 1)
vecX2O = matrix(basis[2,]+3*basis[3,], nrow = 1)
vecY = matrix(1/sqrt(2)*basis[1,]+1/sqrt(2)*basis[2,]+4*basis[3,], nrow = 1)
vec = rbind(vecX1O, vecX2O, vecY)
rownames(vec) = c("X1o", "X2o", "Yo")

planes3d(0, 0, 1, 0, col="cornflowerblue", alpha=0.2)
vectors3d(basis, col = "black")
vectors3d(vec, col = "blue")
segments3d(rbind(c(1/sqrt(2),0,0),c(1/sqrt(2),1/sqrt(2),0)))
segments3d(rbind(c(0,1/sqrt(2),0),c(1/sqrt(2),1/sqrt(2),0)))
segments3d(rbind(c(0,1,0),c(0,1,3)))
segments3d(rbind(c(1,0,0),c(1,0,1)))
segments3d(rbind(c(1/sqrt(2),1/sqrt(2),0),c(1/sqrt(2),1/sqrt(2),4)))

open3d()

names(basis) = NULL
vecX1O = matrix(basis[1,]+2*basis[3,], nrow = 1)
vecX2O = matrix(basis[2,]+6*basis[3,], nrow = 1)
vecY = matrix(1/sqrt(2)*basis[1,]+1/sqrt(2)*basis[2,]+8*basis[3,], nrow = 1)
vec = rbind(vecX1O, vecX2O, vecY)
rownames(vec) = c("X1o", "X2o", "Yo")

planes3d(0, 0, 1, 0, col="cornflowerblue", alpha=0.2)
vectors3d(basis, col = "black")
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
rownames(vec) <- c("E1", "E2", "E3")
open3d()
vectors3d(vec, color=rep("black",3), lwd=2, frac.lab =1.2, cex.lab = 2)
# draw the XZ plane, whose equation is Y=0
planes3d(0, 0, 1, 0, col="gray", alpha=0.2)


#####
# First regression 
set.seed(18)
V = diag(3) # variance matrix 
e = mvrnorm(n=50000, mu = c(0,0,0), Sigma = V)

e1 = e[,1]
e2 = e[,2]
e3 = e[,3]

# e1,e2,e3 are observations of 3  random variables that form an orthonormal basis
# of a 3 dimension subspace of RV in |R

# 1st step : generating y
y = e3+e2
reg1 = lm(y~e2-1) # -1 specifies not to add a constant to the model.
reg1
R2(reg1)          # ~0.5 as expected : ||e2||^2/||e2+e3||^2

# 3D representation of reg1

basis <- diag(3)
rownames(basis) <- c("E1", "E2", "E3")

Y = matrix(basis[2,]+ basis[3,], nrow = 1)
rownames(Y) = "Y"
E2 = matrix(basis[2,], nrow = 1)
rownames(E2) = "E2"
open3d()
planes3d(0, 0, 1, 0, col="lightgrey")
vectors3d(basis[c(1,3),], col = "black",  lwd = 2)
vectors3d(Y, col ="red",  lwd = 2, frac.lab = 1.05)
vectors3d(E2, col = "blue4",  lwd = 2, frac.lab = 1.2)
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


# 3D representation of reg2
basis <- diag(3)
rownames(basis) <- c("E1", "E2", "E3")
X1 = matrix(1/sqrt(2)*basis[1,]+1/sqrt(2)*basis[2,], nrow = 1)
rownames(X1) = "X1"

open3d()
planes3d(0, 0, 1, 0, col="cornflowerblue")
 
vectors3d(basis[c(1,3),], col = "black",  lwd = 2)
vectors3d(Y, col ="red",  lwd = 2, frac.lab = 1.05)
vectors3d(E2, col = "blue4",  lwd = 2, frac.lab = 1.2)
vectors3d(X1, col = "blue4", lwd = 2)
segments3d(rbind(c(0,1,0),c(0,1,1)), col = "red")
corner(c(0,0,0),c(0,1,0),c(0,1,1), col = "red")

light3d(0, 0)
light3d(0, 0)
light3d(0, 0)
close3d()
# 3D representation of Reg2 if  coefficient of X1 is not 0:
basis <- diag(3)
rownames(basis) <- c("E1", "E2", "E3")
X1offset =as.matrix(X1/3+E2, nrow = 1)
rownames(X1offset) = "X1"

open3d()
planes3d(0, 0, 1, 0, col="cornflowerblue")

vectors3d(basis[c(1,3),], col = "black",  lwd = 2)
vectors3d(Y, col ="red",  lwd = 2, frac.lab = 1.05)
vectors3d(E2, col = "blue4",  lwd = 2, frac.lab = 1.2)
vectors3d(X1offset, col = "blue4", lwd = 2, origin = E2, frac.lab = 0.5)
segments3d(rbind(c(0,1,0),c(0,1,1)), col = "red", lwd = 2) # eta
corner(c(0,0,0),c(0,1,0),c(0,1,1), col = "red")
segments3d(rbind(X1offset,c(0,1,1)), col = "darkorange2", lwd = 2) #eta'

light3d(0, 0)
light3d(0, 0)
light3d(0, 0)
close3d()


# 3rd step : A control orthogonal to e2 cannot change its coefficient.

xo = 1/sqrt(2)*e1+1/sqrt(2)*e3

cor(xo,e2)  # "Almost" 0

reg3 = lm(y~e2+xo-1)
reg3
R2(reg3)
# 3D representation of reg3
vec =rbind( 1/sqrt(2)*c(1,0,1), c(0,1,0)) %*% basis
rownames(vec) = c("Xo", "E2")
pointfall = t(vec)%*%project(t(vec))%*% matrix(c(0,1,1), nrow = 3)
open3d()
planes3d(0, 0, 1, 0, col="lightgrey")
planes3d(1, 0,-1, 0, col = "cornflowerblue")
vectors3d(basis[c(1,3),], col = "black",  lwd = 2)
vectors3d(vec, col = "blue4", lwd = 2)
vectors3d(Y, col = "red", lwd = 2, frac.lab = 1.05)
segments3d(rbind(c(0,0,0),t(pointfall)), col = "red")
segments3d(rbind(t(pointfall),c(0,1,1)), col = "red")
segments3d(rbind(t(pointfall),c(0,1,0)), col = "black")
corner(c(0,0,0),as.vector(pointfall),c(0,1,1), col = "red")
corner(c(0,0,0),c(0,1,0),as.vector(pointfall), col = "black")

light3d(0,0)
light3d(0,0)

open3d()
planes3d(0, 0, 1, 0, col="lightgrey")

vectors3d(basis[c(1,3),], col = "black",  lwd = 2)
vectors3d(vec, col = "blue4", lwd = 2)
vectors3d(Y, col = "red", lwd = 2, frac.lab = 1.05)

corner(c(0,0,0),c(0,1,0),as.vector(Y), col = "red")
segments3d(rbind(c(0,1,0), as.vector(Y)), col = "red")

light3d(0,0)
light3d(0,0)

close3d()
close3d()


# 4th step : adding a control in vect(e2, eta) induces bias
xc = e3+1/2*e2
reg3 = lm(y~e2+xc-1)
reg3
R2(reg3)
# 3D representation of extra control
eta = matrix(c(0,0,1), nrow = 1)
rownames(eta) = "eta"
D = matrix(c(0,1,0), nrow = 1)
rownames(D) = "D"
X = matrix(c(1/sqrt(3),1/sqrt(3),1/sqrt(3),-1/sqrt(3),1/sqrt(3),1/sqrt(3)), nrow = 2, byrow = TRUE)
rownames(X) = c("X1","X2")
Xc = matrix(c(0,1/sqrt(2), 1/sqrt(2)), nrow = 1)
rownames(Xc) = "Xc"

open3d()
vectors3d(D, lwd = 3, frac.lab = 1.1, col = "dodgerblue")
vectors3d(X, lwd = 3, frac.lab = 0.5, col = "black")
vectors3d(eta, lwd = 3, frac.lab = 0.5, col = "red")
vectors3d(Xc, lwd = 3, frac.lab = 0.8, col = "green3")
planes3d(0,0,1,0, col = "lightgrey")
planes3d(0,1,-1,0, col = "lemonchiffon", alpha = 0.6)
text3d(0,0.6,0.8, "Vect(X1,X2)", col = "lemonchiffon4", font = 2)
light3d(0,0)

close3d()
# 3D representation of reg4
vec =rbind(c(0,0.5,1), c(0,1,0)) %*% basis
rownames(vec) = c("Xc", "e2")

userMat = matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1),byrow = TRUE, nrow = 4) # view from front
userMat2= matrix(c(0,1,0,0,-0.5,0,1,0,1,0,0.5,0,0,0,0,1),byrow = TRUE, nrow = 4) #View from a little above
open3d()
view3d(userMatrix = userMat2)
planes3d(0, 0, 1, 0, col="lightgrey")
planes3d(1, 0, 0, 0, col = "cornflowerblue")
vectors3d(basis[c(1,3),], col = "black",  lwd = 2)
vectors3d(vec, col = "blue4", lwd = 2, frac.lab = 1.1)
vectors3d(Y, col = "red", lwd = 2)


# more cases of reg4


for (i in 1:6){
  open3d()
  
  vecxc =matrix(c(0,0.3*i,1),nrow= 1)
  
  rownames(vecxc) = "Xc"
  vectors3d(basis,lwd=2)
  segments3d(rbind(c(0,1,1),  c(0.0001,0.3*i,1)), col= "blue4", lwd=2)
  if(i<=1/0.3){
    cone3d(tip = c(0.001,0.1,0), base= c(0,0.9,1) ,length = 0.05, radius= 0.02, col = "blue4" )

  }
  else{
    cone3d(tip = c(0.001,-0.1,0), base= c(0,1.1,1),length = 0.05, radius= 0.02, col = "blue4")

  }
  vectors3d(vecxc, col = "blue4",lwd=2, frac.lab = 0.5)
  vectors3d(Y, col = "red", lwd = 2)
  planes3d(1, 0, 0, 0, col = "cornflowerblue")
  planes3d(0,0,1,0, col = "lightgrey")
  view3d(userMatrix = userMat2)

  filename = paste("C:/Users/tayoy/Documents/ENSAE/2A/RBB - Regressions/reg4.",i,".png", sep = "")
  rgl.snapshot(filename, fmt = 'png')
}

close3d()
close3d()
close3d()
close3d()
close3d()
close3d()
close3d()

# Regression 4.7 if beta_x is bounded
open3d()
vecxc = matrix(c(0,1 ,05/18),nrow= 1)
vece2 = matrix(c(0,0.5,0),nrow= 1)
veceta= matrix(c(0,1,1), nrow = 1)
vecy = matrix(c(0,1,1), nrow = 1)
rownames(vecy) = "Y"
rownames(vecxc) = "Xc"
rownames(vece2) = "E2"
rownames(veceta) = "eta"
vectors3d(basis[3,],lwd=2)
vectors3d(vece2,color = "blue4", lwd = 2, frac.lab = 0.5)
vectors3d(vecxc, origin = vece2, color = "blue4", lwd = 2,frac.lab = 0.5)
vectors3d(veceta, origin = vecxc, color = "black", lwd = 2,frac.lab = 0.5)
vectors3d(vecy, col = "red", lwd = 2,frac.lab = 0.5)
planes3d(1, 0, 0, 0.0001, col = "cornflowerblue")
planes3d(0,0,1,0, col = "lightgrey")
view3d(userMatrix = userMat)

filename = paste("C:/Users/tayoy/Documents/ENSAE/2A/RBB - Regressions/reg4.",7,".png", sep = "")
rgl.snapshot(filename, fmt = 'png')
close3d()




######
#Second part : mixed cases

set.seed(18)
V = diag(3) # variance matrix 
E = mvrnorm(n=50000, mu = c(0,0,0), Sigma = V)
e1 = E[,1]
e2 = E[,2]
e3 = E[,3]
y = e3+e2
d = e2
g = e3
x_intermediate = 1/sqrt(3)*e1 + 1/sqrt(3)*e2 + 1/sqrt(3)*e3

#Propagation of bias: 

reg5 = lm(y~d+x_intermediate-1)
reg5
R2(reg5)

#Theoretical counterpart :
A = matrix( c(1/sqrt(3),1/sqrt(3),1/sqrt(3),0,1,0), nrow = 3) 
M = matrix.inverse(t(A)%*%A)%*%t(A)  # Vector-> coordinates in (x1e,e2) basis.
gamma_0 = M%*%c(0,1,1) #theoretical coefficients
rownames(gamma_0) = c("X1e", "e2")
gamma_0

reg5.1 = lm(y~d + x_intermediate + g -1)
reg5.1
