
trial  <- function(n=2000, dep = 1000, len = 1000) {
  
  error <- rnorm(n, mean = 0, sd = 1)
  
  
  cartel <- c(rep(0, n/2), rep(1,n/2))
  char <- rep(0,n)
  char[dep:(dep+len)] = 1
  
  
  price <- 10*cartel + 5*char -5*(cartel*char)+ error
  
  df <- list(price = price, cartel = cartel, char = char)
  
  return(df)
}

df = trial(n=2000, dep = 1, len = 1099)


plotter<-function(df,lmod = NULL){
  price_regular = df$price
  price_regular[df$char ==1] = NA
  price_char = df$price
  price_char[df$char ==0] = NA
  
  plot(price_regular, ylim = c(-10,30))
  points(price_char, col="red")
  abline(v = 1000, lty = 2)
  abline(h = 0, lty = 2)
  if (is.null(lmod)){
    return()
  }
  M = rbind (rep(1,2000), df$cartel, df$char)
  for (i in 1:3){
    M[i,] = M[i,] * lmod$coefficients[i]
  }
  lines(M[1,], col = "black", lty = 2)
  lines(M[2,], col = "orange")
  lines(M[3,], col = "blue")
  lines(lmod$fitted.values, col = "green", lwd = 2)
}

plotter(df)

lmod = lm(price~cartel+char, data = df)
m = matrix(0,3,0)
for (i in 1:901){
  df = trial(n=2000, dep = i, len = 1099)
  m = cbind(m, as.matrix(lm(price~cartel+char, data = df)$coefficients, nrow = 3))
  if (any(c(1,500,901) == i)){
    dev.new()
    plotter(df, lm(price~cartel+char, data = df))
  }
}
plot(m[1,]-100, ylim =c(-10, 100), col = "green")
points(m[2,], col = "red")
points(m[3,], col = "blue")


dev.off()

# Repeat
n_trials = 10000 # 10 thousand trials
estimates <- vector("numeric", length = n_trials)
