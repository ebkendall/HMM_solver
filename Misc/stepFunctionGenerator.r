# Let dis1 be the true discretization, and dis2 is the estimated on
stepFunction <- function(dis1, dis2, x_low, x_up) {
  
  set.seed(10)
  
  #creating xlim based on the larger discritization
  maxDis= max(c(dis1, dis2))
  xlimit = 1
  
  dis1seq = seq(0, xlimit, by = dis1)
  dis2seq = seq(0, xlimit, by = dis2)
  
  x = y1 = y2 = runif(15, min = 0, max = xlimit)
  
  for(i in 1:length(y1)) {
    indDis = max(which(dis1seq <= x[i])) #selects which month to go to
    y1[i] = dis1seq[indDis]
  }
  
  for(i in 1:length(y2)) {
    indDis = max(which(dis2seq <= x[i])) #selects which month to go to
    y2[i] = dis2seq[indDis]
  }
  
  mod1 = lm(y1 ~ x)
  mod2 = lm(y2 ~ x)
  
  xlimit = 1 #ceiling(10 * maxDis)
  
  x1_o = y1_o = dis1seq
  x2_o = y2_o = dis2seq
  
  ylimit = c(min(coefficients(mod1)[1], coefficients(mod2)[1]), 1)
  
  
  plot(x1_o, y1_o, xlim = c(0,xlimit), ylim = ylimit, xaxt='n', yaxt='n', bty="n",
       xlab = "Time", ylab = "Discretization of Time", col = "white")
  points(x2_o, y2_o, col = "white")
  abline(h=0)
  abline(v=0)
  
  plot(stepfun(x1_o, c(0,y1_o)), col = "red", add = T, do.points = F)
  plot(stepfun(x2_o, c(0,y2_o)), col = "blue", add = T, do.points = F)
  
  abline(mod1, col = "red")
  abline(mod2, col = "blue")
  
  points(x, y1, col = "red")
  points(x, y2, col = "blue")
  
  print(paste0("Intercept of red curve is: ", coefficients(mod1)[1], " with slope ", coefficients(mod1)[2]))
  print(paste0("Intercept of blue curve is: ", coefficients(mod2)[1], " with slope ", coefficients(mod2)[2]))
  
}

# ---------------------------------------------------------------
# - NEW STEP FUNCTION THAT IS MINIMIZING THE SQUARED ERROR LOSS -
# ---------------------------------------------------------------

t_i = seq(-1000,999.99, by = 0.01)

y_i = floor(t_i) 

X = cbind(1, t_i)

b = solve(t(X) %*% X) %*% t(X) %*% y_i

Xhat = cbind(1, seq(-1000, 1000, by=0.01))
yhat = Xhat %*% b


dis1seq = seq(-1000, 1000, by = 0.5)
t_i_half = t_i
for(i in 1:length(t_i)) {
  indDis = max(which(dis1seq <= t_i[i])) #selects which month to go to
  t_i_half[i] = dis1seq[indDis]
}

X_half = cbind(1, t_i_half)

b_half = solve(t(X_half) %*% X_half) %*% t(X_half) %*% y_i


yhat_half = X_half %*% b_half

#plot(stepfun(X_half[,2], c(min(yhat_half), yhat_half)), col = "red", do.points = F)

# ------------ PLOTTING ------------
plot(x = NULL, y = NULL, ylim = c(-1,5), xlim = c(-1,5))
lines(t_i, y_i)
abline(h = 0, lty = 2); abline(v=0, lty = 2)

lines(Xhat[,2],yhat, col= "blue")

lines(X_half[,2],yhat_half, col = "red")

points(t_i_half, y_i, col = "red")

legend("topleft", inset = 0.05, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8,
       title="Line types", text.font=4, bg='lightblue')


# Testing the heights
temp = c(unique(yhat_half))
seq1 = kronecker(c(-1000:999),c(1,1))
temp = temp - seq1; temp
sum((temp^2))

# Note: it is minimizing over the entire region.