t_i = seq(-1000,999.99, by = 0.01)

# y_i = floor(t_i) 
y_i = t_i

X_hat = cbind(1, t_i - (t_i %% 1))
b = solve(t(X_hat) %*% X_hat) %*% t(X_hat) %*% y_i
y_hat = X_hat %*% b

# X = cbind(1, t_i)
# 
# b = solve(t(X) %*% X) %*% t(X) %*% y_i
# 
# Xhat = cbind(1, seq(-1000, 1000, by=0.01))
# yhat = Xhat %*% b


# dis1seq = seq(-1000, 1000, by = 0.5)
# t_i_half = t_i
# for(i in 1:length(t_i)) {
#   indDis = max(which(dis1seq <= t_i[i])) #selects which month to go to
#   t_i_half[i] = dis1seq[indDis]
# }
# 
# X_half = cbind(1, t_i_half)
# 
# b_half = solve(t(X_half) %*% X_half) %*% t(X_half) %*% y_i
# 
# 
# yhat_half = X_half %*% b_half

# ------------ PLOTTING ------------
plot(x = NULL, y = NULL, ylim = c(-1,5), xlim = c(-1,5))
lines(t_i, y_i)
abline(h = 0, lty = 2); abline(v=0, lty = 2)

temp = stepfun(t_i[-1], y_hat)
plot(temp, verticals = F, col = "red", add = T)
# lines(X_half[,2],yhat_half, col = "red")

# points(t_i_half, y_i, col = "red")

legend("topleft", inset = 0.05, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8,
       title="Line types", text.font=4, bg='lightblue')


# Testing the heights
temp = c(unique(yhat_half))
seq1 = kronecker(c(-1000:999),c(1,1))
temp = temp - seq1; temp
sum((temp^2))

# Note: it is minimizing over the entire region.


# ---------------------------------------------------------------
# - NEW STEP FUNCTION THAT IS MINIMIZING THE SQUARED ERROR LOSS -
# ---------------------------------------------------------------
final_plot_fnc <- function(d1, d2, d3) {

  t_i = seq(-1000,999.99, by = 0.01)
  y_i = t_i
  
  X_hat_1 = cbind(1, t_i - (t_i %% d1))
  b1 = solve(t(X_hat_1) %*% X_hat_1) %*% t(X_hat_1) %*% y_i
  y_hat_1 = X_hat_1 %*% b1
  
  X_hat_2 = cbind(1, t_i - (t_i %% d2))
  b2 = solve(t(X_hat_2) %*% X_hat_2) %*% t(X_hat_2) %*% y_i
  y_hat_2 = X_hat_2 %*% b2
  
  X_hat_3 = cbind(1, t_i - (t_i %% d3))
  b3 = solve(t(X_hat_3) %*% X_hat_3) %*% t(X_hat_3) %*% y_i
  y_hat_3 = X_hat_3 %*% b3
  
  # ------------ Expectations ------------
  W = cbind(1, t_i)
  beta = matrix(c(0,1), nrow = 2)
  e1 = solve(t(X_hat_1) %*% X_hat_1) %*% t(X_hat_1) %*% W %*% beta
  e2 = solve(t(X_hat_2) %*% X_hat_2) %*% t(X_hat_2) %*% W %*% beta
  e3 = solve(t(X_hat_3) %*% X_hat_3) %*% t(X_hat_3) %*% W %*% beta
  
  e1 = round(e1, digits = 5)
  e2 = round(e2, digits = 5)
  e3 = round(e3, digits = 5)
  
  # ------------ PLOTTING ------------
  plot(x = NULL, y = NULL, ylim = c(-1,5), xlim = c(-1,5),
       xlab = "t_i", ylab = "Y_i", 
       main = "Illustration of Theorem 1")
  lines(t_i, y_i)
  abline(h = 0, lty = 2); abline(v=0, lty = 2)
  
  temp1 = stepfun(t_i[-1], y_hat_1)
  temp2 = stepfun(t_i[-1], y_hat_2)
  temp3 = stepfun(t_i[-1], y_hat_3)
  
  plot(temp1, verticals = F, col = "red", add = T)
  plot(temp2, verticals = F, col = "blue", add = T)
  plot(temp3, verticals = F, col = "purple", add = T)
  
  legend("topleft", inset = 0.05, 
         legend=c(paste0(d1, " unit(s)"), 
                  paste0(d2, " unit(s)"),
                  paste0(d3, " unit(s)")),
         col=c("red", "blue", "purple"), lty=1, cex=0.8,
         title="Discretization Type", text.font=4, bg='lightblue')
  
  legend("bottomright", inset = 0.05,
         legend=c(paste0("b0: ", formatC(e1[1,1],format = "f", digits=5), 
                         ", b1: ", formatC(e1[2,1],format = "f", digits=5)), 
                  paste0("b0: ", formatC(e2[1,1],format = "f", digits=5),
                         ", b1: ", formatC(e2[2,1],format = "f", digits=5)),
                  paste0("b0: ", formatC(e3[1,1],format = "f", digits=5),
                         ", b1: ", formatC(e3[2,1],format = "f", digits=5))),
         col=c("red", "blue", "purple"), lty=1, cex=0.8,
         title="Expected Values", text.font=4, bg='lightblue')

  
}
