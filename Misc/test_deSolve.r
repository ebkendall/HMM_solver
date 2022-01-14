# **************************************************************************
# This works for representing time inhomogeneous as well as time homogeneous
# **************************************************************************

library(deSolve)
library(expm)

Q <- function(iyears,sex,betaMat){
  
  q1  = exp( c(1,iyears,sex) %*% betaMat[1,] )  # Transition from state 1 to state 2.
  q2  = exp( c(1,iyears,sex) %*% betaMat[2,] )  # Transition from state 2 to state 3.
  q3  = exp( c(1,iyears,sex) %*% betaMat[3,] )  # Transition from state 1 to death.
  q4  = exp( c(1,iyears,sex) %*% betaMat[4,] )  # Transition from state 2 to death.
  q5  = exp( c(1,iyears,sex) %*% betaMat[5,] )  # Transition from state 3 to death.
  
  qmat = matrix(c( 0,q1, 0,q2,
                   0, 0,q3,q4,
                   0, 0, 0,q5,
                   0, 0, 0, 0),nrow=4,byrow=TRUE) 
  diag(qmat) = -rowSums(qmat)
  
  return(qmat)
}

out_mat <- function(t, out) {
  myMat <- matrix(c(out[t,"p1"], out[t,"p2"], out[t,"p3"], out[t,"p4"],
                           0, out[t,"p5"], out[t,"p6"], out[t,"p7"],
                           0,  0, out[t,"p8"], out[t,"p9"],
                           0,  0,  0,  1), nrow = 4, byrow = T)
  return(myMat)
}


#out_mat = Vectorize(out_mat, vectorize.args = "t", SIMPLIFY = T)

# Model dependent on time ----------------------------------------------

model_t <- function(t,p,b) {
  
  betaMat <- matrix(b, ncol = 3, byrow = F)
  
  q1  = exp( c(1,t,1) %*% betaMat[1,] )  # Transition from state 1 to state 2.
  q2  = exp( c(1,t,1) %*% betaMat[2,] )  # Transition from state 2 to state 3.
  q3  = exp( c(1,t,1) %*% betaMat[3,] )  # Transition from state 1 to death.
  q4  = exp( c(1,t,1) %*% betaMat[4,] )  # Transition from state 2 to death.
  q5  = exp( c(1,t,1) %*% betaMat[5,] )  # Transition from state 3 to death.
  
  dP = rep(1,9) # this is the vector with all diffEqs
  
  dP[1] = p[1]*(-q1-q2)
  dP[2] = p[1]*q1 + p[2]*(-q3-q4)
  dP[3] = p[2]*q3 - p[3]*q5
  dP[4] = p[1]*q2 + p[2]*q4 + p[3]*q5
  dP[5] = p[5]*(-q3-q4)
  dP[6] = p[5]*q3 - p[6]*q5
  dP[7] = p[5]*q4 + p[6]*q5
  dP[8] = -p[8]*q5
  dP[9] = p[8]*q5
  
  return(list(dP)) 
  
}

trueValues <- c(c(matrix(c(-2.54, 0.11, -0.56,
                           -2.94, -0.24,  0.15,
                           -1.10, -0.15, -0.03,
                           -3.92, 0.23,  0.21,
                           -2.12, 0.08,  1.17), ncol=3, byrow=T)),
                c(  -4.59512, -1.15268, -2.751535, -2.090741),
                c( -3.178054, -4.59512))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

parms_t <- trueValues[par_index$beta]
p <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=1,p9=0)
times <- c(5,10)

out_time <- ode(p, times, model_t, parms = parms_t)


t_int <- seq(5,10, 0.001)
P = diag(4)
for(i in t_int) {
  
  tempQ <- Q(i, 1, betaMat_t)
  P = P %*% expm(0.001* tempQ)
  
}


