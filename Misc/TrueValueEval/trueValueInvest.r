source("mcmc_routine_deSolve.r")

load("Model_out/New/mcmc_out_1.rda")

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)
      
prior_par = data.frame( prior_mean=rep( 0, ncol(mcmc_out$chain)),
                        prior_sd=rep( 20, ncol(mcmc_out$chain)))

temp = msm::cav
cavData <- data.frame(  "ptnum" = temp$PTNUM,
                        "years" = temp$years,
                        "disc_time" = temp$years,
                        "sex" = temp$sex, 
                        "state" = temp$state,
                        "obstrue" = rep(0, nrow(temp)))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"] 
x = temp_data[, c("years", "sex"), drop=F]
t = temp_data[,"years"]

max_like = fn_log_post(mcmc_out$chain[1,], prior_par, par_index, x, y, t, id) 
max_ind = 1
print(max_like)

for(i in 2:nrow(mcmc_out$chain)) {
	print(i)
	temp = fn_log_post(mcmc_out$chain[i,], prior_par, par_index, x, y, t, id)	
	
	if(temp > max_like) {
		max_ind = i
		max_like = temp
		print(max_like)
	}

}
print(paste0("Row of MCMC chain: ", max_ind))
print(paste0("Likelihood evaluation: ", max_like))
print("MLE parameters:")
print(mcmc_out$chain[max_ind,])
