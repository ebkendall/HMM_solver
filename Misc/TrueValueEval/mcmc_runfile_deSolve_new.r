source("mcmc_routine_deSolve.r")

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(ind)
init_par <- c(c(matrix(c(-2.54, 0.11, -0.56,
                           -2.94, -0.24,  0.15,
                           -1.10, -0.15, -0.03,
                           -3.92, 0.23,  0.21,
                           -2.12, 0.08,  1.17), ncol=3, byrow=T)),
                c(  -4.59512, -1.15268, -2.751535, -2.090741),
                c( -3.178054, -4.59512))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

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
steps = 20000
burnin = 5000
n_cores = 16

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/New/mcmc_out_", ind, ".rda"))
