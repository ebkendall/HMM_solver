source("mcmc_routine.r")

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(ind)

trialNum = 1

init_par = trueValues = c(c(matrix(c(-3, 0, 0,
                                     -3, 0, 0,
                                     -3, 0, 0,
                                     -3, 0, 0,
                                     -3, 0, 0), ncol=3, byrow=T)),
                  c(  -5, -5, -5, -5),
                  c( -5, -5))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load(paste0("Data/cavData_cont.rda"))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, c("disc_time", "sex"), drop=F]
t = temp_data[,"years"]
steps = 30000
burnin = 5000
n_cores = 16
disc = F

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, disc)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/deSolve/mcmc_out_", ind, "_", trialNum, ".rda"))
