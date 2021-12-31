source("mcmc_routine_expm.r")

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(ind)

trueValues <- c(c(matrix(c(-2.54, -0.56,
                           -2.94,  0.15,
                           -1.10, -0.03,
                           -3.92,  0.21,
                           -2.12,  1.17), ncol=2, byrow=T)),
                c(  -4.59512, -1.15268, -2.751535, -2.090741),
                c( -3.178054, -4.59512))

# Add this second row when doing time-inhomogeneous
# 0.11
# -0.24
# -0.15
# 0.23
# 0.08



init_par = trueValues
# init_par = c(c(matrix(c(-2, 0,
#                         -2, 0,
#                         -2, 0,
#                         -2, 0,
#                         -2, 0), ncol=2, byrow=T)),
#                 rep(-3,4), c( -4.5, -5))

par_index = list( beta=1:10, misclass=11:14, pi_logit=15:16)
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

# The cav data here is time homogeneous
load(paste0("Data/cavData", ind, ".rda"))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, 4, drop=F] #dont need years because time-homogeneous
t = temp_data[,"years"]
steps = 10000
burnin = 5000
n_cores = 8

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("MCMC_local/mcmc_out_", ind, ".rda"))
