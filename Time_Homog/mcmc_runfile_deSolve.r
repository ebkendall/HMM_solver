source("mcmc_routine.r")

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(ind)

init_par = trueValues = c(c(matrix(c(-2.26568339, -0.49991746,
                                     -1.22022878, -0.82779213,
                                     -1.56180104,  0.73838829,
                                     -2.20978996, -1.83682627,
                                     -2.41222255,  1.63135439), ncol=2, byrow=T)),
                  c(  -5.73343061, -0.78623894, -2.52747176, -2.12144526),
                  c( -6.52842355, -6.15970066))

par_index = list( beta=1:10, misclass=11:14, pi_logit=15:16)
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

# The cav data here is time homogeneous
load(paste0("Data/cavData", ind, ".rda"))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, "sex", drop=F] #dont need years because time-homogeneous
t = temp_data[,"years"]
steps = 20000
burnin = 5000
n_cores = 16
disc = F

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, disc)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/deSolve/mcmc_out_", ind, ".rda"))
