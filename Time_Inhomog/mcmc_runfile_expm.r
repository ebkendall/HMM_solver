source("mcmc_routine_expm.r")

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(ind)

trueValues <- c(c(matrix(c(-2.54,  0.11, -0.56,
                           -2.94, -0.24,  0.15,
                           -1.10, -0.15, -0.03,
                           -3.92,  0.23,  0.21,
                           -2.12,  0.08,  1.17), ncol=3, byrow=T)),
                c(  -4.59512, -1.15268, -2.751535, -2.090741),
                c( -3.178054, -4.59512))

init_par = trueValues

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

betaMat <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

data_files <- c("Year", "YearTwo", "Month")

for (folder in 1:3) {

    Dir <- paste0("../DataGeneration/DataOut/", data_files[folder], "/")
    load(paste0(Dir,'cavData',seedInd,'.rda'))

    temp_data = as.matrix(cavData); rownames(temp_data) = NULL
    id = temp_data[,"ptnum"]
    y = temp_data[,"state"]
    x = temp_data[, c("disc_time", "sex"), drop=F]
    t = temp_data[,"years"]
    steps = 20000
    burnin = 5000
    n_cores = 16

    s_time = Sys.time()

    mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
                 steps, burnin, n_cores)

    e_time = Sys.time() - s_time; print(e_time)

    save(mcmc_out, file = paste0('Model_out/expm/', data_files[folder] ,'/mcmc_out_',seedInd,'.rda'))

}
