source("mcmc_routine.r")

ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(ind)

init_par = trueValues = c(c(matrix(c(-2.26568339, 3 *  0.08766060, -0.49991746,
                           -1.22022878, 3 * -4.44888558, -0.82779213,
                           -1.56180104, 3 * -0.08262607,  0.73838829,
                           -2.20978996, 3 *  0.05404948, -1.83682627,
                           -2.41222255, 3 *  0.10833734,  1.63135439), ncol=3, byrow=T)),
                  c(  -5.73343061, -0.78623894, -2.52747176, -2.12144526),
                  c( -6.52842355, -6.15970066))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21)

prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

data_files <- c("Year/", "YearTwo/", "Month/")


folder = 2

Dir <- paste0("../DataGeneration/DataOut/", data_files[folder])
load(paste0(Dir,'cavData',ind,'.rda'))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y = temp_data[,"state"]
x = temp_data[, c("disc_time", "sex"), drop=F]
t = temp_data[,"years"]
steps = 20000
burnin = 5000
n_cores = 16
disc = T

s_time = Sys.time()

mcmc_out = mcmc_routine(y, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores, disc)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0('Model_out/expm/', data_files[folder] ,'mcmc_out_',ind,'.rda'))
