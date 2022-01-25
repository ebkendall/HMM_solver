source("mcmc_routine_new.r")

# ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
ind = 10
set.seed(ind)

trueValues <- c(c(matrix(c(-2.29709805,  0.09266760, -0.56262135,
                         -1.17308794, -5.10636947, -0.96162312,
                         -1.71474254, -0.04338819,  0.83882558,
                         -2.08300714,  0.03824367, -2.75345311,
                         -2.42208380,  0.11315485,  1.76897841), ncol=3, byrow=T)),
                         c(  -5.60251814, -0.84455697, -2.56906519, -2.12629033),
                         c( -6.95125291, -7.07504453),
                         c(10, 20, 30, 40),
                         rep(1,4))

init_par = trueValues

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21, mu = 22:25, sigma = 26:29)

betaMat <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)


prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))

load(paste0("DataOut/Continuous/cavData", ind, ".rda"))

temp_data = as.matrix(cavData); rownames(temp_data) = NULL
id = temp_data[,"ptnum"]
y_1 = temp_data[,"state"]
y_2 = temp_data[,"cont_resp"]
x = temp_data[, c("years", "sex"), drop=F]
t = temp_data[,"years"]
steps = 10000
burnin = 5000
n_cores = 16

s_time = Sys.time()

mcmc_out = mcmc_routine(y_1, y_2, x, t, id, init_par, prior_par, par_index,
             steps, burnin, n_cores)

e_time = Sys.time() - s_time; print(e_time)

save(mcmc_out, file = paste0("Model_out/Old/mcmc_out_", ind, ".rda"))
