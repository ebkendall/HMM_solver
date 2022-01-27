# This script file produces trace plots and histograms of the mcmc output files
# To run from command line...
# $ Rscript out_file_mcmc.r num_seeds directory simulation
requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

args = commandArgs(TRUE)
# num_seeds = 7 #args[1] # Only used in this script to locate and load files
folder = as.numeric(args[1])

model_name = c('deSolve', 'expm')
dir = paste0('Model_out/', model_name[folder], '/')
simulation = T #as.logical(args[3]) # true or false

#data_names <- c(rep('orig',3), 'half', 'unhidden')

# Size of posterior sample from mcmc chains
n_post = 5000
# Step number at 3ich the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps =  20000
# Matrix row indices for the posterior sample to use for GFF computation
index_post = (steps - burnin - n_post + 1):(steps - burnin)

par_index = list(beta=1:10, misclass=11:14, pi_logit=15:16)
# par_index = list( beta=1:15, pi_logit=16:17)

true_par = c(c(matrix(c(-2.26568339, -0.49991746,
                        -1.22022878, -0.82779213,
                        -1.56180104,  0.73838829,
                        -2.20978996, -1.83682627,
                        -2.41222255,  1.63135439), ncol=2, byrow=T)),
                  c(  -5.73343061, -0.78623894, -2.52747176, -2.12144526),
                  c( -6.52842355, -6.15970066))

true_par_inv = true_par

true_par_inv[par_index$pi_logit] =
  exp(true_par_inv[par_index$pi_logit])/(1 + exp(true_par_inv[par_index$pi_logit]))
true_par_inv[par_index$misclass[1]] =
  exp(true_par_inv[par_index$misclass[1]])/(1 + exp(true_par_inv[par_index$misclass[1]]))
true_par_inv[par_index$misclass[2:3]] =
  exp(true_par_inv[par_index$misclass[2:3]])/sum(c(1, exp(true_par_inv[par_index$misclass[2:3]])))
true_par_inv[par_index$misclass[4]] =
  exp(true_par_inv[par_index$misclass[4]])/(1 + exp(true_par_inv[par_index$misclass[4]]))

true_par_inv = round(true_par_inv, digits = 6)

labels <- c('b.l. S1 (well)   --->   S2 (mild)',
            'b.l. S1 (well)   --->   S4 (dead)',
            'b.l. S2 (mild)   --->   S3 (severe)',
            'b.l. S2 (mild)   --->   S4 (dead)',
            'b.l. S3 (severe)   --->   S4 (dead)',
            'sex S1 (well)   --->   S2 (mild)',
            'sex S1 (well)   --->   S4 (dead)',
            'sex S2 (mild)   --->   S3 (severe)',
            'sex S2 (mild)   --->   S4 (dead)',
            'sex S3 (severe)   --->   S4 (dead)',
            'P( obs. S2 | true S1 )',
            'P( obs. S1 | true S2 )',
            'P( obs. S3 | true S2 )',
            'P( obs. S2 | true S3 )',
            'P( init S2 )','P( init S3 )')

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Calculating Credible Sets ---------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

cred_set = vector(mode = 'list', length = length(true_par))
for(i in 1:length(cred_set)) cred_set[[i]] = data.frame('lower' = c(-1),
                                                        'upper' = c(-1))
ind = 0
index_seeds = 1:100

for (i in index_seeds) {
    file_name = paste0(dir,'mcmc_out_',toString(i),'.rda')
    if(file.exists(file_name)){
	    load(file_name) # changed toString to 3
        ind = ind + 1

        for(j in 1:length(true_par)) {
            cred_set[[j]][ind,1] =  round(quantile( mcmc_out$chain[index_post,j],
                                        prob=.025), 4)
            cred_set[[j]][ind,2] =  round(quantile( mcmc_out$chain[index_post,j],
                                        prob=.975), 4)
        }
  }
}

save(cred_set, file = paste0('Plots/cred_set_', model_name[folder], '.rda'))

# -----------------------------------------------------------------------------
# Calculating Coverage --------------------------------------------------------
# -----------------------------------------------------------------------------
cov_df <- c()
for(i in 1:length(true_par)) {
    val = true_par[i]
    top = length(which(cred_set[[i]]$lower <= val & val <= cred_set[[i]]$upper))
    bot = nrow(cred_set[[i]])
    covrg = top/bot
    cov_df[i] = covrg
    print(paste0("Coverage for parameter ", val, " is: ", covrg))
}

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------


post_means = matrix(nrow = length(index_seeds), ncol = length(labels))
chain_list <- NULL
ind = 0
for(seed in index_seeds){
  file_name = paste0(dir,'mcmc_out_',toString(seed),'.rda')
	if(file.exists(file_name)){
	    load(file_name) # changed toString to 3
        ind = ind + 1
        print(mcmc_out$accept)

      	chain_list[[ind]] = mcmc_out$chain[index_post,]

        # Putting the violin plots in terms of the probabilities
        mcmc_out$chain[,par_index$pi_logit] =
            exp(mcmc_out$chain[,par_index$pi_logit])/(1 + exp(mcmc_out$chain[,par_index$pi_logit]))
        mcmc_out$chain[,par_index$misclass[1]] =
            exp(mcmc_out$chain[,par_index$misclass[1]])/(1 + exp(mcmc_out$chain[,par_index$misclass[1]]))
        mcmc_out$chain[,par_index$misclass[2:3]] = exp(mcmc_out$chain[,par_index$misclass[2:3]])/(1 +
                                                   exp(mcmc_out$chain[,par_index$misclass[2]]) +
                                                   exp(mcmc_out$chain[,par_index$misclass[3]]))
        mcmc_out$chain[,par_index$misclass[4]] =
            exp(mcmc_out$chain[,par_index$misclass[4]])/(1 + exp(mcmc_out$chain[,par_index$misclass[4]]))

    	post_means[ind,] <- colMeans(mcmc_out$chain[index_post,])
  }
}

# Plot and save the mcmc trace plots and histograms.
stacked_chains = do.call( rbind, chain_list)
par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))
pdf(paste0('Plots/mcmc_', model_name[folder], '.pdf'))
par(mfrow=c(4, 2))
VP <- vector(mode="list", length = length(labels))

for(r in 1:length(labels)){

	plot( NULL, xlab=NA, ylab=NA, main=labels[r], xlim=c(1,length(index_post)),
		    ylim=range(stacked_chains[,r]) )

	for(seed in 1:length(chain_list)) lines( chain_list[[seed]][,r], type='l', col=seed)

	par_mean[r] = round( mean(stacked_chains[,r]), 4)
	par_median[r] = round( median(stacked_chains[,r]), 4)
	upper[r] = round( quantile( stacked_chains[,r], prob=.975), 4)
	lower[r] = round( quantile( stacked_chains[,r], prob=.025), 4)

	hist( stacked_chains[,r], breaks=sqrt(nrow(stacked_chains)), ylab=NA, main=NA,
	      freq=F, xlab=paste0('Mean = ',toString(par_mean[r]),
				                    ' Median = ',toString(par_median[r])))
	abline( v=upper[r], col='red', lwd=2, lty=2)
	abline( v=true_par[r], col='green', lwd=2, lty=2)
	abline( v=lower[r], col='purple', lwd=2, lty=2)

  # Adding the boxplots
  plot_df = data.frame(yVar = post_means[,r])
  VP[[r]] = ggplot(plot_df, aes(x="", y = yVar)) +
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1) +
    ggtitle(labels[r]) +
    ylab(paste0("Parameter Value: ", true_par_inv[r])) +
    xlab(paste0("Coverage is: ", cov_df[r])) +
    geom_hline(yintercept=true_par_inv[r], linetype="dashed", color = "red") +
    theme(text = element_text(size = 7))
  #boxplot(post_means[,r], main=labels[r], ylab=NA, xlab = true_par[r])
  #abline(h=true_par[r], col="red")

}
grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]],
             VP[[14]], VP[[15]], VP[[16]], ncol=3, nrow =3)
dev.off()
