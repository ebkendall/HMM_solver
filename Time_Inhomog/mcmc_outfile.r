# This script file produces trace plots and histograms of the mcmc output files
# To run from command line...
requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

args = commandArgs(TRUE)
# num_seeds = 7 #args[1] # Only used in this script to locate and load files
folder = as.numeric(args[1])

model_name = c('deSolve', 'expm')
sub_folder = c('Year/', 'YearTwo/') #'Month/',
for_length = c(1, 2)

dir = paste0('Model_out/', model_name[folder], '/')

# Size of posterior sample from mcmc chains
n_post = 5000
# Step number at 3ich the adaptive tuning scheme was frozen
burnin = 5000
# Total number of steps the mcmc algorithm is computed for
steps =  10000
# Matrix row indices for the posterior sample to use for GFF computation
index_post = (steps - burnin - n_post + 1):(steps - burnin)

par_index = list(beta=1:15, misclass=16:19, pi_logit=20:21)

true_par = c(c(matrix(c(-2.54,  0.11, -0.56,
                        -2.94, -0.24,  0.15,
                        -1.10, -0.15, -0.03,
                        -3.92,  0.23,  0.21,
                        -2.12,  0.08,  1.17), ncol=3, byrow=T)),
                  c(  -4.59512, -1.15268, -2.751535, -2.090741),
                  c( -3.178054, -4.59512))

# Doing the inverse logit for true_par
true_par[par_index$pi_logit] =
    exp(true_par[par_index$pi_logit])/(1 + exp(true_par[par_index$pi_logit]))
true_par[par_index$misclass[1]] =
    exp(true_par[par_index$misclass[1]])/(1 + exp(true_par[par_index$misclass[1]]))
true_par[par_index$misclass[2:3]] =
    exp(true_par[par_index$misclass[2:3]])/sum(c(1, exp(true_par[par_index$misclass[2:3]])))
true_par[par_index$misclass[4]] =
    exp(true_par[par_index$misclass[4]])/(1 + exp(true_par[par_index$misclass[4]]))



labels <- c('b.l. S1 (well)   --->   S2 (mild)',
            'b.l. S1 (well)   --->   S4 (dead)',
            'b.l. S2 (mild)   --->   S3 (severe)',
            'b.l. S2 (mild)   --->   S4 (dead)',
            'b.l. S3 (severe)   --->   S4 (dead)',
            'iyears S1 (well)   --->   S2 (mild)',
            'iyears S1 (well)   --->   S4 (dead)',
            'iyears S2 (mild)   --->   S3 (severe)',
            'iyears S2 (mild)   --->   S4 (dead)',
            'iyears S3 (severe)   --->   S4 (dead)',
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
for(i in 1:length(cred_set)) {
    if(folder == 1) {
        cred_set[[i]] = data.frame('lower' = c(-1), 'upper' = c(-1))
    } else {
        # Now there are 3 credible sets per parameter for different discretizations
        cred_set[[i]] = vector(mode = 'list', length = length(sub_folder))
        for(j in 1:length(cred_set)) {cred_set[[i]][[j]] = data.frame('lower' = c(-1), 'upper' = c(-1))}
    }

}

for (i in 1:50) {
    file_name = paste0(dir,'mcmc_out_',toString(i),'.rda')
    # If its deSolve, for_length[folder] = 1
    # If its expm, for_length[folder] = 3 b/c of the 3 discretizations
    loopLength = for_length[folder]
    for(w in 1:loopLength) {
        if (folder == 2) {file_name = paste0(dir,sub_folder[w],'mcmc_out_',toString(i),'.rda')}

	    load(file_name)

        # Inverse logit to convert back to probabilities
        mcmc_out$chain[,par_index$pi_logit] =
            exp(mcmc_out$chain[,par_index$pi_logit])/(1 + exp(mcmc_out$chain[,par_index$pi_logit]))
        mcmc_out$chain[,par_index$misclass[1]] =
            exp(mcmc_out$chain[,par_index$misclass[1]])/(1 + exp(mcmc_out$chain[,par_index$misclass[1]]))
        mcmc_out$chain[,par_index$misclass[2:3]] = exp(mcmc_out$chain[,par_index$misclass[2:3]])/(1 +
                                                   exp(mcmc_out$chain[,par_index$misclass[2]]) +
                                                   exp(mcmc_out$chain[,par_index$misclass[3]]))
        mcmc_out$chain[,par_index$misclass[4]] =
            exp(mcmc_out$chain[,par_index$misclass[4]])/(1 + exp(mcmc_out$chain[,par_index$misclass[4]]))

        for(j in 1:length(true_par)) {
            if(folder == 2) {
                cred_set[[j]][[w]][i,1] =  round(quantile( mcmc_out$chain[index_post,j],
                                            prob=.025), 4)
                cred_set[[j]][[w]][i,2] =  round(quantile( mcmc_out$chain[index_post,j],
                                            prob=.975), 4)
            } else {
                cred_set[[j]][i,1] =  round(quantile( mcmc_out$chain[index_post,j],
                                            prob=.025), 4)
                cred_set[[j]][i,2] =  round(quantile( mcmc_out$chain[index_post,j],
                                            prob=.975), 4)
            }
        }
    }
}

save(cred_set, file = paste0('Plots/cred_set_', model_name[folder], '.rda'))

# -----------------------------------------------------------------------------
# Calculating Coverage --------------------------------------------------------
# -----------------------------------------------------------------------------
if (folder == 1) {
    cov_df <- c()
    for(i in 1:length(true_par)) {
        val = true_par[i]
        top = length(which(cred_set[[i]]$lower <= val & val <= cred_set[[i]]$upper))
        bot = nrow(cred_set[[i]])
        covrg = top/bot
        cov_df[i] = covrg
        print(paste0("Coverage for parameter ", val, " is: ", covrg))
    }
} else {
    cov_df = matrix(ncol = length(sub_folder)); colnames(cov_df) = sub_folder
    for(i in 1:length(true_par)) {
        val = true_par[i]
        for(j in 1:ncol(cov_df)) {
            top = length(which(cred_set[[i]][[j]]$lower <= val & val <= cred_set[[i]][[j]]$upper))
            bot = nrow(cred_set[[i]][[j]])
            covrg = top/bot
            cov_df[i,j] = covrg
            print(paste0("Coverage for parameter ", val, " is: ", covrg))
        }
    }
}

# -----------------------------------------------------------------------------
# Create mcmc trace plots and histograms
# -----------------------------------------------------------------------------

index_seeds = 1:50
post_means = vector(mode = "list", length = length(sub_folder))
chain_list <- vector(mode = "list", length = length(sub_folder))

# If folder == 1, then we will only populate the first entry of the list
for(i in 1:length(chain_list)) {chain_list[[i]] = vector(mode = "list", length = length(index_seeds))}
for(i in 1:length(post_means)) {post_means[[i]] = matrix(nrow = length(index_seeds), ncol = length(labels))}

ind = 0
for(seed in index_seeds){
    file_name = paste0(dir,'mcmc_out_',toString(seed),'.rda')

    for(w in 1:for_length[folder]) {
        if (folder == 2) {file_name = paste0(dir,sub_folder[w],'mcmc_out_',toString(i),'.rda')}

        load(file_name)
        ind = ind + 1

        print(mcmc_out$accept)

        # Inverse logit to convert back to probabilities
        mcmc_out$chain[,par_index$pi_logit] =
            exp(mcmc_out$chain[,par_index$pi_logit])/(1 + exp(mcmc_out$chain[,par_index$pi_logit]))
        mcmc_out$chain[,par_index$misclass[1]] =
            exp(mcmc_out$chain[,par_index$misclass[1]])/(1 + exp(mcmc_out$chain[,par_index$misclass[1]]))
        mcmc_out$chain[,par_index$misclass[2:3]] = exp(mcmc_out$chain[,par_index$misclass[2:3]])/(1 +
                                                   exp(mcmc_out$chain[,par_index$misclass[2]]) +
                                                   exp(mcmc_out$chain[,par_index$misclass[3]]))
        mcmc_out$chain[,par_index$misclass[4]] =
            exp(mcmc_out$chain[,par_index$misclass[4]])/(1 + exp(mcmc_out$chain[,par_index$misclass[4]]))


      	chain_list[[w]][[seed]] = mcmc_out$chain[index_post,]
    	post_means[[w]][seed,] <- colMeans(mcmc_out$chain[index_post,])
    }
}

# Plot and save the mcmc trace plots and histograms.
pdf(paste0('Plots/mcmc_', model_name[folder], '.pdf'))
par(mfrow=c(4, 2))

for (q in 1:for_length[folder]) {
    stacked_chains = do.call( rbind, chain_list[[q]])
    par_mean = par_median = upper = lower = rep( NA, ncol(stacked_chains))
    VP <- vector(mode="list", length = length(labels))

    for(r in 1:length(labels)){

    	plot( NULL, xlab=NA, ylab=NA, main=labels[r], xlim=c(1,length(index_post)),
    		    ylim=range(stacked_chains[,r]) )

    	for(seed in 1:length(chain_list[[q]])) lines( chain_list[[q]][[seed]][,r], type='l', col=seed)

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

    }

}

for(r in 1:length(labels)) {
    # Adding the boxplots
    yVar = disc_type = x_label = NULL
    if(folder == 2) {
        yVar = c(post_means[[1]][,r], post_means[[2]][,r]) #, post_means[[3]][,r])
        disc_type = c(rep("Month", nrow(post_means[[1]])),
                      rep("Year", nrow(post_means[[2]])),
                      rep("YearTwo", nrow(post_means[[3]])))
        x_label = paste0("Month: ", cov_df[i,1], ", Year: ", cov_df[i,2]) #, ", YearTwo: ", cov_df[i,3])
    } else {
        yVar = post_means[[1]][,r]
        disc_type = rep("Continuous", nrow(post_means[[1]]))
        x_label = paste0("Coverage is: ", cov_df[r])
    }

    plot_df = data.frame(yVar = yVar, disc_type = disc_type)
    VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[r]) +
      ylab(paste0("Parameter Value: ", true_par[r])) +
      xlab(x_label) +
      geom_hline(yintercept=true_par[r], linetype="dashed", color = "red") +
      theme(text = element_text(size = 7))
}

grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]], VP[[14]],
             VP[[15]], VP[[16]], VP[[17]], VP[[18]], ncol=3, nrow =3)
grid.arrange(VP[[19]], VP[[20]], VP[[21]], ncol=3, nrow =3)

dev.off()
