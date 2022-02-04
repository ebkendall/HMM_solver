requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

true_par = c(c(matrix(c(-2.26568339, 3 *  0.08766060, -0.49991746,
                  -1.22022878, 3 * -4.44888558, -0.82779213,
                  -1.56180104, 3 * -0.08262607,  0.73838829,
                  -2.20978996, 3 *  0.05404948, -1.83682627,
                  -2.41222255, 3 *  0.10833734,  1.63135439), ncol=3, byrow=T)),
          c(  -5.73343061, -0.78623894, -2.52747176, -2.12144526),
          c( -6.52842355, -6.15970066))

par_index = list(beta=1:15, misclass=16:19, pi_logit=20:21)

# Doing the inverse logit for true_par
true_par[par_index$pi_logit] =
    exp(true_par[par_index$pi_logit])/(1 + exp(true_par[par_index$pi_logit]))
true_par[par_index$misclass[1]] =
    exp(true_par[par_index$misclass[1]])/(1 + exp(true_par[par_index$misclass[1]]))
true_par[par_index$misclass[2:3]] =
    exp(true_par[par_index$misclass[2:3]])/sum(c(1, exp(true_par[par_index$misclass[2:3]])))
true_par[par_index$misclass[4]] =
    exp(true_par[par_index$misclass[4]])/(1 + exp(true_par[par_index$misclass[4]]))

print(true_par)
print(length(true_par))

labels <- c('baseline: S1 (well)   --->   S2 (mild)',
            'baseline: S1 (well)   --->   S4 (dead)',
            'baseline: S2 (mild)   --->   S3 (severe)',
            'baseline: S2 (mild)   --->   S4 (dead)',
            'baseline: S3 (severe)   --->   S4 (dead)',
            'time S1 (well)   --->   S2 (mild)',
            'time S1 (well)   --->   S4 (dead)',
            'time S2 (mild)   --->   S3 (severe)',
            'time S2 (mild)   --->   S4 (dead)',
            'time S3 (severe)   --->   S4 (dead)',
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

# ---------------------------------------------------------
# -------- load and format post_means and coverage --------
# ---------------------------------------------------------
load("Plots/cov_df_deSolve.rda")
cov_final = cov_df
load("Plots/cov_df_expm.rda")
cov_final = cbind(cov_final, cov_df)
load("Plots/cov_df_msm.rda")
cov_final = cbind(cov_final, cov_df)
colnames(cov_final) = c("deSolve", "expm_Month", "expm_Year", "expm_BiYear",
                        "msm_Month", "msm_year","msm_BiYear")
cov_final = round(cov_final, digits = 3)

post_means_final = vector(mode = 'list', length = 7)
load("Plots/post_means_deSolve.rda")
post_means_final[[1]] = post_means[[1]]
load("Plots/post_means_expm.rda")
post_means_final[[2]] = post_means[[1]]; post_means_final[[3]] = post_means[[2]]
post_means_final[[4]] = post_means[[3]]
load("Plots/post_means_msm.rda")
post_means_final[[5]] = post_means[[1]]; post_means_final[[6]] = post_means[[2]]
post_means_final[[7]] = post_means[[3]]

VP <- vector(mode="list", length = length(labels))
pdf("Plots/violinPlots_final_bl.pdf")
for(r in 1:length(labels)) {

    yVar = disc_type = x_label = NULL

    yVar = c(post_means_final[[1]][,r], post_means_final[[2]][,r], post_means_final[[3]][,r],
             post_means_final[[4]][,r], post_means_final[[5]][,r], post_means_final[[6]][,r],
             post_means_final[[7]][,r])

    disc_type = c(rep(paste0("deSolve\n", cov_final[r,1]), nrow(post_means_final[[1]])),
                  rep(paste0("expm_Month\n",cov_final[r,2]), nrow(post_means_final[[2]])),
                  rep(paste0("expm_Year\n", cov_final[r,3]), nrow(post_means_final[[3]])),
                  rep(paste0("expm_Year_2\n",cov_final[r,4]), nrow(post_means_final[[4]])),
                  rep(paste0("msm_Month\n",cov_final[r,5]), nrow(post_means_final[[5]])),
                  rep(paste0("msm_Year\n",cov_final[r,6]), nrow(post_means_final[[6]])),
                  rep(paste0("msm_Year_2\n",cov_final[r,7]), nrow(post_means_final[[7]])))

    Method = c(rep("deSolve", nrow(post_means_final[[1]])),
                  rep("expm", nrow(post_means_final[[2]]) + nrow(post_means_final[[3]]) + nrow(post_means_final[[4]])),
                  rep("msm", nrow(post_means_final[[5]]) + nrow(post_means_final[[6]]) + nrow(post_means_final[[7]])))

    plot_df = data.frame(yVar = yVar, disc_type = disc_type, Method = Method)
    VP[[r]] = ggplot(plot_df, aes(x=disc_type, y = yVar, fill = Method)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[r]) +
      ylab("") + #paste0("Parameter Value: ", round(true_par[r], 3))
      xlab("") +
      geom_hline(yintercept=true_par[r], linetype="dashed", color = "red", size=0.75) +
      theme(text = element_text(size = 10))
}

grid.arrange(VP[[1]] ncol = 1, nrow = 1) #, VP[[2]], VP[[3]], VP[[4]], VP[[5]],
# grid.arrange(VP[[6]], VP[[7]], VP[[8]], VP[[9]], VP[[10]], ncol = 1, nrow = 5)
# grid.arrange(VP[[11]], VP[[12]], VP[[13]], VP[[14]], VP[[15]], ncol = 1, nrow = 5)
# grid.arrange(VP[[16]], VP[[17]], VP[[18]], VP[[19]], ncol = 1, nrow = 4)
# grid.arrange(VP[[20]], VP[[21]], ncol = 1, nrow = 2)

dev.off()
