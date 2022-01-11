requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

nFrames = 50

labels <- c('baseline S1 (well)   --->   S2 (mild)',
            'baseline S1 (well)   --->   S4 (dead)',
            'baseline S2 (mild)   --->   S3 (severe)',
            'baseline S2 (mild)   --->   S4 (dead)',
            'baseline S3 (severe)   --->   S4 (dead)',
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
            'P( observed state 2 | true state 1 )',
            'P( observed state 1 | true state 2 )',
            'P( observed state 3 | true state 2 )',
            'P( observed state 2 | true state 3 )',
            'P( initial state 2 )','P( initial state 3 )')


month_data = matrix(data=-1, nrow = nFrames, ncol = 21)
year_data = matrix(data=-1, nrow = nFrames, ncol = 21)
year_2_data = matrix(data=-1, nrow = nFrames, ncol = 21)

for (i in 1:nFrames) {

  print(i)

  load(paste0("Model_out/msm/Month/Output_msm", i, ".rda"))
  month_data[i,] = c(Output_msm$opt$par)

  load(paste0("Model_out/msm/Year/Output_msm", i, ".rda"))
  year_data[i,] = c(Output_msm$opt$par)

  load(paste0("Model_out/msm/YearTwo/Output_msm", i, ".rda"))
  year_2_data[i,] = c(Output_msm$opt$par)

}

par_index = list(beta=1:15, misclass=16:19, pi_logit=20:21)

trueValues <- c(c(matrix(c(-2.54,  0.11, -0.56,
                           -2.94, -0.24,  0.15,
                           -1.10, -0.15, -0.03,
                           -3.92,  0.23,  0.21,
                           -2.12,  0.08,  1.17), ncol=3, byrow=T)),
                  c(  -4.59512, -1.15268, -2.751535, -2.090741),
                  c( -3.178054, -4.59512))

# ------------------------------------------------------------------------------
# Doing the inverse logit ------------------------------------------------------
# ------------------------------------------------------------------------------
# trueValues transformation
trueValues[par_index$pi_logit] =
  exp(trueValues[par_index$pi_logit])/(1 + exp(trueValues[par_index$pi_logit]))
trueValues[par_index$misclass[1]] =
  exp(trueValues[par_index$misclass[1]])/(1 + exp(trueValues[par_index$misclass[1]]))
trueValues[par_index$misclass[2:3]] =
  exp(trueValues[par_index$misclass[2:3]])/sum(c(1, exp(trueValues[par_index$misclass[2:3]])))
trueValues[par_index$misclass[4]] =
  exp(trueValues[par_index$misclass[4]])/(1 + exp(trueValues[par_index$misclass[4]]))

# Month transformation
month_data[,par_index$pi_logit] =
    exp(month_data[,par_index$pi_logit])/(1 + exp(month_data[,par_index$pi_logit]))
month_data[,par_index$misclass[1]] =
    exp(month_data[,par_index$misclass[1]])/(1 + exp(month_data[,par_index$misclass[1]]))
month_data[,par_index$misclass[2:3]] = exp(month_data[,par_index$misclass[2:3]])/(1 +
                                           exp(month_data[,par_index$misclass[2]]) +
                                           exp(month_data[,par_index$misclass[3]]))
month_data[,par_index$misclass[4]] =
    exp(month_data[,par_index$misclass[4]])/(1 + exp(month_data[,par_index$misclass[4]]))

# Year transformation
year_data[,par_index$pi_logit] =
    exp(year_data[,par_index$pi_logit])/(1 + exp(year_data[,par_index$pi_logit]))
year_data[,par_index$misclass[1]] =
    exp(year_data[,par_index$misclass[1]])/(1 + exp(year_data[,par_index$misclass[1]]))
year_data[,par_index$misclass[2:3]] = exp(year_data[,par_index$misclass[2:3]])/(1 +
                                           exp(year_data[,par_index$misclass[2]]) +
                                           exp(year_data[,par_index$misclass[3]]))
year_data[,par_index$misclass[4]] =
    exp(year_data[,par_index$misclass[4]])/(1 + exp(year_data[,par_index$misclass[4]]))

# YearTwo transformation
year_2_data[,par_index$pi_logit] =
    exp(year_2_data[,par_index$pi_logit])/(1 + exp(year_2_data[,par_index$pi_logit]))
year_2_data[,par_index$misclass[1]] =
    exp(year_2_data[,par_index$misclass[1]])/(1 + exp(year_2_data[,par_index$misclass[1]]))
year_2_data[,par_index$misclass[2:3]] = exp(year_2_data[,par_index$misclass[2:3]])/(1 +
                                           exp(year_2_data[,par_index$misclass[2]]) +
                                           exp(year_2_data[,par_index$misclass[3]]))
year_2_data[,par_index$misclass[4]] =
    exp(year_2_data[,par_index$misclass[4]])/(1 + exp(year_2_data[,par_index$misclass[4]]))

# ------------------------------------------------------------------------------


VP <- vector(mode="list", length = length(labels))
for(i in 1:length(trueValues)) {
	print(trueValues[i])
	m = data.frame(y = month_data[,i], type = rep("Month", nrow(month_data)))
	y = data.frame(y = year_data[,i], type = rep("Year", nrow(year_data)))
	y2 = data.frame(y = year_2_data[,i], type = rep("Year_2", nrow(year_2_data)))

    plot_df = rbind(m,y,y2)

    VP[[i]] = ggplot(plot_df, aes(x = type, y = y)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[i]) +
      ylab('') +
      xlab(trueValues[i]) +
      geom_hline(yintercept=trueValues[i], linetype="dashed", color = "red") +
      theme(text = element_text(size = 7))

}

pdf("Plots/msm.pdf", onefile = T)
grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]], VP[[14]],
             VP[[15]], VP[[16]], VP[[17]], VP[[18]], ncol=3, nrow =3)
grid.arrange(VP[[19]], VP[[20]], VP[[21]], ncol=3, nrow =3)
dev.off()
