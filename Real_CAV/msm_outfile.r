requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

nFrames = 100

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

parEst_msm_MONTH = upper_msm_MONTH = lower_msm_MONTH = 
    parEst_msm_YEAR = upper_msm_YEAR = lower_msm_YEAR = 
    parEst_msm_YEAR2 = upper_msm_YEAR2 = lower_msm_YEAR2 = NULL

for (i in 1:nFrames) {

  print(i)

  # ---------------------------------------------------------------------------
  load(paste0("Model_out/Month_msm/Output_msm", i, ".rda"))
  month_data[i,] = c(Output_msm$opt$par)
  # saving the 95% confidence interval
  temp <- NULL
  temp_upper_msm <- NULL
  temp_lower_msm <- NULL
  for(r in c(5,13,10,14,15)){  
    for(l in 1:3){  
      temp <- c( temp, Output_msm$Qmatrices[[l]][r])  
      temp_upper_msm <- c( temp_upper_msm, 
                  Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r])
      temp_lower_msm <- c( temp_lower_msm, 
                  Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r])
    }  
  }
  for(r in c(5,2,10,7)){
      temp <- c( temp, Output_msm$Ematrices[[1]][r])  
      temp_upper_msm <- c( temp_upper_msm, Output_msm$EmatricesU[[1]][r])
      temp_lower_msm <- c( temp_lower_msm, Output_msm$EmatricesL[[1]][r])
  }
  temp <- c( temp, Output_msm$opt$par[20:21])
  temp_upper_msm <- c( temp_upper_msm, Output_msm$ci[36:37,2])
  temp_lower_msm <- c( temp_lower_msm, Output_msm$ci[36:37,1])

  parEst_msm_MONTH <- rbind( parEst_msm_MONTH, temp)
  upper_msm_MONTH <- rbind( upper_msm_MONTH, temp_upper_msm)
  lower_msm_MONTH <- rbind( lower_msm_MONTH, temp_lower_msm)

  # ---------------------------------------------------------------------------
  load(paste0("Model_out/Year_msm/Output_msm", i, ".rda"))
  year_data[i,] = c(Output_msm$opt$par)

  # saving the 95% confidence interval
  temp <- NULL
  temp_upper_msm <- NULL
  temp_lower_msm <- NULL
  for(r in c(5,13,10,14,15)){  
    for(l in 1:3){  
      temp <- c( temp, Output_msm$Qmatrices[[l]][r])  
      temp_upper_msm <- c( temp_upper_msm, 
                  Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r])
      temp_lower_msm <- c( temp_lower_msm, 
                  Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r])
    }  
  }
  for(r in c(5,2,10,7)){
      temp <- c( temp, Output_msm$Ematrices[[1]][r])  
      temp_upper_msm <- c( temp_upper_msm, Output_msm$EmatricesU[[1]][r])
      temp_lower_msm <- c( temp_lower_msm, Output_msm$EmatricesL[[1]][r])
  }
  temp <- c( temp, Output_msm$opt$par[20:21])
  temp_upper_msm <- c( temp_upper_msm, Output_msm$ci[36:37,2])
  temp_lower_msm <- c( temp_lower_msm, Output_msm$ci[36:37,1])

  parEst_msm_YEAR <- rbind( parEst_msm_YEAR, temp)
  upper_msm_YEAR  <- rbind( upper_msm_YEAR, temp_upper_msm)
  lower_msm_YEAR  <- rbind( lower_msm_YEAR, temp_lower_msm)

  # ---------------------------------------------------------------------------
  load(paste0("Model_out/YearTwo_msm/Output_msm", i, ".rda"))
  year_2_data[i,] = c(Output_msm$opt$par)

  # saving the 95% confidence interval
  temp <- NULL
  temp_upper_msm <- NULL
  temp_lower_msm <- NULL
  for(r in c(5,13,10,14,15)){  
    for(l in 1:3){  
      temp <- c( temp, Output_msm$Qmatrices[[l]][r])  
      temp_upper_msm <- c( temp_upper_msm, 
                  Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r])
      temp_lower_msm <- c( temp_lower_msm, 
                  Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r])
    }  
  }
  for(r in c(5,2,10,7)){
      temp <- c( temp, Output_msm$Ematrices[[1]][r])  
      temp_upper_msm <- c( temp_upper_msm, Output_msm$EmatricesU[[1]][r])
      temp_lower_msm <- c( temp_lower_msm, Output_msm$EmatricesL[[1]][r])
  }
  temp <- c( temp, Output_msm$opt$par[20:21])
  temp_upper_msm <- c( temp_upper_msm, Output_msm$ci[36:37,2])
  temp_lower_msm <- c( temp_lower_msm, Output_msm$ci[36:37,1])

  parEst_msm_YEAR2 <- rbind( parEst_msm_YEAR2, temp)
  upper_msm_YEAR2 <- rbind( upper_msm_YEAR2, temp_upper_msm)
  lower_msm_YEAR2 <- rbind( lower_msm_YEAR2, temp_lower_msm)

}
msm_info = list(par = parEst_msm_MONTH, upper = upper_msm_MONTH, lower = lower_msm_MONTH)
save(msm_info, file = "Plots/msm_month.rda")
msm_info = list(par = parEst_msm_YEAR, upper = upper_msm_YEAR, lower = lower_msm_YEAR)
save(msm_info, file = "Plots/msm_year.rda")
msm_info = list(par = parEst_msm_YEAR2, upper = upper_msm_YEAR2, lower = lower_msm_YEAR2)
save(msm_info, file = "Plots/msm_yearTwo.rda")

trueValues <- c(c(matrix(c(-2.54,  0.11, -0.56,
                           -2.94, -0.24,  0.15,
                           -1.10, -0.15, -0.03,
                           -3.92,  0.23,  0.21,
                           -2.12,  0.08,  1.17), ncol=3, byrow=T)),
                  c(  -4.59512, -1.15268, -2.751535, -2.090741),
                  c( -3.178054, -4.59512))

VP <- vector(mode="list", length = length(labels))
for(i in 1:length(trueValues)) {
	# print(trueValues[i])
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
      # geom_hline(yintercept=trueValues[i], linetype="dashed", color = "red") +
      theme(text = element_text(size = 7))

}

pdf("Plots/msm.pdf", onefile = T)
grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]],
             VP[[14]], VP[[15]], VP[[16]], ncol=3, nrow =3)
dev.off()


load_names = c("Plots/msm_month.rda", "Plots/msm_year.rda", "Plots/msm_yearTwo.rda")
for(i in 1:3) {
  print(load_names[i])
  load(load_names[i])

  upper_means = round(colMeans(msm_info[[2]]), digits=2)
  lower_means = round(colMeans(msm_info[[3]]), digits=2)

  for(j in 0:4) {
    print(paste0("[ ", lower_means[3*j + 1], ", ",
                  upper_means[3*j+1], " ]"))
  }
}