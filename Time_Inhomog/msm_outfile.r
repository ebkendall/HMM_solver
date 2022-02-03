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

par_index = list(beta=1:15, misclass=16:19, pi_logit=20:21)

trueValues <- c(c(matrix(c(-2.26568339, 3 *  0.08766060, -0.49991746,
                           -1.22022878, 3 * -4.44888558, -0.82779213,
                           -1.56180104, 3 * -0.08262607,  0.73838829,
                           -2.20978996, 3 *  0.05404948, -1.83682627,
                           -2.41222255, 3 *  0.10833734,  1.63135439), ncol=3, byrow=T)),
                  c(  -5.73343061, -0.78623894, -2.52747176, -2.12144526),
                  c( -6.52842355, -6.15970066))

# trueValues <-  c(c(matrix(c(-2.54,  0.11, -0.56,
#                             -2.94, -0.24,  0.15,
#                             -1.10, -0.15, -0.03,
#                             -3.92,  0.23,  0.21,
#                             -2.12,  0.08,  1.17), ncol=3, byrow=T)),
#                  c(  -4.59512, -1.15268, -2.751535, -2.090741),
#                  c( -3.178054, -4.59512))

# trueErrors <- c( .01, .24, .06, .11)
# trueInitProbs <- c( .04, .01)

month_data = matrix(data=-1, nrow = nFrames, ncol = 21)
year_data = matrix(data=-1, nrow = nFrames, ncol = 21)
year_2_data = matrix(data=-1, nrow = nFrames, ncol = 21)

cred_set = vector(mode = "list", length = length(trueValues))

# For each parameter, there are 3 discretizations to consider
for(i in 1:length(cred_set)) {
    cred_set[[i]] = vector(mode = "list", length = 3)
    cred_set[[i]][[1]] = cred_set[[i]][[2]] = cred_set[[i]][[3]] =
            data.frame('lower' = c(-1), 'upper' = c(-1))
}

# -----------------------------------------------------------------------------
# Load Data and Find Credible Sets --------------------------------------------
# -----------------------------------------------------------------------------

row_ind1 = row_ind2 = row_ind3 = 1

for (i in 1:nFrames) {

  load(paste0("Model_out/msm/Month/Output_msm", i, ".rda"))
  if (length(Output_msm$QmatricesSE) != 0) { # means a non poisitive definite matrix
      month_data[row_ind1,] = c(Output_msm$opt$par)

      # Calculating credible sets (Month)
      ind = 1
      for(l in 1:3) {
          for(r in c(5,13,10,14,15)){
              cred_set[[ind]][[1]][row_ind1,1] = Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r]
              cred_set[[ind]][[1]][row_ind1,2] = Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r]
              ind = ind + 1
          }
      }
      # Probabilities
      ind = 1
      lowerBounds = Output_msm$EmatricesL[[1]][c(5,2,10,7)]
      upperBounds = Output_msm$EmatricesU[[1]][c(5,2,10,7)]
      for(l in par_index$misclass) {
          cred_set[[l]][[1]][row_ind1,1] = lowerBounds[ind]
          cred_set[[l]][[1]][row_ind1,2] = upperBounds[ind]
          ind = ind + 1
      }
      cred_set[[20]][[1]][row_ind1,1] = Output_msm$ci[36,1]; cred_set[[20]][[1]][row_ind1,2] = Output_msm$ci[36,2]
      cred_set[[21]][[1]][row_ind1,1] = Output_msm$ci[37,1]; cred_set[[21]][[1]][row_ind1,2] = Output_msm$ci[37,2]

      row_ind1 = row_ind1 + 1
  } else {
      print(paste0("Month Issue: ", i))
  }


 # -----------------------------------------------------------------------------

  load(paste0("Model_out/msm/Year/Output_msm", i, ".rda"))
  if (length(Output_msm$QmatricesSE) != 0) {
      year_data[row_ind2,] = c(Output_msm$opt$par)

      # Calculating credible sets (Year)
      ind = 1
      for(l in 1:3) {
          for(r in c(5,13,10,14,15)){
              cred_set[[ind]][[2]][row_ind2,1] = Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r]
              cred_set[[ind]][[2]][row_ind2,2] = Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r]
              ind = ind + 1
          }
      }
      # Probabilities
      ind = 1
      lowerBounds = Output_msm$EmatricesL[[1]][c(5,2,10,7)]
      upperBounds = Output_msm$EmatricesU[[1]][c(5,2,10,7)]
      for(l in par_index$misclass) {
          cred_set[[l]][[2]][row_ind2,1] = lowerBounds[ind]
          cred_set[[l]][[2]][row_ind2,2] = upperBounds[ind]
          ind = ind + 1
      }
      cred_set[[20]][[2]][row_ind2,1] = Output_msm$ci[36,1]; cred_set[[20]][[2]][row_ind2,2] = Output_msm$ci[36,2]
      cred_set[[21]][[2]][row_ind2,1] = Output_msm$ci[37,1]; cred_set[[21]][[2]][row_ind2,2] = Output_msm$ci[37,2]

      row_ind2 = row_ind2 + 1
  } else {
      print(paste0("Year Issue: ", i))
  }
 # -----------------------------------------------------------------------------

  load(paste0("Model_out/msm/YearTwo/Output_msm", i, ".rda"))
  if (length(Output_msm$QmatricesSE) != 0) {
      year_2_data[row_ind3,] = c(Output_msm$opt$par)

      # Calculating credible sets (YearTwo)
      ind = 1
      for(l in 1:3) {
          for(r in c(5,13,10,14,15)){
              cred_set[[ind]][[3]][row_ind3,1] = Output_msm$Qmatrices[[l]][r] - 1.96*Output_msm$QmatricesSE[[l]][r]
              cred_set[[ind]][[3]][row_ind3,2] = Output_msm$Qmatrices[[l]][r] + 1.96*Output_msm$QmatricesSE[[l]][r]
              ind = ind + 1
          }
      }

      # Probabilities
      ind = 1
      lowerBounds = Output_msm$EmatricesL[[1]][c(5,2,10,7)]
      upperBounds = Output_msm$EmatricesU[[1]][c(5,2,10,7)]
      for(l in par_index$misclass) {
          cred_set[[l]][[3]][row_ind3,1] = lowerBounds[ind]
          cred_set[[l]][[3]][row_ind3,2] = upperBounds[ind]
          ind = ind + 1
      }
      cred_set[[20]][[3]][row_ind3,1] = Output_msm$ci[36,1]; cred_set[[20]][[3]][row_ind3,2] = Output_msm$ci[36,2]
      cred_set[[21]][[3]][row_ind3,1] = Output_msm$ci[37,1]; cred_set[[21]][[3]][row_ind3,2] = Output_msm$ci[37,2]

      row_ind3 = row_ind3 + 1
  } else {
      print(paste0("YearTwo Issue: ", i))
  }
}

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

print(trueValues)

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

# -----------------------------------------------------------------------------
# Calculating Coverage --------------------------------------------------------
# -----------------------------------------------------------------------------

cov_df = matrix(ncol = 3, nrow = length(trueValues))
colnames(cov_df) = c("Month", "Year", "YearTwo")

for(i in 1:length(trueValues)) {
    val = trueValues[i]
    for(j in 1:ncol(cov_df)) {
        top = length(which(cred_set[[i]][[j]]$lower <= val & val <= cred_set[[i]][[j]]$upper))
        bot = nrow(cred_set[[i]][[j]])
        covrg = top/bot
        cov_df[i,j] = covrg
    }
    print(paste0("Coverage for ", round(val, 3), " is: "))
    cat('\t', '\t', '\t', "Month:   ", cov_df[i,1], '\n')
    cat('\t', '\t', '\t', "Year:    ", cov_df[i,2], '\n')
    cat('\t', '\t', '\t', "YearTwo: ", cov_df[i,3], '\n')
}

# ------------------------------------------------------------------------------

VP <- vector(mode="list", length = length(labels))
for(i in 1:length(trueValues)) {
	print(round(trueValues[i], 3))
	m = data.frame(y = month_data[,i], type = rep("Month", nrow(month_data)))
	y = data.frame(y = year_data[,i], type = rep("Year", nrow(year_data)))
	y2 = data.frame(y = year_2_data[,i], type = rep("Year_2", nrow(year_2_data)))

    plot_df = rbind(m,y,y2)
    xlabel = paste0("Month: ", cov_df[i,1], ", Year: ", cov_df[i,2], ", YearTwo: ", cov_df[i,3])

    VP[[i]] = ggplot(plot_df, aes(x = type, y = y)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[i]) +
      ylab('') +
      xlab(xlabel) +
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
