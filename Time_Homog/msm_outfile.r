requiredPackages = c('tidyverse','gridExtra')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

nFrames = 50

labels <- c('baseline State 1 (well)   --->   State 2 (mild)',
            'baseline State 1 (well)   --->   State 4 (dead)',
            'baseline State 2 (mild)   --->   State 3 (severe)',
            'baseline State 2 (mild)   --->   State 4 (dead)',
            'baseline State 3 (severe)   --->   State 4 (dead)',
            # 'iyears State 1 (well)   --->   State 2 (mild)',
            # 'iyears State 1 (well)   --->   State 4 (dead)',
            # 'iyears State 2 (mild)   --->   State 3 (severe)',
            # 'iyears State 2 (mild)   --->   State 4 (dead)',
            # 'iyears State 3 (severe)   --->   State 4 (dead)',
            'sex State 1 (well)   --->   State 2 (mild)',
            'sex State 1 (well)   --->   State 4 (dead)',
            'sex State 2 (mild)   --->   State 3 (severe)',
            'sex State 2 (mild)   --->   State 4 (dead)',
            'sex State 3 (severe)   --->   State 4 (dead)',
            'P( observed state 2 | true state 1 )',
            'P( observed state 1 | true state 2 )',
            'P( observed state 3 | true state 2 )',
            'P( observed state 2 | true state 3 )',
            'P( initial state 2 )','P( initial state 3 )')


day_data = matrix(data=-1, nrow = nFrames, ncol = 16)

for (i in 1:nFrames) {

  load(paste0("Model_out/msm/Output_msm", i, ".rda"))
  day_data[i,] = Output_msm$opt$par

}

trueValues <- c(c(matrix(c(-2.54, -0.56,
                           -2.94,  0.15,
                           -1.10, -0.03,
                           -3.92,  0.21,
                           -2.12,  1.17), ncol=2, byrow=T)),
                  c(  -4.59512, -1.15268, -2.751535, -2.090741),
                  c( -3.178054, -4.59512))
# 0*0.11,
# 0*-0.24,
# 0*-0.15,
# 0*0.23 ,
# 0*0.08 ,

# for(i in 1:5) {
#   day_data[,i] = day_data[,i] - day_data[,22]*trueValues[i+5]
#   month_data[,i] = month_data[,i] - month_data[,22]*trueValues[i+5]
#   year_data[,i] = year_data[,i] - year_data[,22]*trueValues[i+5]
#   year_2_data[,i] = year_2_data[,i] - year_2_data[,22]*trueValues[i+5]
# }



pdf("Plots/msm.pdf", onefile = T)
VP <- vector(mode="list", length = length(labels))
for(i in 1:length(trueValues)) {
    plot_df = data.frame(yVar = day_data[,i])
    VP[[i]] = ggplot(plot_df, aes(x="", y = yVar)) +
      geom_violin(trim=FALSE) +
      geom_boxplot(width=0.1) +
      ggtitle(labels[i]) +
      ylab('') +
      xlab(trueValues[i]) +
      geom_hline(yintercept=trueValues[i], linetype="dashed", color = "red") +
      theme(text = element_text(size = 7))

}
grid.arrange(VP[[1]], VP[[2]], VP[[3]], VP[[4]], VP[[5]],
             VP[[6]], VP[[7]], VP[[8]], VP[[9]], ncol=3, nrow =3)
grid.arrange(VP[[10]], VP[[11]], VP[[12]], VP[[13]],
             VP[[14]], VP[[15]], VP[[16]], ncol=3, nrow =3)
dev.off()
