library(msm)
temp = msm::cav
cavData <- data.frame(  "ptnum" = temp$PTNUM,
                        "years" = temp$years,
                        "disc_time" = temp$years,
                        "sex" = temp$sex,
                        "state" = temp$state,
                        "obstrue" = rep(0, nrow(temp)))
