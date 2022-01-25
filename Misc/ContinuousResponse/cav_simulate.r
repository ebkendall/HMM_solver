##############################
#### LOADING LIB AND DATA ####
##############################
package.install = function(pack) {
  local({r <- getOption("repos");r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)})

  # name of package to install / load
  pack = pack

  if (pack %in% rownames(installed.packages())) {
    library(pack, character.only=T)
  } else {
    if (pack %in% rownames(installed.packages(lib.loc='/blue/jantonelli/emmett.kendall/Packages/R_4_0'))) {
      library(pack, lib.loc='/blue/jantonelli/emmett.kendall/Packages/R_4_0', character.only=T)
    } else {
      install.packages(pack, lib='/blue/jantonelli/emmett.kendall/Packages/R_4_0')
      library(pack, lib.loc='/blue/jantonelli/emmett.kendall/Packages/R_4_0', character.only=T)
    }
  }
}

floor_new <- function(t,p) {
    new_time = NULL
    if(p == 1) {
      new_time = t
    } else if (p==2) {
      monthSeq = seq(0,1,1/6)
      yearNum = floor(t)
      timeInd = max(which(monthSeq <= (t - yearNum))) #selects which month to go to
      new_time = yearNum + monthSeq[timeInd]
    } else if (p==3) {
      new_time = floor(t)
    } else if (p==4) {
      if(floor(t) %% 2 == 0) { #divisible by 2
        new_time = floor(t)
      } else {
        temp = t - 1
        new_time = floor(temp)
      }
    } else {
      print("Invalid input for p")
    }
    return(new_time)
}

censor_times <- function(t, p) {
  min_t = 0
  max_t = floor_new(max(t), p)
  new_time = c()
  if(p == 1) {
    new_time = t
  } else if (p==2) {
    new_time = seq(min_t, max_t, by = 1/6)
  } else if (p==3) {
    new_time = seq(min_t, max_t, by = 1)
  } else {
    new_time = seq(min_t, max_t, by = 2)
  }
  return(new_time)
}

package.install("msm")
# library(msm)

num_iter = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# p = 1 --> update continuous
# p = 2 --> update every month
# p = 3 --> update every year
# p = 4 --> update every twoYear

set.seed(num_iter)
folder_name = c("Continuous", "Month", "Year", "YearTwo")


# Set the sample size.  Note that the cav data set has 622 subjects.
N <- 2000
# Choose the discretization of time.
dt <- 1/365


# The years and disc_time columns were both centered at round(mean(years),0) = 4 in the cav data set.
# These are the true parameter values for the uncentered data ( intercept - coef*mean ).

trueValues <- c(c(matrix(c(-2.29709805,  0.09266760, -0.56262135,
                         -1.17308794, -5.10636947, -0.96162312,
                         -1.71474254, -0.04338819,  0.83882558,
                         -2.08300714,  0.03824367, -2.75345311,
                         -2.42208380,  0.11315485,  1.76897841), ncol=3, byrow=T)),
                         c(  -5.60251814, -0.84455697, -2.56906519, -2.12629033),
                         c( -6.95125291, -7.07504453),
                         c(10, 20, 30, 40),
                         rep(1,4))

par_index = list( beta=1:15, misclass=16:19, pi_logit=20:21, mu = 22:25, sigma = 26:29)

betaMat <- matrix(trueValues[par_index$beta], ncol = 3, byrow = F)


errorMat_temp = matrix(c(1, exp(trueValues[par_index$misclass][1]), 0, 0,
                         exp(trueValues[par_index$misclass][2]), 1, exp(trueValues[par_index$misclass][3]), 0,
                         0, exp(trueValues[par_index$misclass][4]), 1, 0,
                         0,   0,   0, 1), ncol=4, byrow=TRUE)

errorMat = errorMat_temp / rowSums(errorMat_temp)


initProbs_temp = c( 1, exp(trueValues[par_index$pi_logit][1]), exp(trueValues[par_index$pi_logit][2]), 0)
initProbs = initProbs_temp / sum(initProbs_temp)


#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Collect information about the real data set.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------



ptnum <- unique(cav$PTNUM)
N_cav <- length(ptnum)
interObsTime <- NULL
propMale <- 0
propDeaths <- 0
NumObs <- rep(0,N_cav)
for(i in 1:N_cav){
  subject <- cav[cav$PTNUM==ptnum[i],,drop=FALSE]

  # The number of observations for each subject.
  NumObs[i] <- nrow(subject)

  # The times between observations.
  if(!(4 %in% subject$state)){  interObsTime <- c( interObsTime, round( diff(subject$years), 3))  }

  # Determine whether the subject is male.
  propMale <- propMale + as.integer(subject$sex[1]==1)

  # Determine whether the subject's death was observed.
  propDeaths <- propDeaths + as.integer(4 %in% subject$state)
}
propMale <- propMale / N_cav
propDeaths <- propDeaths / N_cav



#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Fill in the states using the transition rate matrix and error matrix estimated on the real data set.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

Q <- function(time,sex,betaMat){

  q1  = exp( c(1,time,sex) %*% betaMat[1,] )  # Transition from state 1 to state 2.
  q2  = exp( c(1,time,sex) %*% betaMat[2,] )  # Transition from state 2 to state 3.
  q3  = exp( c(1,time,sex) %*% betaMat[3,] )  # Transition from state 1 to death.
  q4  = exp( c(1,time,sex) %*% betaMat[4,] )  # Transition from state 2 to death.
  q5  = exp( c(1,time,sex) %*% betaMat[5,] )  # Transition from state 3 to death.

  qmat = matrix(c( 0,q1, 0,q2,
                   0, 0,q3,q4,
                   0, 0, 0,q5,
                   0, 0, 0, 0),nrow=4,byrow=TRUE)
  diag(qmat) = -rowSums(qmat)

  return(qmat)
}


rawData <- NULL
propDeaths_sim <- 0
NumObs_sim <- NULL
for(i in 1:N){

  print(i)

  # Sample the gender, as proportional to the cav data set.
  sex <- as.integer(runif(1,0,1) < propMale)

  # Sample for an initial state.
  trueState <- sample(1:4, size=1, prob=initProbs)

  # Sample an initial continuous mean
  m1 = trueValues[par_index$mu][trueState]
  s1 = trueValues[par_index$sigma][trueState]
  cont_resp = rnorm(1, m1, s1)

  # Sample the remaining states until death.
  years <- 0
  time1 <- 0
  s <- trueState
  r <- cont_resp
  while(s < 4){

    # Infinitesimal transition rates.
    qmat <- Q(time1,sex,betaMat)

    # Possible next states.
    moveToStates <- which(qmat[s,] > 0)

    # Sample the wait times before transition to each of the next possible states.
    waitTimes <- rexp( n=length(moveToStates), rate= qmat[s,moveToStates])

    # If any of the wait times are smaller than dt, then transition to the state with the minimum wait time.
    min_waitTime <- min(waitTimes)
    if(min_waitTime < dt){  s <- moveToStates[ which(waitTimes == min_waitTime) ]  }

    my_mean = trueValues[par_index$mu][s]
    my_sd = trueValues[par_index$sigma][s]
    r = rnorm(1, my_mean, my_sd)

    time1 <- time1 + dt

    years <- c( years, time1)
    trueState <- c( trueState, s)
    cont_resp <- c( cont_resp, r)

  }
  timeOfDeath <- tail(years,1)

  # Sample inter-observation times from the cav data set.  Maximum of 20 years in study.
  visitTimes <- NULL
  time2 <- 0

  while(time2 < min( 20, timeOfDeath)){

    visitTimes <- c( visitTimes, time2)
    time2 <- time2 + sample( interObsTime, size=1)
  }

  # If first visit time occurs after death, then subject is NOT entered into the study.
  if( !is.null(visitTimes) ){

    # If death occured before the study ended, then record the time of death.
    if( timeOfDeath < 20 ){  visitTimes <- c( visitTimes, timeOfDeath) }


    n_i <- length(visitTimes)
    state <- NULL
    resp <- NULL
    for(k in 1:n_i){
        state <- c( state, tail( trueState[ years <= visitTimes[k] ], 1))
        resp <-  c( resp, tail( cont_resp[ years <= visitTimes[k] ], 1))
    }

    ptnum <- rep(i,n_i)
    years <- visitTimes
    rawData <- rbind( rawData, data.frame(ptnum,years,sex,state, resp) )

    if(4 %in% state){  propDeaths_sim <- propDeaths_sim + 1  }
    NumObs_sim <- c( NumObs_sim, n_i)
  }

}

colnames(rawData) <- c('ptnum','years','sex','state','cont_resp')
N <- length(unique(rawData$ptnum))
propDeaths_sim <- propDeaths_sim / N


# Add noise to the states.
for(i in 1:nrow(rawData)){	rawData$state[i] <- sample(1:4, size=1, prob=errorMat[rawData$state[i],])  }

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Add censored rows.
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Key
# p = 1: continuous (no censor)
# p = 2: monthly censor
# p = 3: yearly censor
# p = 4: two year censor
for(p in 1:4) {

    disc_time <- sapply(rawData$years, floor_new, p = p)

    obstrue <- rep(0,nrow(rawData))

    hold <- cbind(rawData,obstrue,disc_time)
    hold <- hold[,c('ptnum','years','disc_time','sex','state','cont_resp','obstrue')]

    tempRow <- rep(0,ncol(hold))
    names(tempRow) <- c('ptnum','years','disc_time','sex','state','cont_resp','obstrue')

    num <- 1
    cavData <- NULL
    if(p != 1) {
        for(i in unique(rawData$ptnum)){

          current <- NULL
          subject <- hold[hold$ptnum==i,,drop=FALSE]

          #------------------------------------
          censoredAges <- censor_times(subject$years, p)

          for(t in censoredAges ){

            # Rounding t, subject$years, & subject$disc_time to make sure we have equality
            t_round = round(t, digits = 5)
            yrs_round = round(subject$years, digits = 5)
            disc_round = round(subject$disc_time, digits = 5)

            # If 't' corresponds to an observed age, then the next row will include the observed clinical visit data.
            if(t_round %in% yrs_round){
              current <- rbind( current, subject[disc_round==round(floor_new(t,p), digits=5),])
            } else{

              # Create a CENSORED row for each subject at each discritezed time.
              tempRow['ptnum'] <- i
              tempRow['years'] <- t
              tempRow['disc_time'] <- t
              tempRow['sex'] <- subject$sex[1]
              tempRow['state'] <- 99
              tempRow['cont_resp'] <- 99
              tempRow['obstrue'] <- 1

              current <- rbind( current, tempRow)

              # If 't' corresponds to an observed INTEGER years, then the subject was observed some time during this years.  According, the next row will include the observed clinical visit data.  Recall that integer years is simply the floor(years).
              if(t_round %in% disc_round){ current <- rbind( current, subject[disc_round==t_round,])}
            }

          }
          #------------------------------------

          cavData <- rbind( cavData, current)
          #print(num)
          num <- num+1
        }
    } else {
        cavData = hold
    }
    colnames(cavData) <- c('ptnum','years','disc_time','sex','state','cont_resp','obstrue')
    rownames(cavData) <- NULL

    save(cavData, file=paste("DataOut/", folder_name[p], "/cavData", num_iter, ".rda", sep=''))

    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    # Compute the frequencies of observed transitions.
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------

    # Transition frequencies for the true cav data set.
    if (p == 1) {
        nTrans <- rep(0,5)
        for(i in 1:N_cav){

        	subject <- cav[cav$PTNUM==unique(cav$PTNUM)[i],,drop=FALSE]

        	if(4 %in% subject$state){
        		if( max(setdiff(subject$state,4)) == 1){nTrans[2] = nTrans[2] +1}
        		if( max(setdiff(subject$state,4)) == 2){nTrans[4] = nTrans[4] +1}
        		if( max(setdiff(subject$state,4)) == 3){nTrans[5] = nTrans[5] +1}
        	}

        	if(   1 %in% subject$state   &   2 %in% subject$state   ){ nTrans[1] = nTrans[1] +1 }
        	if(   2 %in% subject$state   &   3 %in% subject$state   ){ nTrans[3] = nTrans[3] +1 }
        }

        cat('Cav data set sample size                               = ', N_cav,'\n')
        cat('Cav data set transition fequencies                     = ', nTrans / sum(nTrans),'\n')
        cat('Cav data set transition counts                         = ', nTrans,'\n')
        cat('Cav data set proportion of observed deaths             = ', propDeaths,'\n')
        cat('Cav data set quantiles of number of observations       = ','\n')
        print(quantile(NumObs))
        cat('\n')
    }
    # Transition frequencies for the simulated data set.
    obs_cavData <- cavData[cavData$obstrue == 0, ]
    nTrans_sim <- rep(0,5)
    for(i in unique(obs_cavData$ptnum)){

    	subject <- obs_cavData[obs_cavData$ptnum==i,,drop=FALSE]

    	if(4 %in% subject$state){
    		if( max(setdiff(subject$state,4)) == 1){nTrans_sim[2] = nTrans_sim[2] +1}
    		if( max(setdiff(subject$state,4)) == 2){nTrans_sim[4] = nTrans_sim[4] +1}
    		if( max(setdiff(subject$state,4)) == 3){nTrans_sim[5] = nTrans_sim[5] +1}
    	}

    	if(   1 %in% subject$state   &   2 %in% subject$state   ){ nTrans_sim[1] = nTrans_sim[1] +1 }
    	if(   2 %in% subject$state   &   3 %in% subject$state   ){ nTrans_sim[3] = nTrans_sim[3] +1 }
    }
    cat('Simulated data set sample size                         = ', N,'\n')
    cat('Simulated data set transition fequencies               = ', nTrans_sim / sum(nTrans_sim),'\n')
    cat('Simulated data set transition counts                   = ', nTrans_sim,'\n')
    cat('Simulated data set proportion of observed deaths       = ', propDeaths_sim,'\n')
    cat('Simulated data set quantiles of number of observations = ','\n')
    print(quantile(NumObs_sim))
    print(sum(cavData$obstrue == 0))

}
