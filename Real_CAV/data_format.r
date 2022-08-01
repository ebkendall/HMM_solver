library(msm)

# Function for creating the new data with censored rows
new_data_gen <- function(d, dat_fr) {
    dat_fr$disc_time = dat_fr$years - (dat_fr$years %% d)

    newDat = NULL
    num = 0
    # Adding censored rows
    for(i in unique(dat_fr$ptnum)){
        
        current <- NULL
        subject <- dat_fr[dat_fr$ptnum==i,,drop=FALSE]
        
        #------------------------------------
        censoredAges <- seq(0, max(subject$years), d)  
        yrs_round = round(subject$years, digits = 5)
        disc_round = round(subject$disc_time, digits = 5)
        for(t in censoredAges){
            tempRow = NULL

            t_round = round(t, digits = 5)

            # If 't' corresponds to an observed age, then the next row will 
            # include the observed clinical visit data.
            if(t_round %in% yrs_round){	
                current <- rbind( current, subject[subject$disc_time==t_round,]) 
            } else{
            
                # Create a CENSORED row for each subject at each INTEGER year of years.
                tempRow['ptnum'] <- i
                tempRow['years'] <- t
                tempRow['sex'] <- subject$sex[1] 
                tempRow['state'] <- 99
                tempRow['obstrue'] <- 1  
                tempRow['disc_time'] <- t
                
                current <- rbind( current, tempRow)
                
                # If 't' corresponds to an observed INTEGER years, then the 
                # subject was observed some time during this years.  According, 
                # the next row will include the observed clinical visit data.  
                # Recall that integer years is simply the floor(years).
                if(t_round %in% disc_round){ current <- rbind( current, subject[disc_round==t_round,]) }
            }

        }
        #------------------------------------
        
        newDat <- rbind( newDat, current)
        rownames(newDat) = NULL
        num <- num+1
    }

    dat_fr = newDat
    return(dat_fr)
}

real_data = msm::cav
real_data = real_data[,c('PTNUM', 'years', 'sex', 'state')]
real_data$obstrue = 0
colnames(real_data) = c('ptnum', 'years', 'sex', 'state', 'obstrue')

# There will be four different scenarios on which I need to format the data
# 1) no censored rows
# 2) censoring every year
# 3) censoring every two years
# 4) censoring every two months

# ----------------------------------- Case 1 -----------------------------------
# colnames(cavData)= c("ptnum", "years", "sex", "state", "obstrue", "disc_time")
cavData = real_data
cavData$disc_time = cavData$years
save(cavData, file = 'Data/cavData_cont.rda')

# ----------------------------------- Case 2 -----------------------------------
# ---------------------------------- (yearly) ----------------------------------
d = 1
cavData = new_data_gen(d, real_data)
save(cavData, file = 'Data/cavData_year.rda')

# ----------------------------------- Case 3 -----------------------------------
# --------------------------------- (2 years) ----------------------------------
d = 2
cavData = new_data_gen(d, real_data)
save(cavData, file = 'Data/cavData_twoYear.rda')

# ----------------------------------- Case 4 -----------------------------------
# --------------------------------- (2 months) ---------------------------------
d = 1/6
cavData = new_data_gen(d, real_data)
save(cavData, file = 'Data/cavData_month.rda')
