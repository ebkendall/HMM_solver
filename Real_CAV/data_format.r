library(msm)

real_data = msm::cav

# There will be four different scenarios on which I need to format the data
# 1) no censored rows
# 2) censoring every year
# 3) censoring every two years
# 4) censoring every two months

# ----------------------------------- Case 1 -----------------------------------

# ----------------------------------- Case 2 -----------------------------------
# ---------------------------------- (yearly) ----------------------------------
d = 1

# ----------------------------------- Case 3 -----------------------------------
# --------------------------------- (2 years) ----------------------------------
d = 2

# ----------------------------------- Case 4 -----------------------------------
# --------------------------------- (2 months) ---------------------------------
d = 1/6


mice_format$disc_time = mice_format$t1 - (mice_format$t1 %% d)
mice_format$obstrue = 0

newDat = NULL
num = 0
# Adding censored rows
for(i in unique(mice_format$ptnum)){
	
	current <- NULL
	subject <- mice_format[mice_format$ptnum==i,,drop=FALSE]
	
	#------------------------------------
	#censoredAges <- unique( c( min(subject$age), ceiling(min(subject$age)):max(subject$age)) )
	censoredAges <- seq(0, max(subject$t1), d)  # updates every 4 seconds
	for(t in censoredAges){
    tempRow = NULL
    d_floor_t = t - (t %% d)
		# If 't' corresponds to an observed age, then the next row will include the observed clinical visit data.
		if(t %in% subject$t1){	
			current <- rbind( current, subject[subject$disc_time==d_floor_t,]) 
		} else{
		
			# Create a CENSORED row for each subject at each INTEGER year of years.
			tempRow['t1'] <- t 
            tempRow['t2'] <- NA
			tempRow['state'] <- -99
            tempRow['delta'] = tempRow['theta'] = tempRow['alpha'] = tempRow['beta'] = NA
			tempRow['ptnum'] <- i
			tempRow['disc_time'] <- t
			tempRow['obstrue'] <- 1  
			
			current <- rbind( current, tempRow)
			
			# If 't' corresponds to an observed INTEGER years, then the subject was observed some time during this years.  According, the next row will include the observed clinical visit data.  Recall that integer years is simply the floor(years).
			if(t %in% subject$disc_time){ current <- rbind( current, subject[subject$disc_time==d_floor_t,]) }
		}

	}
	#------------------------------------
	
	newDat <- rbind( newDat, current)
  rownames(newDat) = NULL
	print(num)
	num <- num+1
}

mice_format = newDat
save(mice_format, file = 'Data_format/mice_format_sub_total_split_expm.rda')
