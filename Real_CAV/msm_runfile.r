library("msm")

seedInd = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seedInd)

data_files <- c("Year_msm", "YearTwo_msm", "Month_msm")
data_names <- c('Data/cavData_year.rda', 
				'Data/cavData_twoYear.rda', 
				'Data/cavData_month.rda')

for (folder in 1:3) {

	print(folder)

	load(data_names[folder])

	qmat <- matrix(c( 0,exp(-2),      0,exp(-2),
	                  0,      0,exp(-2),exp(-2),
	                  0,      0,      0,exp(-2),
	                  0,      0,      0,      0), ncol=4, byrow=TRUE)
	dimnames(qmat) <- list( c('Well', 'Mild','Severe','Death'), c('Well', 'Mild','Severe','Death'))

	#----------------------------------------------------------------------------------------------------------------
	# Run the msm implementation ------------------------------------------------------------------------------------
	#----------------------------------------------------------------------------------------------------------------

	emat = matrix(c(    1, exp(-3),       0, 0,
	       		  exp(-3),       1, exp(-3), 0,
						0, exp(-3),       1, 0,
				        0,       0,       0, 1), ncol=4, byrow=TRUE)
	emat = emat / rowSums(emat)
	dimnames(emat) <- list( c('Well','Mild','Severe','Death'), c('Well','Mild','Severe','Death'))

	Output_msm <- msm(state ~ years, subject=ptnum, data=cavData, qmatrix=qmat, covariates= ~ 1 + disc_time + sex, 
					  center=FALSE, covinits=list(disc_time=c(0,0,0,0,0),sex=c(0,0,0,0,0)), obstrue=obstrue, 
					  ematrix=emat, initprobs=c(1, exp(-4.5), exp(-5), 0), est.initprobs=TRUE, deathexact=4, 
					  censor=99, censor.states=1:3, method='BFGS', control=list(fnscale=4000, maxit=10000))   
	
	save(Output_msm,file=paste0('Model_out/', data_files[folder] ,'/Output_msm',seedInd,'.rda'))

}