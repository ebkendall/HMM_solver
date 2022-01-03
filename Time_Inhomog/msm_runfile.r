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

package.install("msm")

args = commandArgs(TRUE)
folder = as.numeric(args[1])

seedInd = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seedInd)

data_files <- c("Month", "Year", "YearTwo")

Dir <- paste0("../DataGeneration/DataOut/", data_files[folder], "/")
load(paste0(Dir,'cavData',seedInd,'.rda'))

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

Output_msm <- msm(state ~ years, subject=ptnum, data=cavData, qmatrix=qmat, covariates= ~ 1 + years + sex,
				  center=FALSE, covinits=list(sex=c(0,0,0,0,0)), obstrue=obstrue,
				  ematrix=emat, initprobs=c(1, exp(-4.5), exp(-5), 0), est.initprobs=TRUE, deathexact=4,
				  censor=99, censor.states=1:3, method='BFGS', control=list(fnscale=4000, maxit=10000))

save(Output_msm,file=paste0(Dir,'Output_msm',seedInd,'.rda'))
