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

# This script contains the code for the mcmc and its helper functions

package.install("mvtnorm")
package.install("foreach")
package.install("doParallel")
package.install("msm")
package.install("deSolve")
package.install("expm")
package.install("foreach")
package.install("doParallel")


# library(mvtnorm, quietly=T)
# library(foreach, quietly=T)
# library(doParallel, quietly=T)
#
# library("msm")
# library("deSolve")
# library("expm")
# library("foreach")
# library("doParallel")

Q <- function(x_ik,beta){

  betaMat = matrix(beta, ncol = 3, byrow = F) # ncol is 3 for time inhomogenous
  q1  = exp( c(1,x_ik) %*% betaMat[1,] )  # Transition from state 1 to state 2.
  q2  = exp( c(1,x_ik) %*% betaMat[2,] )  # Transition from state 2 to state 3.
  q3  = exp( c(1,x_ik) %*% betaMat[3,] )  # Transition from state 1 to death.
  q4  = exp( c(1,x_ik) %*% betaMat[4,] )  # Transition from state 2 to death.
  q5  = exp( c(1,x_ik) %*% betaMat[5,] )  # Transition from state 3 to death.

  qmat = matrix(c( 0,q1, 0,q2,
                   0, 0,q3,q4,
                   0, 0, 0,q5,
                   0, 0, 0, 0),nrow=4,byrow=TRUE)
  diag(qmat) = -rowSums(qmat)

  return(qmat)
}

fn_log_post <- function(pars, prior_par, par_index, x, y, t, id) {

  init_logit = c( 1, exp(pars[par_index$pi_logit][1]), exp(pars[par_index$pi_logit][2]), 0)
  init = init_logit / sum(init_logit)
  resp_fnc = matrix(c(1, exp(pars[par_index$misclass][1]), 0, 0,
                      exp(pars[par_index$misclass][2]), 1, exp(pars[par_index$misclass][3]), 0,
                      0, exp(pars[par_index$misclass][4]),1, 0,
                      0,   0,   0, 1), ncol=4, byrow=TRUE)

  resp_fnc = resp_fnc / rowSums(resp_fnc)

  beta <- pars[par_index$beta]
  p_ic <- c(p1=1,p2=0,p3=0,p4=0,p5=1,p6=0,p7=0,p8=1,p9=0)

  log_total_val = foreach(i=unique(id), .combine='+', .export = "Q", .packages = "expm") %dopar% {

    y_i = y[id == i]
    x_i = x[id == i,,drop = F]
    t_i = t[id == i]
    Q_i = Q(x_i[1,], beta) # x_i is a 2D vector in time inhomogeneous

    f_i = init %*% diag(resp_fnc[, y_i[1]])
	log_norm = 0

    for(k in 2:length(t_i)) {

      P = expm((t_i[k] - t_i[k-1]) * Q_i)

  	  if(y_i[k] < 4) {
          val = f_i %*% P %*% diag(resp_fnc[, y_i[k]])
        } else {
          val = f_i %*% P %*% Q(x_i[k,], beta) %*% diag(resp_fnc[, y_i[k]])
        }

      norm_val = sqrt(sum(val^2))
  	  f_i = val / norm_val
  	  log_norm = log_norm + log(norm_val)
    }

    return(log(sum(f_i)) + log_norm)
  }

    mean = prior_par$prior_mean
    sd = diag(prior_par$prior_sd)
    log_prior_dens = dmvnorm( x=pars, mean=mean, sigma=sd, log=T)

  return(log_total_val + log_prior_dens)

}

# -----------------------------------------------------------------------------
# The mcmc routine for samping the parameters
# -----------------------------------------------------------------------------
mcmc_routine = function( y, x, t, id, init_par, prior_par, par_index,
                         steps, burnin, n_cores){

  cl <- makeCluster(n_cores, outfile="")
  registerDoParallel(cl)

  pars = init_par
  n = length(y)
  n_par = length(pars)
  chain = matrix( 0, steps, n_par)

  group = list(c(par_index$beta, par_index$misclass, par_index$pi_logit))
  n_group = length(group)

  pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))
  pscale = rep( .0001, n_group)
  accept = rep( 0, n_group)

  # Evaluate the log_post of the initial pars
  log_post_prev = fn_log_post( pars, prior_par, par_index, x, y, t, id)
  if(!is.finite(log_post_prev)){
    print("Infinite log-posterior; choose better initial parameters")
    break
  }

  # Begin the MCMC algorithm --------------------------------------------------
  chain[1,] = pars
  for(ttt in 2:steps){
    for(j in 1:n_group){

      # Propose an update
      ind_j = group[[j]]
      proposal = pars
      proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])

      # Compute the log density for the proposal
      log_post = fn_log_post(proposal, prior_par, par_index, x, y, t, id)
      # print("Likelihood Evaluation:")
      # print(log_post)

      # Only propose valid parameters during the burnin period
      if(ttt < burnin){
        while(!is.finite(log_post)){
          print('bad proposal')
          proposal = pars
          proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                     sigma=pcov[[j]]*pscale[j])
          log_post = fn_log_post(proposal, prior_par, par_index, x, y, t, id)
        }
      }

      # Evaluate the Metropolis-Hastings ratio
      if( log_post - log_post_prev > log(runif(1,0,1)) ){
        log_post_prev = log_post
        pars[ind_j] = proposal[ind_j]
        accept[j] = accept[j] +1
      }
      chain[ttt,ind_j] = pars[ind_j]
      # print("Parameters Accepted")
      # print(pars)

      # Proposal tuning scheme ------------------------------------------------
      if(ttt < burnin){
        # During the burnin period, update the proposal covariance in each step
        # to capture the relationships within the parameters vectors for each
        # transition.  This helps with mixing.
        if(ttt == 100)  pscale[j] = 1

        if(100 <= ttt & ttt <= 2000){
          temp_chain = chain[1:ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])

        } else if(2000 < ttt){
          temp_chain = chain[(ttt-2000):ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
        }
        if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )

        # Tune the proposal covariance for each transition to achieve
        # reasonable acceptance ratios.
        if(ttt %% 30 == 0){
          if(ttt %% 480 == 0){
            accept[j] = 0

          } else if( accept[j] / (ttt %% 480) < .4 ){ #change to 0.3
            pscale[j] = (.75^2)*pscale[j]

          } else if( accept[j] / (ttt %% 480) > .5 ){ #change to 0.4
            pscale[j] = (1.25^2)*pscale[j]
          }
        }
      }
      # -----------------------------------------------------------------------
    }
    # Restart the acceptance ratio at burnin.
    if(ttt == burnin)  accept = rep( 0, n_group)

    if(ttt%%1==0)  cat('--->',ttt,'\n')
  }
  # ---------------------------------------------------------------------------

  stopCluster(cl)
  print(accept/(steps-burnin))
  return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
               pscale=pscale))
}
# -----------------------------------------------------------------------------
