#Calculates mu1|2
CalcCondMean <- function(mu1,sig1,mu2,sig2,bivar_corr,obs2){
  return(mu1 + bivar_corr*(sig1/sig2)*(obs2-mu2))
}

#Calculates sig1|2
CalcCondSig <- function(sig1,bivar_corr){
  return(sig1*sqrt(1-bivar_corr^2))
}

#Calcuates case where both activity and light are below LoD
Case4 <- function(act_obs,mu_act,sig_act,mu_light,sig_light,bivar_corr,light_LOD){

  mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs)
  sig_light_cond <- CalcCondSig(sig_light,bivar_corr)

  lognorm_dens <- dnorm(act_obs,mu_act,sig_act) *
    pnorm(light_LOD,mu_light_cond,sig_light_cond)
  return(lognorm_dens)
}

#Calculates case 4 for given normal parameters and store in matrix for access later
CalcLintegralMat <- function(emit_act,emit_light,corr_mat,lod_act,lod_light){
  mix_num <- dim(emit_act)[3]
  if (is.na(mix_num)){mix_num <- 1}

  lintegral_mat <- array(NA,dim = c(mix_num,2,2))
  #j is week/weekend
  for (j in 1:2){
    for (i in 1:mix_num){

      lintegral_mat[i,1,j] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                          emit_act[1,1,i,j],emit_act[1,2,i,j],
                                          emit_light[1,1,i,j],emit_light[1,2,i,j],
                                          corr_mat[i,1,j],lod_light)[[1]])

      lintegral_mat[i,2,j] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                          emit_act[2,1,i,j],emit_act[2,2,i,j],
                                          emit_light[2,1,i,j],emit_light[2,2,i,j],
                                          corr_mat[i,2,j],lod_light)[[1]])
    }
  }

  #work on log scale so -9999 is effectively -Inf
  lintegral_mat[lintegral_mat == -Inf] <- -9999

  return(lintegral_mat)
}

#above parameters are for debugging
#calculates likelihood of emission dist
#used in direct optimization
EmitLogLike <- function(act,light,mu_act,sig_act,mu_light,sig_light,bivar_corr,lod_act,lod_light,vcovar_mat,vcovar_ind,weights_mat){

  #lower should theoretically be -Inf, but had some divergence issues
  lb <- mu_act - 5*sig_act
  lb <- min(lb,-10)
  lintegral <- log(integrate(Case4,lower = lb,upper = lod_act,
                             mu_act,
                             sig_act,
                             mu_light,
                             sig_light,
                             bivar_corr,lod_light)[[1]])

  if (lintegral == -Inf){
    lintegral <- -9999
  }

  vcovar_vec <- as.vector(vcovar_mat)
  vcovar_vec_indicator <- vcovar_vec == vcovar_ind

  act_vec <- as.vector(act)[vcovar_vec_indicator]
  light_vec <- as.vector(light)[vcovar_vec_indicator]

  log_like <- logClassificationCTobit(act_vec,light_vec,
                                 mu_act,
                                 sig_act,
                                 mu_light,
                                 sig_light,
                                 lod_act,lod_light,bivar_corr,lintegral)

  log_like[log_like == -Inf] <- -9999

  return(-sum(log_like * weights_mat[vcovar_vec_indicator]))
}

#optimizes activity mean
#all emission dist param calculated this way
#comment out which parameter currently being optimized
CalcActMean <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){

  mu_act <- optimize(EmitLogLike, c(-10,10), act = act, light = light,
                     # mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcActSig <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){

  mu_act <- optimize(EmitLogLike, c(0.1,10), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     # sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcLightMean <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){

  mu_act <- optimize(EmitLogLike, c(-30,10), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     # mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum

  return(mu_act)
}

CalcLightSig <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){

  mu_act <- optimize(EmitLogLike, c(0.01,20), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     #sig_light = emit_light[mc_state,2,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcBivarCorr <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){

  mu_act <- optimize(EmitLogLike, c(-.999,.999), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     # bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

#Highest level for optimizing emission dist
#takes function as input, easier to process this way
UpdateNorm <- function(FUN,mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat){
  opt_param_mat <- matrix(0,mix_num,2)

  if(mc_state == 1){weights_array <- weights_array_wake}
  if(mc_state == 2){weights_array <- weights_array_sleep}

  for (re_ind in 1:dim(emit_act)[3]){
    for(vcovar_ind in 1:2){
      opt_param_mat[re_ind,vcovar_ind] <- FUN(mc_state,vcovar_ind,
                                              act,light,
                                              emit_act,emit_light,
                                              corr_mat,lod_act,lod_light,
                                              weights_array,re_ind,vcovar_mat)
    }
  }

  return(opt_param_mat)
}
