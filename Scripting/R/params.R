CreateDefaultParams <- function(param_mix_num, vcovar_num){
  init <- matrix(NA,ncol = 2,nrow = param_mix_num)
  init[,1] <- seq(.1,.9,length.out = param_mix_num)
  init[,2] <- 1 - init[,1]

  params_tran_week <- matrix(rep(c(0,0,0,0,0,0),param_mix_num),ncol = 6,byrow = T)
  params_tran_week[,1] <- seq(-2.9,-2.3,length.out = param_mix_num)
  params_tran_week[,2] <- seq(1.6,1.2,length.out = param_mix_num)
  params_tran_week[,3] <- seq(.3,.8,length.out = param_mix_num)
  params_tran_week[,4] <- seq(-1.8,-2.2,length.out = param_mix_num)
  params_tran_week[,5] <- seq(-1.6,-1.2,length.out = param_mix_num)
  params_tran_week[,6] <- seq(-.5,-.7,length.out = param_mix_num)

  params_tran_weekend <- matrix(rep(c(0,0,0,0,0,0),param_mix_num),ncol = 6,byrow = T)
  params_tran_weekend[,1] <- seq(-3.1,-2.5,length.out = param_mix_num)
  params_tran_weekend[,2] <- seq(1,.8,length.out = param_mix_num)
  params_tran_weekend[,3] <- seq(.7,.9,length.out = param_mix_num)
  params_tran_weekend[,4] <- seq(-1.8,-2.4,length.out = param_mix_num)
  params_tran_weekend[,5] <- seq(-1.3,-.9,length.out = param_mix_num)
  params_tran_weekend[,6] <- seq(-.9,-1.1,length.out = param_mix_num)

  params_tran_array <- array(NA,dim = c(param_mix_num,6,vcovar_num))
  params_tran_array[,,1] <- params_tran_week
  params_tran_array[,,2] <- params_tran_weekend

  emit_act_week <- array(NA, c(2,2,param_mix_num))
  emit_act_week[1,1,] <- seq(4,5,length.out = param_mix_num)
  emit_act_week[1,2,] <- seq(2,3,length.out = param_mix_num)
  emit_act_week[2,1,] <- seq(2,1,length.out = param_mix_num)
  emit_act_week[2,2,] <- seq(3,2,length.out = param_mix_num)

  emit_act_weekend <- array(NA, c(2,2,param_mix_num))
  emit_act_weekend[1,1,] <- seq(5,6,length.out = param_mix_num)
  emit_act_weekend[1,2,] <- seq(3,4,length.out = param_mix_num)
  emit_act_weekend[2,1,] <- seq(2,2,length.out = param_mix_num)
  emit_act_weekend[2,2,] <- seq(3,2,length.out = param_mix_num)

  emit_act <- array(NA, c(2,2,param_mix_num,2))
  emit_act[,,,1] <- emit_act_week
  emit_act[,,,2] <- emit_act_weekend

  emit_light_week <- array(NA, c(2,2,param_mix_num))
  emit_light_week[1,1,] <- seq(-2,-1,length.out = param_mix_num)
  emit_light_week[1,2,] <- seq(8,6,length.out = param_mix_num)
  emit_light_week[2,1,] <- seq(-19,-19,length.out = param_mix_num)
  emit_light_week[2,2,] <- seq(13,15,length.out = param_mix_num)

  emit_light_weekend <- array(NA, c(2,2,param_mix_num))
  emit_light_weekend[1,1,] <- seq(-3,-2,length.out = param_mix_num)
  emit_light_weekend[1,2,] <- seq(9,7,length.out = param_mix_num)
  emit_light_weekend[2,1,] <- seq(-21,-21,length.out = param_mix_num)
  emit_light_weekend[2,2,] <- seq(15,18,length.out = param_mix_num)

  emit_light <- array(NA, c(2,2,param_mix_num,2))
  emit_light[,,,1] <- emit_light_week
  emit_light[,,,2] <- emit_light_weekend

  corr_mat <- array(NA, c(param_mix_num,2,2))
  corr_mat[,1,1] <- seq(.1,.2,length.out = param_mix_num)
  corr_mat[,2,1] <- seq(.2,.3,length.out = param_mix_num)
  corr_mat[,1,2] <- seq(.3,.4,length.out = param_mix_num)
  corr_mat[,2,2] <- seq(.4,.5,length.out = param_mix_num)

  beta_vec <- seq(0,3,length.out = param_mix_num)
  beta_age <- 0.07

  nu <- c(0,seq(-.04,.05,length.out = (param_mix_num-1)))
  nu2 <- c(0,seq(.0005,-.0015,length.out = (param_mix_num-1)))
  nu_stat <- c(0,seq(.01,-.025,length.out = (param_mix_num-1)))
  nu2_stat <- c(0,seq(-.001,.002,length.out = (param_mix_num-1)))

  nu_mat <- matrix(NA,nrow = 4, ncol = param_mix_num)
  nu_mat[1,] <- nu
  nu_mat[2,] <- nu2
  nu_mat[3,] <- nu_stat
  nu_mat[4,] <- nu2_stat

  lambda_act_mat <- array(NA,dim = c(param_mix_num,2,2))
  lambda_act_mat[,1,1] <- seq(.01,.05,length.out = param_mix_num)
  lambda_act_mat[,2,1] <- seq(.3,.6,length.out = param_mix_num)
  lambda_act_mat[,1,2] <- seq(.01,.1,length.out = param_mix_num)
  lambda_act_mat[,2,2] <- seq(.3,.7,length.out = param_mix_num)

  lambda_light_mat <- array(NA,dim = c(param_mix_num,2,2))
  lambda_light_mat[,1,1] <- seq(.01,.15,length.out = param_mix_num)
  lambda_light_mat[,2,1] <- seq(.2,.6,length.out = param_mix_num)
  lambda_light_mat[,1,2] <- seq(.01,.2,length.out = param_mix_num)
  lambda_light_mat[,2,2] <- seq(.3,.8,length.out = param_mix_num)

  list(init = init,
       params_tran_array = params_tran_array,
       emit_act = emit_act,
       emit_light = emit_light,
       corr_mat = corr_mat,
       nu_mat = nu_mat,
       beta_vec = beta_vec,
       beta_age = beta_age,
       lambda_act_mat = lambda_act_mat,
       lambda_light_mat = lambda_light_mat)
}
