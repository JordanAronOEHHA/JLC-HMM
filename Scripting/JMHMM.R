################## Intro ################## 

library(Rcpp)
library(RcppArmadillo)
library(matrixStats)
library(MASS)
library(survival)
library(dplyr)
library(numDeriv)
library(Matrix)
library(Hmisc)
library(survex)
library(tidyverse)

source_jmhmm_module <- function(module_file){
  candidate_paths <- c(file.path("Scripting","R",module_file),
                       file.path("R",module_file),
                       file.path("..","Rcode","R",module_file))
  for (candidate_path in candidate_paths){
    if (file.exists(candidate_path)){
      source(candidate_path, local = parent.frame())
      return(invisible(candidate_path))
    }
  }
  stop(paste("Could not find module:",module_file))
}

for (module_file in c("constants.R","saved_results.R","validation.R",
                      "settings.R","params.R","transitions.R",
                      "emissions_tobit.R","forward_backward.R",
                      "data_simulation.R", "helpers.R",
                      "survival.R","diagnostics.R")){
  source_jmhmm_module(module_file)
}

settings <- build_settings()
settings$model_name <- build_model_name(settings)

cli_args <- settings$command_args
sim_num <- settings$sim_num

# Compatibility aliases: most of the legacy script still reads these names.
list2env(settings, envir = environment())

print("Command line arguments:")
print(cli_args)
print("Run settings:")
print(settings)

if (set_seed){set.seed(sim_num)}

print(paste("Sim Seed:",sim_num,"Fit HMM Num:",fit_mix_num,"True HMM Num:",true_mix_num))
print(model_name)

################## EM Setup ################## 

readCpp( "Scripting/cFunctions.cpp" )
readCpp( "../Rcode/cFunctions.cpp" )

###### True Settings ###### 

#Sets up simulation sizing
sim_config <- SIM_SCENARIOS[[as.character(sim_size)]]
if (is.null(sim_config)){
  stop(paste("Unknown sim_scenario:",sim_size))
}
day_length <- period_len * sim_config$days
num_of_people <- sim_config$num_people
missing_perc <- sim_config$missing_perc



true_param_list <- CreateDefaultParams(true_mix_num, vcovar_num)
validate_param_list(true_param_list,true_mix_num,vcovar_num,"true_param_list")
init_true <- true_param_list$init
params_tran_array_true <- true_param_list$params_tran_array
emit_act_true <- true_param_list$emit_act
emit_light_true <- true_param_list$emit_light
corr_mat_true <- true_param_list$corr_mat
nu_mat_true <- true_param_list$nu_mat
beta_vec_true <- true_param_list$beta_vec
beta_age_true <- true_param_list$beta_age
lambda_act_mat_true <- true_param_list$lambda_act_mat
lambda_light_mat_true <- true_param_list$lambda_light_mat

fit_param_list <- CreateDefaultParams(fit_mix_num, vcovar_num)
validate_param_list(fit_param_list,fit_mix_num,vcovar_num,"fit_param_list")
init_start <- fit_param_list$init
params_tran_array_start <- fit_param_list$params_tran_array
emit_act_start <- fit_param_list$emit_act
emit_light_start <- fit_param_list$emit_light
corr_mat_start <- fit_param_list$corr_mat
nu_mat_start <- fit_param_list$nu_mat
beta_vec_start <- fit_param_list$beta_vec
lambda_act_mat_start <- fit_param_list$lambda_act_mat
lambda_light_mat_start <- fit_param_list$lambda_light_mat
#loads data in for a hot start
if (load_data){
  model_name_loadin <- "JMHMM"
  load_mix_num <- if (real_data){fit_mix_num} else {true_mix_num}
  folder_name <- paste(load_mix_num)
  if (incl_surv == MODEL_TYPE_CODES[["two_stage"]]){
    model_name_loadin <- paste0(model_name_loadin,"NoSurv")
    folder_name <- paste0("NS",folder_name)
  }


  legacy_model_name_loadin <- paste0(model_name_loadin,"Mix",load_mix_num,"Seed",".rda")
  model_name_loadin <- paste0(model_name_loadin,"FitMix",load_mix_num,"Seed",".rda")
  print(paste("Loading",model_name_loadin))
  setwd("Data")
  if (!file.exists(model_name_loadin)){
    model_name_loadin <- legacy_model_name_loadin
    print(paste("Loading legacy",model_name_loadin))
  }
  load(model_name_loadin)
  setwd("..")
  validate_saved_results(to_save,required_sections = c("est_params"),
                         source_name = model_name_loadin)

  # # model_name_loadin <- paste0(model_name_loadin,"Mix",mix_num,"Seed",sim_num,".rda")
  # # print(paste("Loading",model_name_loadin))
  #
  # model_name_loadin <- paste0("Inter",model_name)
  # print(paste("Loading",model_name_loadin))
  # # setwd(folder_name)
  # load(model_name_loadin)
  # setwd("..")
 
  loaded_est_params <- get_saved_section(to_save,"est_params")
  loaded_init <- get_saved_param(loaded_est_params,"init")
  loaded_params_tran_array <- get_saved_param(loaded_est_params,"params_tran_array")
  loaded_emit_act <- get_saved_param(loaded_est_params,"emit_act")
  loaded_emit_light <- get_saved_param(loaded_est_params,"emit_light")
  loaded_corr_mat <- get_saved_param(loaded_est_params,"corr_mat")
  loaded_nu_mat <- get_saved_param(loaded_est_params,"nu_mat")
  loaded_beta_vec <- get_saved_param(loaded_est_params,"beta_vec")
  loaded_lambda_act_mat <- lambda_act_mat_true
  loaded_lambda_light_mat <- lambda_light_mat_true
  surv_coef_true <- get_saved_param(loaded_est_params,"surv_coef")

  re_prob_true <- get_saved_param(loaded_est_params,"re_prob")
  re_prob <- re_prob_true

  loaded_lambda_act_mat <- get_saved_param(loaded_est_params,"lambda_act_mat",
                                           default = loaded_lambda_act_mat,
                                           required = FALSE)
  loaded_lambda_light_mat <- get_saved_param(loaded_est_params,"lambda_light_mat",
                                             default = loaded_lambda_light_mat,
                                             required = FALSE)

  if (!real_data){
    init_true <- loaded_init
    params_tran_array_true <- loaded_params_tran_array
    emit_act_true <- loaded_emit_act
    emit_light_true <- loaded_emit_light
    corr_mat_true <- loaded_corr_mat
    nu_mat_true <- loaded_nu_mat
    beta_vec_true <- loaded_beta_vec
    lambda_act_mat_true <- loaded_lambda_act_mat
    lambda_light_mat_true <- loaded_lambda_light_mat
  }

  if (dim(loaded_init)[1] == fit_mix_num){
    init_start <- loaded_init
    params_tran_array_start <- loaded_params_tran_array
    emit_act_start <- loaded_emit_act
    emit_light_start <- loaded_emit_light
    corr_mat_start <- loaded_corr_mat
    nu_mat_start <- loaded_nu_mat
    beta_vec_start <- loaded_beta_vec
    lambda_act_mat_start <- loaded_lambda_act_mat
    lambda_light_mat_start <- loaded_lambda_light_mat
  }
  
}

###### Simulate Data ###### 
if (!real_data){
  lod_act_true <- -5.809153
  lod_light_true <- -1.560658
  
  lod_act <- lod_act_true
  lod_light <- lod_light_true
  
  beta_covar_sim <- c(0,.6,-.5)
  
  simulated_hmm <- SimulateHMM(day_length,num_of_people,
                               init=init_true,params_tran_array = params_tran_array_true,
                               emit_act = emit_act_true,emit_light = emit_light_true,
                               corr_mat = corr_mat_true,
                               lod_act = lod_act_true,lod_light = lod_light_true,
                               nu_mat = nu_mat_true,
                               beta_age_true = beta_age_true,beta_covar_sim = beta_covar_sim,
                               missing_perc = missing_perc, beta_vec_true = beta_vec_true,
                               lambda_act_mat = lambda_act_mat_true,lambda_light_mat = lambda_light_mat_true,
                               true_mix_num = true_mix_num)
  mc <- simulated_hmm$mc
  act <- simulated_hmm$act
  light <- simulated_hmm$light
  mixture_mat <- simulated_hmm$mixture_mat
  age_vec <- simulated_hmm$age_vec
  nu_covar_mat <- simulated_hmm$nu_covar_mat
  vcovar_mat <-  simulated_hmm$vcovar_mat
  surv_list <- simulated_hmm$survival
  surv_covar_sim <- simulated_hmm$surv_covar_sim
  
  id_sim <- cbind(age_vec,surv_covar_sim-1)
  surv_covar <- list(age_vec,Vec2Mat(surv_covar_sim))
  surv_coef <- list(beta_age_true,beta_covar_sim)
  surv_coef_true <- surv_coef
  combined_covar_mat <- matrix(surv_covar_sim-1,nrow = num_of_people)
  combined_covar_mat <- as.factor(combined_covar_mat)
  
  surv_time <- surv_list$time
  surv_event <- surv_list$event
  
  #in simulated data sample weights are set to 0
  log_sweights_vec <- numeric(dim(act)[2])
  
  
} 
###### Read in Data ###### 

if (real_data) {
  #loads in NHANES
  setwd("Data/")
  load("NHANES_2011_2012_2013_2014.rda")
  nhanes1 <- NHANES_mort_list[[1]] %>% filter(eligstat == 1)
  nhanes2 <- NHANES_mort_list[[2]] %>% filter(eligstat == 1)
  lmf_data <- rbind(nhanes1,nhanes2)
  
  #depending on period len, loads in different data
  if (period_len == HOURLY_PERIODS_PER_DAY){
    load("Wavedata24_G.rda")
    load("Wavedata24_H.rda")
  } else if(period_len == DEFAULT_PERIODS_PER_DAY){
    load("Wavedata_G.rda")
    load("Wavedata_H.rda")
  } else if(period_len == MINUTE_PERIODS_PER_DAY){
    load("Wavedata1440_G.rda")
    load("Wavedata1440_H.rda")
  }
  
  setwd("..")
  
  #preps and combines 2 wakes of act data
  act_G <- wave_data_G[[1]]
  act_H <- wave_data_H[[1]]
  act <- rbind(act_G,act_H)
  act <- t(act[,-1])
  act0 <- act == 0
  act <- log(act)
  lod_act <- min(act[act!=-Inf],na.rm = T) - LOD_OFFSET
  act[act0] <- lod_act
  
  
  #preps and combines 2 wakes of light data
  light_G <- wave_data_G[[2]]
  light_H <- wave_data_H[[2]]
  light <- rbind(light_G,light_H)
  light <- t(light[,-1])
  light0 <- light == 0
  light <- log(light)
  lod_light <- min(light[light!=-Inf],na.rm = T) - LOD_OFFSET
  light[light0] <- lod_light
  
  id_G <- wave_data_G[[3]]
  id_H <- wave_data_H[[3]]
  id <- rbind(id_G,id_H)
  
  mims_G <- wave_data_G[[4]]
  mims_H <- wave_data_H[[4]]
  mims <- rbind(mims_G,mims_H)
  
  #matches actigraphy data to public mortality data
  seqn_com_id <- id$SEQN %in% lmf_data$seqn
  seqn_com_lmf <- lmf_data$seqn %in% id$SEQN
  
  id <- id[seqn_com_id,]
  act <- act[,seqn_com_id]
  light <- light[,seqn_com_id]
  mims <- mims[seqn_com_id,]
  
  lmf_data <- lmf_data[seqn_com_lmf,]
  
  #sanity check
  if (sum(id$SEQN - lmf_data$seqn) != 0){print("LMF NOT LINKED CORRECTLY")}
  
  #sample weights for 2 waves
  log_sweights_vec <- log(id$sweights/NHANES_NUM_WAVES)
  
  id <- id %>% mutate(age_disc = case_when(age <=30 ~ 1,
                                           age <=50 & age > 30 ~ 2,
                                           age <=65 & age > 50 ~ 3,
                                           age > 65 ~ 4))
  
  id <- id %>% mutate(pov_disc = floor(poverty)+1)
  
  id$modact <- id$modact - 1
  
  surv_event <- lmf_data$mortstat
  surv_time <- lmf_data$permth_exm
  
  
    
  
  #resample data for variance estimation
  if (bootstrap){
    boot_inds <- sample(dim(act)[2],dim(act)[2],T)
    act <- act[,boot_inds]
    light <- light[,boot_inds]
    id <- id[boot_inds,]
    surv_event <- surv_event[boot_inds]
    surv_time <- surv_time[boot_inds]
    
  }
  
  #saves original data before LOCV
  act_old <- act
  light_old <- light
  
  #cross validation
  if (leave_out){
    
    setwd("Data")
    #previously calculated who is left in/out for each seed
    load("LeaveOutMat.rda")
    setwd("..")
    leave_out_inds <- leave_out_mat[sim_num,]
    leave_out_inds <- leave_out_inds[!is.na(leave_out_inds)]
    
    
    act_old <- act
    light_old <- light
    id_old <- id
    surv_event_old <- surv_event
    surv_time_old <- surv_time
    log_sweights_vec_old <- log_sweights_vec
    
    first_day_vec_old <- as.numeric(id_old$PAXDAYWM)
    vcovar_mat_old <- sapply(first_day_vec_old,FirstDay2WeekInd)
    
    surv_covar_old <- list(id_old$age,
                       Vec2Mat(id_old$gender+1),
                       Vec2Mat(id_old$race+1),
                       Vec2Mat(id_old$overall_health+1),
                       Vec2Mat(id_old$education+1),
                       Vec2Mat(id_old$bmi_disc+1),
                       Vec2Mat(id_old$diabetes+1),
                       Vec2Mat(id_old$CHD+1),
                       Vec2Mat(id_old$CHF+1),
                       Vec2Mat(id_old$heart_attack+1),
                       Vec2Mat(id_old$stroke+1),
                       Vec2Mat(id_old$alcohol+1),
                       Vec2Mat(id_old$smoking+1),
                       Vec2Mat(id_old$phyfunc+1))
    
    age_vec_old <-id_old$age
    statact_vec_old <- id_old$statact
    nu_covar_mat_old <- cbind(age_vec_old/10,(age_vec_old/10)^2,statact_vec_old,statact_vec_old^2)
    
    act <- act[,-c(leave_out_inds)]
    light <- light[,-c(leave_out_inds)]
    id <- id[-c(leave_out_inds),]
    surv_event <- surv_event[-c(leave_out_inds)]
    surv_time <- surv_time[-c(leave_out_inds)]
    
    log_sweights_vec <- log(id$sweights/NHANES_NUM_WAVES)
    
  }
  
  first_day_vec <- as.numeric(id$PAXDAYWM)
  vcovar_mat <- sapply(first_day_vec,FirstDay2WeekInd)
  
  #if single day can reduce memory used
  if (single_day != 0){
    single_day_mat <- sapply(first_day_vec,FirstDay2SingleDay,target_day = single_day)
    new_act <- matrix(NA,period_len,dim(act)[2])
    new_light <- matrix(NA,period_len,dim(light)[2])
    vcovar_mat <- matrix(0,period_len,dim(light)[2])
    
    for (i in 1:dim(act)[2]){
      new_act[,i] <- act[,i][single_day_mat[,i]==1]
      new_light[,i] <- light[,i][single_day_mat[,i]==1]
    }
    
    act <- new_act
    light <- new_light
  }
  
  day_length <- dim(act)[1]
  num_of_people <- dim(act)[2]
  
  age_vec <-id$age
  modact_vec <- id$modact
  statact_vec <- id$statact
  nu_covar_mat <- cbind(age_vec/10,(age_vec/10)^2,statact_vec,statact_vec^2)

  #sets of sociodemo covar list
  surv_covar <- list(id$age,
                     Vec2Mat(id$gender+1),
                     Vec2Mat(id$race+1),
                     Vec2Mat(id$overall_health+1),
                     Vec2Mat(id$education+1),
                     Vec2Mat(id$bmi_disc+1),
                     Vec2Mat(id$diabetes+1),
                     Vec2Mat(id$CHD+1),
                     Vec2Mat(id$CHF+1),
                     Vec2Mat(id$heart_attack+1),
                     Vec2Mat(id$stroke+1),
                     Vec2Mat(id$alcohol+1),
                     Vec2Mat(id$smoking+1),
                     Vec2Mat(id$phyfunc+1))

  #initializes survival covariates
  if (!load_data){
    surv_coef_true <- lapply(surv_covar[-1],SurvCovar2Coef)
    surv_coef_true <- append(list(.05),surv_coef_true)
  }
  surv_coef_len <- unlist(lapply(surv_coef_true,length))
  surv_coef <- surv_coef_true
  
  #all sociodemo covar values
  combined_covar_mat <- id %>% dplyr::select(gender,race,overall_health,education,bmi_disc,diabetes,
                                             race,CHD,CHF,heart_attack,stroke,alcohol,smoking,phyfunc)
  combined_covar_mat <- lapply(combined_covar_mat, factor)
  
  if (weekend_only){
    #Only weekend data
    
    act[vcovar_mat == 0] <- NA
    light[vcovar_mat == 0] <- NA
  }
    
}

###### Initial Settings ###### 

##########
#if doing cv, load in full data values for hot start
if (leave_out){
  
  model_name_loadin <- "JMHMM"
  if (incl_surv == MODEL_TYPE_CODES[["two_stage"]]){model_name_loadin <- paste0(model_name_loadin,"NoSurv")}
  legacy_model_name_loadin <- paste0(model_name_loadin,"Mix",fit_mix_num,"Seed",".rda")
  model_name_loadin <- paste0(model_name_loadin,"FitMix",fit_mix_num,"Seed",".rda")

  print(paste("Loading",model_name_loadin))
  setwd("Data")
  if (!file.exists(model_name_loadin)){
    model_name_loadin <- legacy_model_name_loadin
    print(paste("Loading legacy",model_name_loadin))
  }
  load(model_name_loadin)
  setwd("..")
  validate_saved_results(to_save,required_sections = c("est_params"),
                         source_name = model_name_loadin)
  full_data_est_params <- get_saved_section(to_save,"est_params")
  full_data_re_prob <- get_saved_param(full_data_est_params,"re_prob")
  mix_assignment_true <- apply(full_data_re_prob,1,which.max)
  post_decode_collapsed_true <- get_saved_param(full_data_est_params,"post_decode")[,leave_out_inds]
  
  
  # used if loading non-standard data
  model_name_loadin <- "JMHMMLeaveOut"
  if (incl_surv == MODEL_TYPE_CODES[["joint"]]){
    # foldername <- paste0("LO",mix_num)
  } else {
    model_name_loadin <- paste0(model_name_loadin,"NoSurv")
    # foldername <- paste0("LONS",mix_num)
  }
  legacy_model_name_loadin <- paste0(model_name_loadin,"Mix",fit_mix_num,"Seed",sim_num,"len96.rda")
  model_name_loadin <- paste0(model_name_loadin,"FitMix",fit_mix_num,"Seed",sim_num,"len96.rda")

  print(paste("Loading",model_name_loadin))
  setwd("LO")
  setwd(paste0(mix_num))
  if (!file.exists(model_name_loadin)){
    model_name_loadin <- legacy_model_name_loadin
    print(paste("Loading legacy",model_name_loadin))
  }
  load(model_name_loadin)
  setwd("..")
  setwd("..")
  validate_saved_results(to_save,required_sections = c("est_params"),
                         source_name = model_name_loadin)


  leave_out_est_params <- get_saved_section(to_save,"est_params")
  init_start <- get_saved_param(leave_out_est_params,"init")
  params_tran_array_start <- get_saved_param(leave_out_est_params,"params_tran_array")
  emit_act_start <- get_saved_param(leave_out_est_params,"emit_act")
  emit_light_start <- get_saved_param(leave_out_est_params,"emit_light")
  corr_mat_start <- get_saved_param(leave_out_est_params,"corr_mat")
  nu_mat_start <- get_saved_param(leave_out_est_params,"nu_mat")
  pi_l_true <- CalcPi(nu_mat_start,nu_covar_mat)
  beta_vec_start <- get_saved_param(leave_out_est_params,"beta_vec")
  surv_coef_true <- get_saved_param(leave_out_est_params,"surv_coef")
  re_prob <- get_saved_param(leave_out_est_params,"re_prob")

  lambda_act_mat_start <- get_saved_param(leave_out_est_params,"lambda_act_mat")
  lambda_light_mat_start <- get_saved_param(leave_out_est_params,"lambda_light_mat")
  
}

################## EM ##################

##### randomize starting parameters
# init <- matrix(rep(.5,mix_num*2),ncol = 2)
init <- init_start

params_tran_array <- params_tran_array_start + runif(unlist(length(params_tran_array_start)),-randomize_init*2,randomize_init*2)


emit_act <- emit_act_start + runif(length(unlist(emit_act_start)),-randomize_init,randomize_init)
emit_act[,2,,] <- abs(emit_act[,2,,])

emit_light <- emit_light_start + runif(length(unlist(emit_light_start)),-randomize_init*2,randomize_init*2)
emit_light[,2,,] <- abs(emit_light[,2,,])

#makes sure correlation makes sense
corr_mat <- corr_mat_start + runif(length(unlist(corr_mat_start)),-randomize_init/5,randomize_init/5)
corr_mat[corr_mat>.99] <- .99
corr_mat[corr_mat<-.99] <- -.99

#makes sure first val is always reference
beta_vec <- beta_vec_start + runif(mix_num,-randomize_init,randomize_init)
beta_vec[1] <- 0

for (i in 1:length(surv_coef)){
  surv_coef[[i]] <-surv_coef_true[[i]]  +  runif(length(surv_coef_true[[i]]),-randomize_init/10,randomize_init/10)
  if (length(surv_coef[[i]]) != 1){surv_coef[[i]][1] <- 0} 
}
surv_coef[[1]] <- surv_coef_true[[1]] + runif(1,-randomize_init/100,randomize_init/100)

#dont randomize as these are very sensitive
nu_mat <- nu_mat_start
lambda_act_mat <- lambda_act_mat_start
lambda_light_mat <- lambda_light_mat_start

start_params <- make_start_param_list(init = init,
                                      params_tran_array = params_tran_array,
                                      emit_act = emit_act,
                                      emit_light = emit_light,
                                      corr_mat = corr_mat,
                                      nu_mat = nu_mat,
                                      beta_vec = beta_vec,
                                      surv_coef = surv_coef,
                                      lambda_act_mat = lambda_act_mat,
                                      lambda_light_mat = lambda_light_mat)
validate_param_list(start_params,fit_mix_num,vcovar_num,"start_params")

time_vec <- c()
pi_l <- CalcPi(nu_mat,nu_covar_mat)

#sets some controls so matrix sizing lines up
if (!leave_out & !load_data){re_prob <- pi_l}
if (load_data & !real_data){re_prob <- pi_l}
if (load_data & period_len != 96){re_prob <- pi_l}
if (!is.null(dim(re_prob)) && ncol(re_prob) != fit_mix_num){re_prob <- pi_l}
if (is.null(dim(re_prob))){re_prob <- matrix(re_prob,ncol = 1)}

validate_hmm_data(act,light,vcovar_mat)
validate_survival_inputs(surv_time,surv_event,surv_covar,num_of_people)
survival_context <- make_survival_context(surv_time,surv_event,surv_covar,
                                          re_prob,fit_mix_num)

surv_coef_len <- unlist(lapply(surv_coef,length))
surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)

if(!incl_light){
  light <- matrix(NA,dim(light)[1],dim(light)[2])
  light_old <- matrix(NA,dim(light_old)[1],dim(light_old)[2])
}
if(!incl_act){
  act <- matrix(NA,dim(act)[1],dim(act)[2])
  act_old <- matrix(NA,dim(act_old)[1],dim(act_old)[2])
}

#calculates baseline hazards
bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,survival_context$re_prob,
                      surv_covar_risk_vec,survival_context$surv_event,
                      survival_context$surv_time,survival_context$surv_covar)
bline_vec <- bhaz_vec[[1]]
cbline_vec <- bhaz_vec[[2]]

#caluclates case4 probabilities ahead of time and transition list
lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)

print("Pre Alpha")
alpha <- Forward(act = act,light = light,
         init = init,tran_list = tran_list,
         emit_act = emit_act,emit_light = emit_light,
         lod_act = lod_act, lod_light = lod_light, 
         corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef,surv_covar_risk_vec = surv_covar_risk_vec,
         event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
         lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec,
         surv_covar = surv_covar, vcovar_mat = vcovar_mat,
         lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit,incl_surv = incl_surv)

beta <- Backward(act = act,light = light, tran_list = tran_list,
                 emit_act = emit_act,emit_light = emit_light,
                  lod_act = lod_act, lod_light =  lod_light, 
                  corr_mat = corr_mat,lintegral_mat = lintegral_mat,vcovar_mat = vcovar_mat,
                  lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit)
         
print("Post Beta")
new_likelihood <- CalcLikelihood(alpha,pi_l)
if (incl_surv == MODEL_TYPE_CODES[["joint"]] & beta_bool == 0){new_likelihood <- new_likelihood - SurvLike(beta_vec,surv_covar_risk_vec,surv_coef,survival_context)}
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood
#check to make sure all values are the same, simple sanity check
# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)
iter_count <- 1
stop_crit <- BASE_STOP_CRIT
if (!real_data){stop_crit <- stop_crit * SIM_STOP_CRIT_MULTIPLIER}
# if(mix_num > 8){stop_crit <- stop_crit * 10}
# if(mix_num > 12){stop_crit <- stop_crit * 10}
# if(mix_num > 15){stop_crit <- stop_crit * 5}

while((abs(like_diff) > stop_crit | iter_count < MIN_EM_ITERATIONS) & !run_only_surv){
  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  ##### MC Param  #####
  
  #### Mixing Proportion  #####
  re_prob <- CalcProbRE(alpha,pi_l)
  survival_context <- update_survival_context_re_prob(survival_context,
                                                      re_prob,fit_mix_num)
  
  ##### Survival ####
  #need model to fit a bit first otherwise may run into some instability
  if(beta_bool){

    nu_mat  <- CalcNu(nu_mat,re_prob,nu_covar_mat)
    pi_l <- CalcPi(nu_mat,nu_covar_mat)
    re_prob <- CalcProbRE(alpha,pi_l)
    survival_context <- update_survival_context_re_prob(survival_context,
                                                        re_prob,fit_mix_num)
    
    if (incl_surv == MODEL_TYPE_CODES[["joint"]]){
      #calculates survival coef for JM 
      beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef,fit_mix_num)
      beta_surv_coef_se <- CalcBeta(beta_surv_coef,combined_covar_mat,
                                    surv_covar_risk_vec,incl_surv,
                                    survival_context,surv_coef_len,fit_mix_num)
      beta_surv_coef_temp_list <- OutofBetaSurvCoef(beta_surv_coef_se[[1]],
                                                    surv_coef_len,fit_mix_num)
      beta_vec <- beta_surv_coef_temp_list[[1]]
      surv_coef <- beta_surv_coef_temp_list[[2]]
      beta_se <- beta_surv_coef_se[[2]]

      surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)

      bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,survival_context$re_prob,
                            surv_covar_risk_vec,survival_context$surv_event,
                            survival_context$surv_time,survival_context$surv_covar)
      bline_vec <- bhaz_vec[[1]]
      cbline_vec <- bhaz_vec[[2]]
    }
      
  }
  
  #### Weights  #####
  #calculates wake/sleep probabilities, needed for emission dist estimation
  weights_array_list <- CondMarginalize(alpha,beta,pi_l)
  weights_array_wake <- exp(weights_array_list[[1]])
  weights_array_sleep <- exp(weights_array_list[[2]])

  
  ##### Tobit bivariate normal emission update #####
  if(incl_light){

    emit_light[1,1,,] <- UpdateNorm(CalcLightMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
    emit_light[2,1,,] <- UpdateNorm(CalcLightMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

    emit_light[1,2,,] <- UpdateNorm(CalcLightSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
    emit_light[2,2,,] <- UpdateNorm(CalcLightSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
  }

  if (incl_act){
    emit_act[1,2,,] <- UpdateNorm(CalcActSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
    emit_act[2,2,,] <- UpdateNorm(CalcActSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

    emit_act[1,1,,] <- UpdateNorm(CalcActMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
    emit_act[2,1,,] <- UpdateNorm(CalcActMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
  }

  if (incl_act & incl_light){
    corr_mat[,1,] <- UpdateNorm(CalcBivarCorr,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
    corr_mat[,2,] <- UpdateNorm(CalcBivarCorr,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

  }
  
  ###
  #this only relies on normal parameters so calculate it now
  lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
  init <- CalcInit(alpha,beta,pi_l,log_sweights_vec)
  
  #saves old transition values in case likelihood decrease
  params_tran_array_old <- params_tran_array
  #gradient and hessian for tran parameters
  tran_gradhess_list <- CalcTranCHelper(alpha,beta,act,light,params_tran_array,
                                        emit_act,emit_light,corr_mat,
                                        pi_l,lod_act,lod_light,lintegral_mat,vcovar_mat,
                                        lambda_act_mat, lambda_light_mat, tobit, check_tran,likelihood)
  
  params_tran_array <- LM(tran_gradhess_list[[1]],tran_gradhess_list[[2]],params_tran_array,check_tran,likelihood,pi_l)
  
  tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)

  alpha <- Forward(act = act,light = light,
                   init = init,tran_list = tran_list,
                   emit_act= emit_act,emit_light = emit_light,
                   lod_act = lod_act, lod_light = lod_light,
                   corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef,surv_covar_risk_vec = surv_covar_risk_vec,
                   event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                   lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec,
                   surv_covar = surv_covar, vcovar_mat = vcovar_mat,
                   lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat, tobit = tobit,incl_surv = incl_surv*beta_bool)
  
  new_likelihood <- CalcLikelihood(alpha,pi_l)
  
  #if JM but during cold start, dont wan to actually include survial in likelihood yet
  if (incl_surv == MODEL_TYPE_CODES[["joint"]] & beta_bool == 0){new_likelihood <- new_likelihood - SurvLike(beta_vec,surv_covar_risk_vec,surv_coef,survival_context)}
  
  like_diff <- new_likelihood - likelihood
  
  if (like_diff < 0){
    #all other parameters are either
    #1) closed form
    #2) we can quickly calculate likelihood difference
    #tran likelihood requires forward algorithm and thus much slower
    #thus any like decrease is from transition parameters
    #effectively just doesnt optimize tran in this step of EM
    print("Transition Likelihood Decrease")

    params_tran_array <- params_tran_array_old
    tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)

    alpha <- Forward(act = act,light = light,
                     init = init,tran_list = tran_list,
                     emit_act = emit_act,emit_light= emit_light,
                     lod_act = lod_act, lod_light = lod_light,
                     corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef, surv_covar_risk_vec = surv_covar_risk_vec,
                     event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                     lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec,
                     surv_covar = surv_covar, vcovar_mat = vcovar_mat,
                     lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat, tobit = tobit,incl_surv = incl_surv)

    new_likelihood <- CalcLikelihood(alpha,pi_l)
    if (incl_surv == MODEL_TYPE_CODES[["joint"]] & beta_bool == 0){new_likelihood <- new_likelihood - SurvLike(beta_vec,surv_covar_risk_vec,surv_coef,survival_context)}
    like_diff <- new_likelihood - likelihood
  }
  
  print(paste("RE num:",mix_num,"Like:",round(like_diff,6)))
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  
  # after a few iterations don't have any issues with survival estimation stability
  if (iter_count == 3){
    beta_bool <- T
    print("Starting to Est Survival and Age Mixing Effect")
  }
  
  iter_count <- iter_count + 1
  
  
  #saves info every few iterations just in case
  if (iter_count %% INTERIM_SAVE_EVERY == 0){
    
    tran_df <- ParamsArray2DF(params_tran_array)
    if (!real_data & true_mix_num == fit_mix_num){
      tran_df_true <- ParamsArray2DF(params_tran_array_true)
      tran_df_truth <- tran_df_true[,1]
      
      tran_df <- tran_df %>% mutate(truth = tran_df_truth)
      tran_df <- tran_df %>% mutate(resid = prob - truth)
    }
    
    
    true_params <- make_true_param_list(init = init_true,
                                        params_tran_array = params_tran_array_true,
                                        emit_act = emit_act_true,
                                        emit_light = emit_light_true,
                                        corr_mat = corr_mat_true,
                                        nu_mat = nu_mat_true,
                                        beta_vec = beta_vec_true,
                                        beta_age = beta_age_true,
                                        lambda_act_mat = lambda_act_mat_true,
                                        lambda_light_mat = lambda_light_mat_true)

    est_params <- make_est_param_list(init = init,
                                      params_tran_array = params_tran_array,
                                      emit_act = emit_act,
                                      emit_light = emit_light,
                                      corr_mat = corr_mat,
                                      nu_mat = nu_mat,
                                      beta_vec = beta_vec,
                                      surv_coef = surv_coef,
                                      tran_df = tran_df,
                                      re_prob = re_prob,
                                      new_likelihood = new_likelihood,
                                      lambda_act_mat = lambda_act_mat,
                                      lambda_light_mat = lambda_light_mat)
    validate_param_list(est_params,fit_mix_num,vcovar_num,"est_params")

    to_save <- make_saved_results(true_params = true_params,
                                  est_params = est_params,
                                  settings = settings,
                                  start_params = start_params)
    
    if(!leave_out){
      save(to_save,file = paste0("Inter",model_name))
    }
    
      
  }
  
  
  #reorders clusters from best to worst survival
  if ((abs(like_diff) < stop_crit*REORDER_STOP_CRIT_MULTIPLIER) & !relabel_reset & !bootstrap & !leave_out & real_data){
    relabel_reset <- TRUE
    relabel_bool <- 0
    print("Relabelling")
    print("Potential Soft Reset")
    #### Reorder #####
    #Reorder to avoid label switching
    #Cluster means go from small to large by activity
    
    if (incl_surv == MODEL_TYPE_CODES[["joint"]]){
      reord_inds <- order(beta_vec)
    } else if (incl_surv == MODEL_TYPE_CODES[["two_stage"]]){
      beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef,fit_mix_num)
      beta_surv_coef_se <- CalcBeta(beta_surv_coef,combined_covar_mat,
                                    surv_covar_risk_vec,incl_surv,
                                    survival_context,surv_coef_len,fit_mix_num)
      beta_surv_coef_temp_list <- OutofBetaSurvCoef(beta_surv_coef_se[[1]],
                                                    surv_coef_len,fit_mix_num)
      beta_vec_temp <- beta_surv_coef_temp_list[[1]]
      reord_inds <- order(beta_vec_temp)
    }
    
    # reord_inds <- c(0,rev(order(beta_vec[-1])))+1
    if (!all(reord_inds == c(1:mix_num)) & !leave_out){
      print("Swapping Labels")
      relabel_bool <- 1
      emit_act <- emit_act[,,reord_inds,]
      emit_light <- emit_light[,,reord_inds,]
      nu_mat <- nu_mat[,reord_inds]
      nu_mat <- nu_mat - nu_mat[,1]
      params_tran_array <- params_tran_array[reord_inds,,]
      corr_mat <- corr_mat[reord_inds,,]
      init <- init[reord_inds,]
      beta_vec <- beta_vec[reord_inds]
      beta_vec <- beta_vec-min(beta_vec)

      pi_l <- CalcPi(nu_mat,nu_covar_mat)
      re_prob <- re_prob[,reord_inds]
      survival_context <- update_survival_context_re_prob(survival_context,
                                                          re_prob,fit_mix_num)

      bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,survival_context$re_prob,
                            surv_covar_risk_vec,survival_context$surv_event,
                            survival_context$surv_time,survival_context$surv_covar)
      bline_vec <- bhaz_vec[[1]]
      cbline_vec <- bhaz_vec[[2]]

    }
    
    for (re_ind in 1:mix_num){
      for (week_ind in 1:2){
        if (emit_act[2,1,re_ind,week_ind] > emit_act[1,1,re_ind,week_ind]){
          relabel_bool <- 1
          print(paste("Swapping wake/sleep for week_ind",week_ind,"Mixture",re_ind))

          if (week_ind == 1){
            temp <- init[re_ind,1]
            #NO WEEKEND INIT
            #SWAPPING MAY SLIGHTLEY DECREASE LIKE?
            init[re_ind,1] <- init[re_ind,2]
            init[re_ind,2] <- temp
          }

          temp <- emit_act[1,,re_ind,week_ind]
          emit_act[1,,re_ind,week_ind] <- emit_act[2,,re_ind,week_ind]
          emit_act[2,,re_ind,week_ind] <- temp

          temp <- emit_light[1,,re_ind,week_ind]
          emit_light[1,,re_ind,week_ind] <- emit_light[2,,re_ind,week_ind]
          emit_light[2,,re_ind,week_ind] <- temp

          #ISSUE HERE
          temp <- params_tran_array[re_ind,1:3,week_ind]
          params_tran_array[re_ind,1:3,week_ind] <- params_tran_array[re_ind,4:6,week_ind]
          params_tran_array[re_ind,4:6,week_ind] <- temp
          tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)

          temp <- corr_mat[re_ind,1,week_ind]
          corr_mat[re_ind,1,week_ind] <- corr_mat[re_ind,2,week_ind]
          corr_mat[re_ind,2,week_ind] <- temp
        }
      }
    }
    
    if (relabel_bool){
      tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)
      lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
      
      alpha <- Forward(act = act,light = light,
                       init = init,tran_list = tran_list,
                       emit_act= emit_act,emit_light = emit_light,
                       lod_act = lod_act, lod_light = lod_light, 
                       corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef,surv_covar_risk_vec = surv_covar_risk_vec,
                       event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                       lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec,
                       surv_covar = surv_covar, vcovar_mat = vcovar_mat,
                       lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat, tobit = tobit,incl_surv = incl_surv*beta_bool)
      
      new_likelihood <- CalcLikelihood(alpha,pi_l)
      like_diff <- stop_crit * 1.1
    }
      
  }
  
  
  #Finally calls beta at very end of while loop
  beta <- Backward(act = act,light = light, tran_list = tran_list,
                   emit_act = emit_act,emit_light = emit_light,
                   lod_act = lod_act, lod_light =  lod_light, 
                   corr_mat = corr_mat,lintegral_mat = lintegral_mat,vcovar_mat = vcovar_mat,
                   lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit)
  
}

#if 2-stage model, calculate survival here
if (incl_surv != MODEL_TYPE_CODES[["joint"]]){
  survival_context <- update_survival_context_re_prob(survival_context,
                                                      re_prob,fit_mix_num)
  beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef,fit_mix_num)
  beta_surv_coef_se <- CalcBeta(beta_surv_coef,combined_covar_mat,
                                surv_covar_risk_vec,incl_surv,
                                survival_context,surv_coef_len,fit_mix_num)
  beta_surv_coef_temp_list <- OutofBetaSurvCoef(beta_surv_coef_se[[1]],
                                                surv_coef_len,fit_mix_num)
  beta_vec <- beta_surv_coef_temp_list[[1]]
  surv_coef <- beta_surv_coef_temp_list[[2]]
  beta_se <- beta_surv_coef_se[[2]]
  
  surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)
  
  bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,survival_context$re_prob,
                        surv_covar_risk_vec,survival_context$surv_event,
                        survival_context$surv_time,survival_context$surv_covar)
  bline_vec <- bhaz_vec[[1]]
  cbline_vec <- bhaz_vec[[2]]
}

  
    



#if LOCV we don't care about viterbi decoding
if(!leave_out){
  decoded_mat <- Viterbi(act,light,vcovar_mat)
  
  
  weights_array_wake_collapsed <- apply(weights_array_wake,c(1,2),sum)
  weights_array_sleep_collapsed <- apply(weights_array_sleep,c(1,2),sum)
  post_decode <- weights_array_wake_collapsed < .5
  
  # post_decode <- weights_array_wake
  # for (i in 1:mix_num){
  #   post_decode[,,i] <- weights_array_wake[,,i] < weights_array_sleep[,,i]
  # }
  
} else {
  decoded_mat <- matrix(NA,2,2)
  post_decode <- matrix(NA,2,2)
  
}

#transition parameters for later analysis
tran_df <- ParamsArray2DF(params_tran_array)
if (!real_data & true_mix_num == fit_mix_num){
  tran_df_true <- ParamsArray2DF(params_tran_array_true)
  tran_df_truth <- tran_df_true[,1]
  
  tran_df <- tran_df %>% mutate(truth = tran_df_truth)
  tran_df <- tran_df %>% mutate(resid = prob - truth)
}

#concatenates true, starting, and estimated parameters into named lists to save
true_params <- make_true_param_list(init = init_true,
                                    params_tran_array = params_tran_array_true,
                                    emit_act = emit_act_true,
                                    emit_light = emit_light_true,
                                    corr_mat = corr_mat_true,
                                    nu_mat = nu_mat_true,
                                    beta_vec = beta_vec_true,
                                    beta_age = beta_age_true,
                                    lambda_act_mat = lambda_act_mat_true,
                                    lambda_light_mat = lambda_light_mat_true)

est_params <- make_est_param_list(init = init,
                                  params_tran_array = params_tran_array,
                                  emit_act = emit_act,
                                  emit_light = emit_light,
                                  corr_mat = corr_mat,
                                  nu_mat = nu_mat,
                                  beta_vec = beta_vec,
                                  surv_coef = surv_coef,
                                  tran_df = tran_df,
                                  re_prob = re_prob,
                                  new_likelihood = new_likelihood,
                                  decoded_mat = decoded_mat,
                                  lambda_act_mat = lambda_act_mat,
                                  lambda_light_mat = lambda_light_mat,
                                  bline_vec = bline_vec,
                                  cbline_vec = cbline_vec,
                                  beta_se = beta_se,
                                  post_decode = post_decode)
validate_param_list(est_params,fit_mix_num,vcovar_num,"est_params")

#if doing leave 100 out cross validation
#predict cluster assignment using varying levels of information
if (leave_out){
  new_act <- act_old[,leave_out_inds]
  new_light <- light_old[,leave_out_inds]
  new_vcovar_mat <- vcovar_mat_old[,leave_out_inds]
  len <- dim(new_act)[1]
  num_of_people <- dim(new_act)[2]
  new_surv_covar <- SubsetSurvCovar(surv_covar_old,leave_out_inds)
  new_pi_l <- CalcPi(nu_mat,nu_covar_mat_old[leave_out_inds,])
  
  surv_covar_risk_vec_new <- SurvCovarRiskVec(new_surv_covar,surv_coef)

  surv_event_new <- surv_event_old[leave_out_inds]
  surv_time_new <- surv_time_old[leave_out_inds]
  
  log_sweights_vec_new <- log_sweights_vec_old[leave_out_inds]


  
  empty_list <- vector(mode = "list", length = 6)
  for (i in 1:6){
    empty_list[[i]] <- list()
  }
  
  empty_mat_list <- vector(mode = "list", length = 6)
  for (i in 1:6){
    empty_mat_list[[i]] <- matrix(0,mix_num,mix_num)
  }
  
  empty_vec_list <- vector(mode = "list", length = 6)
  
  empty_mat_sublist <- vector(mode = "list", length = 3)
  for (i in 1:3){
    empty_mat_sublist[[i]] <- matrix(0,2,2)
  }
  
  
  conf_mat_list <- empty_mat_list
  cindex_new_list <- empty_vec_list
  ibs_new_list <- empty_vec_list
  ibs2_new_list <- empty_vec_list
  senspec_list <- empty_list
  senspec_mix_list <- empty_list

  
  
  
  #One is only cycle
  #Two is no light
  #Three is no activity
  #Four is no tran
  #Five is only act (no light\tran)
  #Six is standard
  
  for (leave_out_type in 1:6){
    new_act_working <- new_act
    new_light_working <- new_light
    
    if (leave_out_type == 2 | leave_out_type == 5){
      new_light_working <- matrix(NA,nrow = dim(new_act)[1],ncol = dim(new_act)[2])
    } 
    if (leave_out_type == 3){
      new_act_working <- matrix(NA,nrow = dim(new_act)[1],ncol = dim(new_act)[2])
    } 
    
    if (leave_out_type == 4 | leave_out_type == 5){
      tran_list <- GenTranList(array(0,dim = dim(params_tran_array)),c(1:day_length),mix_num,vcovar_num)
    } else {
      tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)
    }
    
    
    
    alpha <- Forward(act = new_act_working,light = new_light_working,
                     init = init,tran_list = tran_list,
                     emit_act = emit_act,emit_light = emit_light,
                     lod_act = lod_act, lod_light = lod_light,
                     corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef, surv_covar_risk_vec = surv_covar_risk_vec_new,
                     event_vec = numeric(100), bline_vec = numeric(100), cbline_vec = numeric(100),
                      lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec_new,
                      surv_covar = new_surv_covar, vcovar_mat = new_vcovar_mat,
                      lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,
                      tobit = T,incl_surv = MODEL_TYPE_CODES[["two_stage"]])
    
    beta <- Backward(act = new_act_working,light = new_light_working, tran_list = tran_list,
                     emit_act = emit_act,emit_light = emit_light,
                     lod_act = lod_act, lod_light =  lod_light, 
                     corr_mat = corr_mat,lintegral_mat = lintegral_mat,vcovar_mat = vcovar_mat,
                     lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit)
    
    
    weights_array_list <- CondMarginalize(alpha,beta,new_pi_l)
    weights_array_wake <- exp(weights_array_list[[1]])
    weights_array_wake_collapsed <- apply(weights_array_wake,c(1,2),sum)
    post_decode_collapsed <- weights_array_wake_collapsed < .5
    
    if (leave_out_type == 1){
      alpha <- ForwardAlt(post_decode_collapsed,init,tran_list,new_vcovar_mat)
    } 
    
    post_decode_collapsed_true_vec <- as.vector(post_decode_collapsed_true)
    post_decode_collapsed_vec <- as.vector(post_decode_collapsed)
    
    re_prob_new <- CalcProbRE(alpha,new_pi_l)
    mix_assignment_pred <- apply(re_prob_new,1,which.max)
    
    senspec_list[[leave_out_type]] <- empty_mat_sublist
    senspec_mix_list[[leave_out_type]] <- vector(mode = "list", length = 3)
    
    if (leave_out_type != 3){
      
      ind_med <- apply(new_act_working,2,median,na.rm = T)
      
      #1 - below
      #2 - above
      #3 - total
      for (lohi_med in 1:3){
        if (lohi_med == 1){
          valid_inds <- t(t(new_act_working) < ind_med)
          pdc_colnames <-  c("Below-Med Pred Wake","Pred Sleep")
        } else if (lohi_med == 2) {
          valid_inds <- t(t(new_act_working) > ind_med)
          pdc_colnames <-  c("Above-Med Pred Wake","Pred Sleep")
        } else {
          valid_inds <- new_act_working > -Inf
          pdc_colnames <-  c("Pred Wake","Pred Sleep")
        }
        valid_inds_vec <- as.vector(valid_inds)
        
        valid_inds_vec[is.na(valid_inds_vec)] <- F
        
        coda <- c(T,F)
        
        pdc_tab <- table(c(post_decode_collapsed_vec[valid_inds_vec],coda),c(post_decode_collapsed_true_vec[valid_inds_vec],coda))
        diag(pdc_tab) <- diag(pdc_tab) - 1
        rownames(pdc_tab) <- pdc_colnames
        colnames(pdc_tab) <- c("True Wake","True Sleep")
        senspec_list[[leave_out_type]][[lohi_med]] <- pdc_tab
        
        senspec_mix_list[[leave_out_type]][[lohi_med]] <- vector(mode = "list", length = mix_num)
        for (curr_class in 1:mix_num){
          # senspec_mix_list[[leave_out_type]][[curr_class]] <- empty_mat_sublist
          valid_inds_vec_mix <- as.vector(valid_inds[,mix_assignment_pred == curr_class])
          valid_inds_vec_mix[is.na(valid_inds_vec_mix)] <- F
          post_decode_collapsed_true_vec_mix <- as.vector(post_decode_collapsed_true[,mix_assignment_pred == curr_class])
          post_decode_collapsed_vec_mix <- as.vector(post_decode_collapsed[,mix_assignment_pred == curr_class])
          pdc_tab_mix <- table(c(post_decode_collapsed_vec_mix[valid_inds_vec_mix],coda),c(post_decode_collapsed_true_vec_mix[valid_inds_vec_mix],coda))
          diag(pdc_tab_mix) <- diag(pdc_tab_mix) - 1
          senspec_mix_list[[leave_out_type]][[lohi_med]][[curr_class]] <- pdc_tab_mix
        }
         
        

    }
      
      
    }
    
    # senspec_list[[leave_out_type]] <- sens_spec_df
   
    
    
    
    mix_assignment_pred <- c(mix_assignment_pred,c(1:mix_num))
    mix_assignment_true_ind <- c(mix_assignment_true[leave_out_inds],c(1:mix_num))
    conf_mat_ind <- table(mix_assignment_pred,mix_assignment_true_ind)
    diag(conf_mat_ind) <- diag(conf_mat_ind) - 1
    
    cindex <- CalcCindex(surv_time_new,surv_event_new,beta_vec,surv_coef,re_prob_new,new_surv_covar,surv_covar_risk_vec_new)
    ibs <- CalcIBS(surv_time_new,surv_event_new,cbline_vec,beta_vec,surv_coef,new_surv_covar,re_prob_new,incl_surv,mix_assignment_pred,surv_covar_risk_vec_new)
    ibs2 <- CalcIBS2(surv_time_new,surv_event_new,cbline_vec,beta_vec,re_prob_new,surv_covar_risk_vec_new)
  
    
    conf_mat_list[[leave_out_type]] <- conf_mat_ind
    cindex_new_list[[leave_out_type]] <-cindex
    ibs_new_list[[leave_out_type]] <- ibs
    ibs2_new_list[[leave_out_type]] <- ibs2

  }
    
  
  
  leave_out_to_save <- make_leave_out_results(leave_out_inds = leave_out_inds,
                                              conf_mat_list = conf_mat_list,
                                              cindex_new_list = cindex_new_list,
                                              ibs_new_list = ibs_new_list,
                                              senspec_list = senspec_list,
                                              ibs2_new_list = ibs2_new_list,
                                              senspec_mix_list = senspec_mix_list)
} else {
  leave_out_to_save <- list()
}

#if not using simulated data or want to save space
if (real_data | save_space){
  simulated_hmm <- list()
}


#mixture predections
mix_assignment <- apply(re_prob,1,which.max)
if (!real_data){
  true_class <- factor(c(as.vector(mixture_mat),seq_len(true_mix_num)), levels = seq_len(true_mix_num))
  fitted_class <- factor(c(mix_assignment,seq_len(fit_mix_num)), levels = seq_len(fit_mix_num))
  tab <- table(true_class,fitted_class)
  diag_ind <- seq_len(min(nrow(tab),ncol(tab)))
  tab[cbind(diag_ind,diag_ind)] <- tab[cbind(diag_ind,diag_ind)] - 1
} else {
  fitted_class <- factor(c(mix_assignment,seq_len(fit_mix_num)), levels = seq_len(fit_mix_num))
  tab <- table(fitted_class,fitted_class)
  diag(tab) <- diag(tab) - 1
}
  
#removes some data from saving
if (save_space){
  est_params$decoded_mat <- 0
  est_params$bline_vec <- 0
  est_params$cbline_vec <- 0
  est_params$post_decode <- 0
}

#diagnostics
ibs2 <- CalcIBS2(surv_time,surv_event,cbline_vec,beta_vec,re_prob,surv_covar_risk_vec)
ibs <- CalcIBS(surv_time,surv_event,cbline_vec,beta_vec,surv_coef,surv_covar,re_prob,incl_surv,mix_assignment_pred,surv_covar_risk_vec)
cindex <- CalcCindex(surv_time,surv_event,beta_vec,surv_coef,re_prob,surv_covar,surv_covar_risk_vec)
diagnostics <- make_diagnostics_list(cindex = cindex,
                                     ibs = ibs,
                                     confusion_table = tab,
                                     ibs2 = ibs2)
#save everything
bic <- CalcBIC(new_likelihood,mix_num,act,light)
to_save <- make_saved_results(true_params = true_params,
                              est_params = est_params,
                              bic = bic,
                              leave_out = leave_out_to_save,
                              simulated_hmm = simulated_hmm,
                              diagnostics = diagnostics,
                              settings = settings,
                              start_params = start_params)
setwd("/gpfs/gsfs12/users/aronjr/JM/Routputs")
#"~/JM/Routputs"
# model_name <- paste0("ReRun",model_name)
save(to_save,file = model_name)

