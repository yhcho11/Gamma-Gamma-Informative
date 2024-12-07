# xij added to all functions 20220222

Simul.Parallel <- function(iter_num,D, alpha, mean.u, var.u, Ni, beta0, beta1, variable_names, L, B, xij){
  
  # cl <- makeCluster(round(detectCores()-1));
  cl <- makeCluster(35)
  clusterEvalQ(cl, { 
    
    all(sapply(c("lme4","sampling","data.table","dplyr","doParallel","fastGHQuad","merDeriv","MASS","numDeriv","Matrix","tidyr","reldist"), require, character.only = TRUE))
    source("fns_gammamodel.R")
    source("Parallel_gamma.R")
    source("fns_MSE_gamma.R")
    
  })
  
  registerDoParallel(cl);
  
  res_list <- foreach( l = (1:iter_num), .errorhandling = "remove") %dopar% { 
    return( sub.fn.genresult(D, alpha, mean.u, var.u, Ni, beta0, beta1,variable_names, L, B, xij) ) };
  
  stopCluster(cl);
  
  return(res_list)
  
}



sub.fn.genresult = function(D, alpha, mean.u, var.u, Ni, beta0, beta1, variable_names, L, B, xij){
  
  h1 = function(y) {mean(y)}
  h2 = function(y) {quantile(y, probs = 0.25)} #Q1
  h3 = function(y) {quantile(y, probs = 0.75)} #Q3
  h4 = function(y) {gini(y)} # changed 08192021   # h7 = function(y) {gini(y)}
  hlist <- list(h1 = h1, h2 = h2, h3 = h3, h4 = h4)
  nh = length(hlist)
  T.par = paste0("h",1:nh)
  
  
  dat.sets = gen.dat(D, alpha, mean.u, var.u, Ni, beta0, beta1, xij)
  for(i in 1:length(dat.sets)){assign(names(dat.sets)[i], dat.sets[[i]])}
  
  #Proposed Method
  ini.values = initials.propsed(datsamp)
  est.par = optim(ini.values, fn = loglik.pars, datsamp = datsamp)$par
  est.par["alpha"] = abs(est.par["alpha"]); est.par["delta"] = abs(est.par["delta"]) ;
  
  # L=200; B=100;
  
  TR <- sapply(1:nh, function(i){ 
    tapply(datpop$yij, datpop$area, hlist[[i]])}) %>%
    as.data.frame(.) %>% setNames(T.par) # D * H4
  
  Bs <- ebs.pred.CI(true.par, datsamp, datpop, nis, Nis, D, L, hlist, nh, TR) ## best predictor is added03012022
  EBs <- ebs.pred.CI(est.par, datsamp, datpop, nis, Nis, D, L, hlist, nh, TR)
  
  
  BP = Bs$estpssampa %>% .[,T.par]  ## best predictor is added03012022
  EB = EBs$estpssampa %>% .[,T.par]
  M1hat = EBs$M1hat %>% .[,T.par]
  
  Naive.90 = EBs$Naive.90
  Naive.95 = EBs$Naive.95
  Naive.99 = EBs$Naive.99
  
  
  ## MSE Estimation:
  ## 1. MSE Bias Correction Estimation
  MSE_BC_term =  MSE_BC(est.par, datsamp, datpop, nis, Nis, D, L, B, hlist,EB)
  M1hat.bt  = MSE_BC_term$M1hat_bt
  M2hat.bt = MSE_BC_term$M2hat_bt
  # ## 1. MSE Bias Correction Estimation
  # M1hat200.bt  = MSE_BC_term$M1hat200_bt ## added 03132022 muted 08122024
  # M2hat200.bt = MSE_BC_term$M2hat200_bt  ## added 03132022 muted 08122024
  
  # theta_rep_boot = MSE_BC_term$theta_rep_boot  ## editted 11282021 
  
  rm(MSE_BC_term) ## 07012021 editted
  gc()
  
  ## True MSE
  
  MSE.True = Tmse(EB, TR, nh)
  M1 = Tmse(BP,TR,nh)
  M2 = Tmse(EB,BP,nh)
  MSE.TRcond = M1 + M2## best predictor is added03012022
  ## MSE without Correction
  MSE.noBC = M1hat+M2hat.bt
  
  MSE.add = (2*M1hat-M1hat.bt)+M2hat.bt
  MSE.mult = M1hat*(M1hat/M1hat.bt)+M2hat.bt
  null.mat = matrix(NA, nrow  = length(nis), ncol = nh)
  null.mat[M1hat<M1hat.bt] = M1hat[M1hat<M1hat.bt]*exp(-(M1hat.bt[M1hat<M1hat.bt]-M1hat[M1hat<M1hat.bt])/M1hat.bt[M1hat<M1hat.bt])
  null.mat[!M1hat<M1hat.bt] = 2*M1hat[!M1hat<M1hat.bt]-M1hat.bt[!M1hat<M1hat.bt]
  MSE.hm = null.mat + M2hat.bt
  rm(null.mat);gc();
  
  null.mat2 = matrix(NA, nrow  = length(nis), ncol = nh)
  null.mat2[M1hat<M1hat.bt] =  M1hat[M1hat<M1hat.bt]*(M1hat[M1hat<M1hat.bt]/M1hat.bt[M1hat<M1hat.bt])
  # null.mat2[M1hat<M1hat.bt] = M1hat[M1hat<M1hat.bt]*exp(-(M1hat.bt[M1hat<M1hat.bt]-M1hat[M1hat<M1hat.bt])/M1hat.bt[M1hat<M1hat.bt])
  null.mat2[!M1hat<M1hat.bt] = 2*M1hat[!M1hat<M1hat.bt]-M1hat.bt[!M1hat<M1hat.bt]
  MSE.comp = null.mat2 + M2hat.bt
  rm(null.mat2);gc();
  
  #### muted B = 200 08122024
  ## MSE Estimation with B=200
  
  # ## 1. MSE Bias Correction Estimation
  # M1hat200.bt  = MSE_BC_term200$M1hat_bt
  # M2hat200.bt = MSE_BC_term200$M2hat_bt
  # # theta_rep_boot = MSE_BC_term$theta_rep_boot  ## editted 11282021
  
  # ## MSE without Correction
  # MSE.noBC200 = M1hat+M2hat200.bt
  # 
  # MSE.add200 = (2*M1hat-M1hat200.bt)+M2hat200.bt
  # MSE.mult200 = M1hat*(M1hat/M1hat200.bt)+M2hat200.bt
  # null.mat = matrix(NA, nrow  = length(nis), ncol = nh)
  # null.mat[M1hat<M1hat200.bt] = M1hat[M1hat<M1hat200.bt]*exp(-(M1hat200.bt[M1hat<M1hat200.bt]-M1hat[M1hat<M1hat200.bt])/M1hat200.bt[M1hat<M1hat200.bt])
  # null.mat[!M1hat<M1hat200.bt] = 2*M1hat[!M1hat<M1hat200.bt]-M1hat200.bt[!M1hat<M1hat200.bt]
  # MSE.hm200 = null.mat + M2hat200.bt
  # rm(null.mat);gc();
  # 
  # null.mat2 = matrix(NA, nrow  = length(nis), ncol = nh)
  # null.mat2[M1hat<M1hat200.bt] =  M1hat[M1hat<M1hat200.bt]*(M1hat[M1hat<M1hat200.bt]/M1hat200.bt[M1hat<M1hat200.bt])
  # # null.mat2[M1hat<M1hat.bt] = M1hat[M1hat<M1hat.bt]*exp(-(M1hat.bt[M1hat<M1hat.bt]-M1hat[M1hat<M1hat.bt])/M1hat.bt[M1hat<M1hat.bt])
  # null.mat2[!M1hat<M1hat200.bt] = 2*M1hat[!M1hat<M1hat200.bt]-M1hat200.bt[!M1hat<M1hat200.bt]
  # MSE.comp200 = null.mat2 + M2hat200.bt
  # rm(null.mat2);gc();
  
  
  ## 2. MSE Fully Parameteric :
  
  MSE.Full <- MSE_Fully(est.par, datsamp,  datpop, nis, Nis, D, samparea, L, hlist, B)
  
  
  MSE.single = MSE.Full$MSE_single
  MSE.double = MSE.Full$MSE_double
  MSE.bc_double = MSE.Full$MSE_bc_double
  
  
  # T.mthd = c("add","mult","hm","comp","noBC","add200", "mult200","hm200","comp200","noBC200","single") # B200  version added 0311/2022
  T.mthd = c("add","mult","hm","comp","noBC","single") # B200  version added 0311/2022
  ncp = c(0.9,0.95,0.99)
  
  for(i in T.mthd){
    for (j in ncp){
      Text1 = paste0("Norm.",i,".",j*100," = CI_norm_Ind(nh, TR, EB, MSE.",i,",",j,")")
      Text2 = paste0("colnames(Norm.",i,".",j*100,") = T.par")
      eval(parse(text = Text1))
      eval(parse(text = Text2))
    }
  }
  
  Ld.M1_h1 =  LeadTerm(true.par, datpop, datsamp)
  Ld.M1hat_h1 = LeadTerm(est.par, datpop, datsamp)
  yetaboot = LeadTerm_BT(est.par, datsamp, datpop, nis, Nis, D, B)
  
  # Ld.M1hat.bt_h1 = apply(yetaboot[,1:(ncol(yetaboot)/2)],1,mean) # muted 08122024
  # Ld.M1hat200.bt_h1 = apply(yetaboot,1,mean) # muted 08122024
  Ld.M1hat.bt_h1 = apply(yetaboot,1,mean)
  
  ### leadingterm bias for bias test
  
  lead.bias = M1hat - M1
  est.lead.bias = M1hat.bt - M1hat
  
  ########################################
  # Save to variabl_names:
  
  EBPs = c("BP","EB","TR")  ## BP : best predictor is added03012022
  # M1s = c("M1","M1hat","M1hat.bt","M1hat200.bt","M2","M2hat.bt","M2hat200.bt","lead.bias","est.lead.bias")  muted 08122024
  M1s = c("M1","M1hat","M1hat.bt","M2","M2hat.bt","lead.bias","est.lead.bias") 
  MSEs = expand.grid("MSE",c("True","TRcond","noBC","add", "mult","hm","comp","single","double","bc_double")) %>%  # B=200  added 03112022
    apply(.,1,function(.){paste0(.,collapse = ".")}) ## TR.cond : best predictor is added03012022
  CIs =  expand.grid("Norm",c("add","mult","hm","comp","noBC","single"),c(90,95,99)) %>% apply(.,1,function(.){paste0(.,collapse = ".")})
  naiveCIs = paste("Naive",c(90,95,99), sep= ".")
  T.variablenames = c(EBPs,M1s,MSEs,CIs, naiveCIs)  # 03012022
  T.par = paste0("h",1:4) # updated h8 gini(abs(y))
  variable_names = expand.grid(T.variablenames,T.par) %>% apply(.,1,function(.){paste0(.,collapse = "_")})

  
  for(i in T.variablenames){
    for(j in T.par){
      Text = paste0(i,"_",j," = ",i,"[,'",j,"']")
      eval(parse(text = Text))
    }
  }
  
  # variable_names2 = c(variable_names,"Ld.M1_h1","Ld.M1hat_h1","Ld.M1hat.bt_h1","Ld.M1hat200.bt_h1")   ## Added 03112022
  # res_df = lapply(variable_names2, function(.){eval(parse(text = .))}) %>% Reduce("rbind",.)
  # colnames(res_df) = as.character(1:D)

  
  ################ added  2022/03/26  Comparison between Methods for closed formula of h1 #################
  MSE.TRcond.CL_h1 = Ld.M1_h1 + M2_h1  # closed formula of h1
  MSE.noBC.CL_h1 = Ld.M1hat_h1 + M2hat.bt_h1 
  MSE.add.CL_h1 = (2*Ld.M1hat_h1-Ld.M1hat.bt_h1)+M2hat.bt_h1
  MSE.mult.CL_h1 = Ld.M1hat_h1*(Ld.M1hat_h1/Ld.M1hat.bt_h1) + M2hat.bt_h1
  null.mat = as.vector(matrix(NA, nrow  = 60, ncol = 1))
  null.mat[Ld.M1hat_h1<Ld.M1hat.bt_h1] = Ld.M1hat_h1[Ld.M1hat_h1<Ld.M1hat.bt_h1]*exp(-(Ld.M1hat.bt_h1[Ld.M1hat_h1<Ld.M1hat.bt_h1]-Ld.M1hat_h1[Ld.M1hat_h1<Ld.M1hat.bt_h1])/Ld.M1hat.bt_h1[Ld.M1hat_h1<Ld.M1hat.bt_h1])
  null.mat[!Ld.M1hat_h1<Ld.M1hat.bt_h1] = 2*Ld.M1hat_h1[!Ld.M1hat_h1<Ld.M1hat.bt_h1]-Ld.M1hat.bt_h1[!Ld.M1hat_h1<Ld.M1hat.bt_h1]
  MSE.hm.CL_h1 = null.mat + M2hat.bt_h1
  rm(null.mat);gc();

  null.mat2 = as.vector(matrix(NA, nrow  = 60, ncol = 1))
  null.mat2[Ld.M1hat_h1<Ld.M1hat.bt_h1] =  Ld.M1hat_h1[Ld.M1hat_h1<Ld.M1hat.bt_h1]*(Ld.M1hat_h1[Ld.M1hat_h1<Ld.M1hat.bt_h1]/Ld.M1hat.bt_h1[Ld.M1hat_h1<Ld.M1hat.bt_h1])
  # null.mat2[M1hat<M1hat.bt] = M1hat[M1hat<M1hat.bt]*exp(-(M1hat.bt[M1hat<M1hat.bt]-M1hat[M1hat<M1hat.bt])/M1hat.bt[M1hat<M1hat.bt])
  null.mat2[!Ld.M1hat_h1<Ld.M1hat.bt_h1] = 2*Ld.M1hat_h1[!Ld.M1hat_h1<Ld.M1hat.bt_h1]-Ld.M1hat.bt_h1[!Ld.M1hat_h1<Ld.M1hat.bt_h1]
  MSE.comp.CL_h1 = null.mat2 + M2hat.bt_h1
 
  variable_names2 = c(variable_names
                      , "Ld.M1_h1","Ld.M1hat_h1","Ld.M1hat.bt_h1"
                      , "MSE.TRcond.CL_h1", "MSE.noBC.CL_h1", "MSE.add.CL_h1", "MSE.mult.CL_h1", "MSE.hm.CL_h1", "MSE.comp.CL_h1")   ## Added 03112022
  res_df = lapply(variable_names2, function(.){eval(parse(text = .))}) %>% Reduce("rbind",.)
  colnames(res_df) = as.character(1:D)
  

  
  return(res_df)
  
}
