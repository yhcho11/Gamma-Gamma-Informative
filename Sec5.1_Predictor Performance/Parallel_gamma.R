# 2022/02/22 : argument "xij" added for all functions!! xij should be fixed!

Simul.Parallel <- function(iter_num,D, alpha, mean.u, var.u, Ni, beta0, beta1, variable_names, L, xij,a,b, phi, pop.model, wt.model){
  
  # cl <- makeCluster(round(detectCores()-1));
  cl <- makeCluster(35)
  
  clusterEvalQ(cl, { 
    
    all(sapply(c("lme4","sampling","data.table","dplyr","doParallel","fastGHQuad","merDeriv","MASS","numDeriv","Matrix","tidyr","reldist"), require, character.only = TRUE))
    source("fns_gammamodel.R")
    source("Parallel_gamma.R")
    
  })
  
  registerDoParallel(cl);
  
  res_list <- foreach( l = (1:iter_num), .errorhandling = "remove") %dopar% { 
    return( sub.fn.genresult(D, alpha, mean.u, var.u, Ni, beta0, beta1,variable_names, L, xij,a,b, phi, pop.model, wt.model) ) };
  
  stopCluster(cl);
  
  return(res_list)
  
}



sub.fn.genresult = function(D, alpha, mean.u, var.u, Ni, beta0, beta1, variable_names, L, xij,a,b, phi, pop.model, wt.model){
  
  h1 = function(y) {mean(y)}
  # h2 = function(y) {mean(exp(y))}
  h2 = function(y) {quantile(y, probs = 0.25)} #Q1
  h3 = function(y) {quantile(y, probs = 0.75)} #Q3
  # h6 = function(y) {
  #   ind <- ifelse(exp(y) < 155, 1, 0)   ## changed from 90 to 100 : 08162021
  #   val <- (1-exp(y)/155)*ind
  #   mean(val)
  # }
  h4 = function(y) {gini(y)} # changed 08192021   # h7 = function(y) {gini(y)}
  hlist <- list(h1 = h1, h2 = h2, h3 = h3, h4 = h4)
  nh = length(hlist)
  T.par = paste0("h",1:nh)
  
  
  dat.sets = gen.dat(D, alpha, mean.u, var.u, Ni, beta0, beta1, xij,a,b,phi, pop.model, wt.model)
  for(i in 1:length(dat.sets)){assign(names(dat.sets)[i], dat.sets[[i]])}
  
  #Proposed Method
  ini.values = initials.propsed(datsamp)
  est.par = optim(ini.values, fn = loglik.pars, datsamp = datsamp)$par
  est.par["alpha"] = abs(est.par["alpha"]); est.par["delta"] = abs(est.par["delta"]) ;
  
  
  #######################################################################################################################
  #### ADDED for Informative Sampling:
  #######################################################################################################################
  hat_coef_ab = wt_par_ab(datsamp, D)
  hat.a = hat_coef_ab[1]; hat.b = hat_coef_ab[2];
  kappa.i = hat_coef_ab[(1:D)+2]
  #######################################################################################################################
  TR <- sapply(1:nh, function(i){ 
    tapply(datpop$yij, datpop$area, hlist[[i]])}) %>%
    as.data.frame(.) %>% setNames(T.par) 
  
  EBs <- ebs.pred(est.par, datsamp, datpop, nis, Nis, D, L, hlist, nh)
  
  EB = EBs$estpssampa %>% .[,T.par]
  M1hat = EBs$M1hat %>% .[,T.par]
  #######################################################################################################################
  #### ADDED for Informative Sampling:
  #######################################################################################################################
  EBs.INFOR <- ebs.pred.INFOR(est.par, datsamp, datpop, nis, Nis, D, L, hlist, nh, hat.a, hat.b)
  EB_infor <- EBs.INFOR$estpssampa[,T.par]
  
  EBs.INFOR.fciyij <- ebs.pred.INFOR.fciyij(est.par, datsamp, datpop, nis, Nis, D, L, hlist, nh, hat.a, hat.b, kappa.i) ## changed to (a,b) to (hat.a, hat.b) 07262024
  EB_infor_fci <- EBs.INFOR.fciyij$estpssampa[,T.par]
  
  # EBs.INFOR.fp <- ebs.pred.INFOR.fp(est.par, datsamp, datpop, nis, Nis, D, L, hlist, nh, a,b)
  # EB_infor_fp <- EBs.INFOR.fp$estpssamp[,T.par]
  ##########################################################################################################    
  
  ##########################################################################################################  
  Bs <- ebs.pred.INFOR(true.par, datsamp, datpop, nis, Nis, D, L, hlist, nh, a, b) # changed 07082022
    # ebs.pred(true.par, datsamp, datpop, nis, Nis, D, L, hlist, nh)
  B = Bs$estpssampa %>% .[,T.par]
  
  ## Closed mean predictor added:
  # B.clsd_h1 <- ebs_clsd.pred(true.par, datsamp, datpop, nis, Nis, D, hlist, nh,TR)%>% dplyr::select(h1) %>% unlist(.)
  # EB.clsd_h1 <- ebs_clsd.pred(est.par, datsamp, datpop, nis, Nis, D,  hlist, nh,TR,hat.b)%>% dplyr::select(h1) %>% unlist(.)
  # TMSE_B.clsd_h1 <- (B.clsd_h1-TR[,"h1"])^2 %>% unlist(.)
  # TMSE_EB.clsd_h1 <- (EB.clsd_h1-TR[,"h1"])^2 %>% unlist(.)
 ##########################################################################################################  
  
 ###############################################################################################
  
  ## Hobza Model :
  
  glme.ini.log = glmer(yij ~ x+ (1|area) , family = Gamma(link="log"),  data = datsamp) # added 20220206 ## use loglink
  # glme.inv = glmer(yij ~ x + (1|area), family = Gamma(link="inverse"),  data = datsamp,
  #                  control = glmerControl(tolPwrss=1e-30),  start = list(c(fixef(glme.inv),vcov(glme.inv)[1,1])))
  
  HZ.est.pars = HZ.est.par(glme.ini.log, datsamp)
  
  EB_HZ = HZ.EBP(HZ.est.pars, datsamp, datpop, nis, Nis, D, L1=L, L2=L,hlist, nh)
  
  EB_HZplug = HZ.PLUG(glme.ini.log, HZ.est.pars, datsamp, datpop, nis, Nis, D, hlist, nh)
  
  EB_HZmarginal = HZ.Marginal(glme.ini.log, HZ.est.pars, datsamp, datpop, nis, Nis, D,hlist, nh,L)
  
  EB_dir = sapply(1:nh, function(i){ tapply(datsamp$yij, datsamp$area, hlist[[i]])}) %>% as.data.frame(.) %>% setNames(paste0("h",1:nh))
  
  
  # Caluculate True MSE :
  
  # T.methods = c("B","EB","EB_HZ","EB_HZplug","EB_HZmarginal", "EB_dir", "EB_infor","EB_infor_fci","EB_infor_fp") # "B" added 05162022
  T.methods = c("B","EB","EB_HZ","EB_HZplug","EB_HZmarginal", "EB_dir", "EB_infor","EB_infor_fci") # "B" added 05162022
  # T.methods = c("B","EB","EB_HZ","EB_HZplug","EB_HZmarginal", "EB_dir", "EB_infor","EB_infor_fp") # "B" added 05162022
  for(i in T.methods){                                                 #  TMSE_EB_HZ = Tmse(EB_HZ, TR, nh)
    Temp.text = paste0("TMSE_",i,"= Tmse(", i,", TR, nh)")
    eval(parse(text = Temp.text))
  }
  
  # Calculate True Bias:
  for(i in T.methods){                                                 #  TMSE_EB_HZ = Tmse(EB_HZ, TR, nh)
    Temp.text = paste0("TBIAS_",i,"= Tbias(", i,", TR, nh)")
    eval(parse(text = Temp.text))
  }  
  
  T.variablenames = c("TR",T.methods,  paste0("TMSE_", T.methods),  paste0("TBIAS_", T.methods))
  T.par = paste0("h",1:nh)
  
  # updated h8 gini(abs(y))
  # variable_names = expand.grid(T.variablenames,T.par) %>% apply(.,1,function(.){paste0(.,collapse = "_")})
  # variable_names = c(variable_names)
  
  for(i in T.variablenames){
    for(j in T.par){
      Text = paste0(i,"_",j," = ",i,"[,'",j,"']")
      eval(parse(text = Text))
    }
  }
  
    
  res_df = lapply(variable_names, function(.){eval(parse(text = .))}) %>% Reduce("rbind",.)
  colnames(res_df) = as.character(1:D)
  
  return(res_df)
  
}
