
CI_Ind = function(nh, popsampa, lcl_mat, ucl_mat){
  res_Ind = sapply(1:nh, function(i){Ind = (popsampa[,i]>=lcl_mat[,i]) & (popsampa[,i]<=ucl_mat[,i])})
  return(res_Ind)
}


CI_norm_Ind = function(nh, popsampa, EBP, MSE_est_mat, ncp=0.95){
  lcl_mat = EBP - qnorm(1-(1-ncp)/2) * sqrt(MSE_est_mat)
  ucl_mat = EBP + qnorm(1-(1-ncp)/2) * sqrt(MSE_est_mat)
  CI_Ind(nh, popsampa, lcl_mat, ucl_mat)
  
}

ebs.pred.CI = function(est.par, datsamp, datpop, nis, Nis, D, L, hlist, nh,TR){
  
  yNs = replicate(L, onecycle.ebs.pred(est.par, datsamp, datpop, nis, Nis, D))  
  T.par = paste0("h",1:nh)
  samparea = unique(datsamp$area)
  
  theta_rep <- sapply(1:nh, function(i) {                                     # changed 08192021 for gini coefficient  abs(y) and exp(y) in gini
    apply(yNs, 2 , function(y){tapply((y), datpop$area, hlist[[i]])})}) %>%  as_tibble() %>%
    setNames(paste0("h", 1:nh)) %>% 
    cbind(area = rep(samparea, times = L))
  
  estpssampa <- theta_rep %>% group_by(area) %>% summarise_all(mean) %>% as.data.frame(.) ## EBUP ####
  ## observed conditional variance => estimator of MSE leading term
  M1hat <- theta_rep %>% group_by(area) %>% summarise_all(var) %>% as.data.frame(.)
  
  ## naive CIs ####
  for(ncp in c(0.9,0.95,0.99)){
    lcl_naive <- theta_rep %>% group_by(area) %>% summarise_all(quantile, prob = (1-ncp)/2) %>% as.matrix(.)  %>% .[,T.par]
    ucl_naive <- theta_rep %>% group_by(area) %>% summarise_all(quantile, prob = 1-(1-ncp)/2) %>% as.matrix(.)  %>% .[,T.par]
    Text1 = paste0("Naive.",ncp*100," = CI_Ind(nh, TR, lcl_naive, ucl_naive)")
    Text2 = paste0("colnames(Naive.",ncp*100,") = T.par")
    eval(parse(text = Text1))
    eval(parse(text = Text2))
    rm(lcl_naive, ucl_naive)
    gc()
    
  }
  
  return(list (estpssampa = estpssampa
               , theta_rep = theta_rep
               , M1hat = M1hat
               , Naive.90 = Naive.90
               , Naive.95 = Naive.95
               , Naive.99 = Naive.99))
}

########################################################################################################################################################################
########################################################################################################################################################################

### Bias- Corrected MSE :

EBP_infor <- function(est.par, datsamp, datpop, nis, Nis, D, L){
  
  ## Genereate Bootsrap Parameter estmates:
  hat.alpha = est.par["alpha"];  hat.delta = est.par["delta"]; 
  hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"]
  
  # Boostrap sample:
  usD.b <- rgamma(D, shape = hat.delta, rate = hat.delta) %>% setNames(1:D)
  betaij.b <- exp(hat.beta0  + hat.beta1*datsamp$x)*usD.b[datsamp$area]
  yij.b <- rgamma(length(betaij.b), shape = hat.alpha, rate = betaij.b)
  
  # Bootstrap parameter estimates: 
  datsamp.b = datsamp; datsamp.b$yij = yij.b;
  ini.values.b = initials.propsed(datsamp.b)
  est.par.b = optim(ini.values.b, fn = loglik.pars, datsamp = datsamp.b)$par

  yNs = replicate(L, onecycle.ebs.pred(est.par.b, datsamp, datpop, nis, Nis, D))
    
  yNs[datpop$Iunits==1] <- datsamp$yij
  
  rm(datpop, datsamp, datsamp.b);gc(); ##updated 08192021
  
  return(t(yNs))
}


MSE_BC <- function(est.par, datsamp, datpop, nis, Nis, D, L=200,B, hlist,EB){ # When B=200 
  nh = length(hlist)
  T.par =  paste0("h",1:nh) 
  
  yetaboot <- sapply(1:B, function(i){EBP_infor(est.par, datsamp, datpop, nis, Nis, D, L)})
  colnames(yetaboot) <-1:B 
  
  yetaboot <- yetaboot %>% cbind.data.frame(
    rep = rep(1:L, times = nrow(datpop)),
    area = rep(datpop$area, each = L)
  )
###########################################################################  
  theta_rep_boot <- lapply(1:B, function(i){
      yetaboot %>% dplyr::select(i, area, rep) %>% 
      setNames(c("theta", "area", "rep")) %>%
      group_by(area, rep) %>%  ## order has been changed
      summarise(h1 = hlist[["h1"]](theta), h2 = hlist[["h2"]](theta), h3 = hlist[["h3"]](theta), h4 = hlist[["h4"]](theta), .groups = 'drop') %>%  # h8 added 08192021
      mutate(boot = i) %>% 
      return(.)
  }) %>% do.call("rbind",.)
  
  rm(yetaboot);gc();
  
  M1hat_boot <- theta_rep_boot %>% dplyr::select(-rep) %>% 
    group_by(boot, area)  %>% summarise_all(var)  %>% group_by(area) %>% summarise_at(vars(contains("h")), mean) %>% 
    dplyr::select(T.par)%>% as.matrix(.)  # return 120*6 matrix
  
  
  ## Boot M2 ####
  estpssampa_boot <- theta_rep_boot %>% dplyr::select(-rep) %>% group_by(area,boot) %>% 
    summarise_all(mean) %>% arrange(boot,area) %>% ungroup()
  # gc()
  
  # rm(theta_rep_boot) ## 07012021 editted ## muted 11282021 
  # gc()
  
  M2hat_boot <- sapply(1:nh, function(hi){
    colname <- paste0("h", hi)
    tapply((unlist(estpssampa_boot[, colname]) - (EB[, colname]))^2, estpssampa_boot$area, mean)
  })
  
  colnames(M2hat_boot) = paste0("h",1:nh)
  
  rm(estpssampa_boot);gc() # editted 08172021
  # B =200 muted 08122024
  # 
  # ################################B=200########################################### added part 03132022 Start:
  # 
  # theta_rep_boot_200 <- lapply(1:B, function(i){
  #   yetaboot %>% dplyr::select(i, area, rep) %>% 
  #     setNames(c("theta", "area", "rep")) %>%
  #     group_by(area, rep) %>%  ## order has been changed
  #     summarise(h1 = hlist[["h1"]](theta), h2 = hlist[["h2"]](theta), h3 = hlist[["h3"]](theta), h4 = hlist[["h4"]](theta), .groups = 'drop') %>%  # h8 added 08192021
  #     mutate(boot = i) %>% 
  #     return(.)
  # }) %>% do.call("rbind",.)
  # 
  # rm(yetaboot);gc();
  # 
  # M1hat200_boot <- theta_rep_boot_200 %>% dplyr::select(-rep) %>% 
  #   group_by(boot, area)  %>% summarise_all(var)  %>% group_by(area) %>% summarise_at(vars(contains("h")), mean) %>% 
  #   dplyr::select(T.par)%>% as.matrix(.)  # return 120*6 matrix
  # 
  # 
  # ## Boot M2 ####
  # estpssampa200_boot <- theta_rep_boot_200 %>% dplyr::select(-rep) %>% group_by(area,boot) %>% 
  #   summarise_all(mean) %>% arrange(boot,area) %>% ungroup()
  # # gc()
  # 
  # # rm(theta_rep_boot) ## 07012021 editted ## muted 11282021 
  # # gc()
  # 
  # M2hat200_boot <- sapply(1:nh, function(hi){
  #   colname <- paste0("h", hi)
  #   tapply((unlist(estpssampa200_boot[, colname]) - (EB[, colname]))^2, estpssampa200_boot$area, mean)
  # })
  # 
  # colnames(M2hat200_boot) = paste0("h",1:nh)
  # 
  # rm(estpssampa200_boot);gc() # editted 08172021
  # ################################B=200########################################### added part 03132022 end:
  
  return(list( M1hat_bt = M1hat_boot,
               M2hat_bt = M2hat_boot
               # M1hat200_bt = M1hat200_boot,
               # M2hat200_bt = M2hat200_boot
               # theta_rep_boot = theta_rep_boot ## editted 11282021 for unconditional CI
  ))
  
}


####################################################################################################################################

one.cycle.fully = function(est.par, datsamp,  datpop, nis, Nis, D, samparea, L, hlist, seed, nh, T.par){
  
  set.seed(seed);
  
  hat.alpha = est.par["alpha"]; hat.delta = est.par["delta"]; 
  hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"];
  
  # Generate yij 
  usD.b <- rgamma(D, shape = hat.delta, rate = hat.delta) %>% setNames(1:D)
  betaij.b <- exp(hat.beta0 + hat.beta1*datpop$x)*usD.b[datpop$area]
  yij.b <- rgamma(length(betaij.b), shape = hat.alpha, rate = betaij.b)  
  
  # level-one : Bootstrap ver. Parameters:
  datpop.b = datpop; datpop.b$yij = yij.b;
  TR.b = sapply(1:nh, function(i){ 
    tapply(datpop.b$yij, datpop.b$area, hlist[[i]])}) %>%
    as.data.frame(.) %>% setNames(T.par) 
  rm(datsamp); gc();
  
  # Bootstrap ver. predictors:
  datsamp.b = datpop.b %>% filter(Iunits==1)
  ini.values.b = initials.propsed(datsamp.b)
  est.par.b = optim(ini.values.b, fn = loglik.pars, datsamp = datsamp.b)$par
  est.par.b["alpha"] = abs(est.par.b["alpha"]); est.par.b["delta"] = abs(est.par.b["delta"]) ;
  
  EB.b = ebs.pred(est.par.b, datsamp.b, datpop, nis, Nis, D, L , hlist, nh)$estpssampa %>% .[,T.par]
  
  rm(datpop.b, datsamp.b); gc();
  
  MSE_bt = ((EB.b - TR.b)^2) 
  
  return(list(est.par.b = est.par.b
              , MSE_bt = MSE_bt))
  
}



EBP_fully <- function(est.par, datsamp,  datpop, nis, Nis, D, samparea, L, hlist, nh, T.par){
 
  
  seeds = sample(1:2^30, 2);
  seed1 = seeds[1]; seed2 = seeds[2];
  
  # level-one Bootstrap :
  lv_one = one.cycle.fully(est.par, datsamp,  datpop, nis, Nis, D, samparea, L, hlist, seed = seed1, nh, T.par)
  lv_one_est.par = lv_one$est.par.b
  lv_one_MSE = lv_one$MSE_bt
  
  # level-two Bootstrap :  
  lv_two = one.cycle.fully(lv_one_est.par, datsamp,  datpop, nis, Nis, D, samparea, L, hlist, seed = seed2, nh, T.par)
  # lv_two_est.par = lv_two$est.par.b
  lv_two_MSE = lv_two$MSE_bt
  
  
  # level-one-two Bootstrap :
  lv_one_two = one.cycle.fully(est.par, datsamp,  datpop, nis, Nis, D, samparea, L, hlist, seed = seed2, nh, T.par)
  # lv_one_two_est.par = lv_one_two$est.par.b
  lv_one_two_MSE = lv_one_two$MSE_bt
  
  
  MSE_single = lv_one_MSE
  
  MSE_double = (2*lv_one_MSE)-lv_two_MSE
  
  MSE_bc_double = lv_one_MSE+lv_one_two_MSE-lv_two_MSE
  
  
  return(list(MSE_single = MSE_single
              , MSE_double =  MSE_double
              , MSE_bc_double = MSE_bc_double))
}


MSE_Fully <- function(est.par, datsamp,  datpop, nis, Nis, D, samparea, L= 200, hlist, B){
  nh = length(hlist)
  T.par = paste0("h",1:nh) 
  
  MSE_ys <- lapply(1:B, function(i){
    EBP_fully(est.par, datsamp,  datpop, nis, Nis, D, samparea, L, hlist, nh, T.par)
  }) 
  
  denom = length(MSE_ys)
  
  MSE_single = lapply(1:denom, function(i){MSE_ys[[i]]$MSE_single})  %>% Reduce("+",.)
  MSE_double = lapply(1:denom, function(i){MSE_ys[[i]]$MSE_double})  %>% Reduce("+",.)
  MSE_bc_double = lapply(1:denom, function(i){MSE_ys[[i]]$MSE_bc_double})  %>% Reduce("+",.)
    
  MSE_single = MSE_single/denom
  MSE_double = MSE_double/denom
  MSE_bc_double = MSE_bc_double/denom
  # numer/denom
  # MSE_ys %>% Reduce("+",.) %>% ./length(MSE_ys)
  rm(datsamp,datpop); gc()
  return( list(MSE_single = MSE_single
               , MSE_double = MSE_double
               , MSE_bc_double = MSE_bc_double) )
}

###################################



## Area 1 : Leading term estimator:


LeadTerm = function(pars,datpop, datsamp){
  hat.alpha = pars["alpha"]; hat.beta0 = pars["beta0"]; hat.beta1 = pars["beta1"]; hat.delta = pars["delta"]; 
  datnonsamp = datpop %>% filter(Iunits==FALSE)
  
  nis = tapply(datsamp$yij, datsamp$area, length)
  Nis = tapply(datpop$yij, datpop$area, length)
  
  u.inv.shape = (hat.alpha*tapply(datsamp$yij, datsamp$area, length))+hat.delta  ## hat.alpha is added!!!!!!!! 03262022 Previousdated version is wrong.
  u.inv.scale = tapply(datsamp$yij*exp(hat.beta0+hat.beta1*datsamp$x), datsamp$area, sum) + hat.delta
  Eu.inv.cond = u.inv.scale/(u.inv.shape-1)
  Vu.inv.cond =  (Eu.inv.cond^2)/(u.inv.shape-2)
  Eu2.inv.cond = Vu.inv.cond+(Eu.inv.cond)^2
  
  EV = hat.alpha * tapply( exp(-2*(hat.beta0+hat.beta1*datnonsamp$x)), datnonsamp$area, sum)*  Eu2.inv.cond
  VE = (tapply(hat.alpha * exp(-hat.beta0-hat.beta1*datnonsamp$x), datnonsamp$area, sum))^2*Vu.inv.cond
  leadterm = (EV+VE)/(Nis^2)
  return(leadterm)
}



LeadTerm_b <- function(est.par, datsamp, datpop, nis, Nis, D){
  
  ## Genereate Bootsrap Parameter estmates:
  hat.alpha = est.par["alpha"];  hat.delta = est.par["delta"]; 
  hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"]
  
  # Boostrap sample:
  usD.b <- rgamma(D, shape = hat.delta, rate = hat.delta) %>% setNames(1:D)
  betaij.b <- exp(hat.beta0  + hat.beta1*datsamp$x)*usD.b[datsamp$area]
  yij.b <- rgamma(length(betaij.b), shape = hat.alpha, rate = betaij.b)
  
  # Bootstrap parameter estimates: 
  datsamp.b = datsamp; datsamp.b$yij = yij.b;
  ini.values.b = initials.propsed(datsamp.b)
  est.par.b = optim(ini.values.b, fn = loglik.pars, datsamp = datsamp.b)$par
  
  
  Ld.M1hat.b = LeadTerm(est.par.b,datpop, datsamp)
  
  rm(datpop, datsamp, datsamp.b);gc(); ##updated 08192021
  
  return(Ld.M1hat.b)
}


LeadTerm_BT <- function(est.par, datsamp, datpop, nis, Nis, D, B){
  
  yetaboot <- sapply(1:B, function(b){LeadTerm_b(est.par, datsamp, datpop, nis, Nis, D)})
  colnames(yetaboot) <-1:B 
  
  return(yetaboot)
  
}





# We decied not to use this method
# #######################################################################################################################
# 
# 
# EBP_jack = function(i, datsamp, datpop, nis, Nis, D, samparea, L){
#   
#   datsamp.i = datsamp[datsamp$area!=samparea[i],]
#   ini.values.i = initials.propsed(datsamp.i)
#   est.par.i = optim(ini.values.i, fn = loglik.pars, datsamp = datsamp.i)$par
#   
#   yNs = replicate(L, onecycle.ebs.pred(est.par.i, datsamp, datpop, nis, Nis, D))
#   
#   yNs[datpop$Iunits==1] <- datsamp$yij
#   rm(datpop, datsamp);gc();
#   return(t(yNs))
# }
# 
# MSE_Jack <- function(datsamp,  datpop, nis, Nis, D, samparea, L, hlist, EB){
#   nh = length(hlist)
#   yetajack <- sapply(1:D, function(d){EBP_jack(d, datsamp, datpop, nis, Nis, D, samparea, L)})
#   
#   colnames(yetajack) <-1:D
#   
#   yetajack <- yetajack %>% cbind.data.frame(
#     area = rep(datpop$area, each = L),
#     rep = rep(1:L, times = nrow(datpop))
#   )
#   
#   theta_rep_jack <- lapply(1:D, function(i){
#     yetajack %>% dplyr::select(i, area, rep) %>% 
#       setNames(c("theta", "area", "rep")) %>%    
#       group_by(area, rep) %>%  
#       summarise(h1 = hlist[["h1"]](theta), h2 = hlist[["h2"]](theta), h3 = hlist[["h3"]](theta), h4 = hlist[["h4"]](theta), .groups = 'drop') %>%
#       # mutate(h1 = hlist[["h1"]](theta), h2 = hlist[["h2"]](theta), h3 = hlist[["h3"]](theta), h4 = hlist[["h4"]](theta)) %>%  # h8 added 08192021
#       # group_by(area, rep) %>%  # muted 20220216
#       # summarise_at(vars(contains("h")), mean) %>%
#       mutate(boot = i) %>%
#       return(.)
#   }) %>% do.call("rbind", .)
#   
#   rm(yetajack);gc();
#   
#   
#   M1hat_jack <- theta_rep_jack %>% group_by(area,boot) %>% summarise_all(var) %>% dplyr::select(-rep)
#   
#   M1hat.jk.Rao = sapply(1:nh, function(i){
#     colnames = paste0("h",i)
#     bias = (D-1)* tapply(pull(M1hat_jack[,colnames]), M1hat_jack$area,mean)
#     bias
#   }) %>% as.data.frame(.) %>% setNames(paste0("h",1:nh))
#   gc()
#   
#   M1hat.jk.Lohr = sapply(1:nh, function(i){
#     colnames=paste0("h",i)
#     M1hat_i = M1hat_jack%>%dplyr::select(area,boot,colnames)%>%filter(.,boot!=area)
#     bias = tapply(pull(M1hat_i[,colnames]), M1hat_i$area,sum)
#     # bias = tapply(pull(M1hat_i[,colnames])-rep(pull(M1hat[,colnames]),each =(sampnum-1)), M1hat_i$area,sum)
#   })%>%as.data.frame(.)%>%setNames(paste0("h",1:nh))
#   gc()
#   
#   ## bootstrap EBUP ####
#   estpssampa_jack <- theta_rep_jack %>% dplyr::select(-rep) %>% group_by(area, boot) %>%
#     summarise_all(mean) %>% arrange(boot, area) %>% ungroup()
#   rm(theta_rep_jack);gc();
#   
#   ## Boot M2 ####
#   M2hat_jack <- sapply(1:nh, function(hi){
#     colname <- paste0("h", hi)
#     tapply((unlist(estpssampa_jack[, colname]) - (EB[, colname]))^2, estpssampa_jack$area, mean)
#   })
#   colnames(M2hat_jack)=paste0("h",1:nh)
#   
#   M2hat.jack <- as.data.frame(M2hat_jack) %>%mutate(area=samparea) %>%.[,c(paste0("h",1:nh))]
#   
#   gc()
#   return(list(M1hat.jk.Rao = M1hat.jk.Rao,
#               M1hat.jk.Lohr = M1hat.jk.Lohr,
#               M2hat.jk = M2hat.jack))
# }
# 


################First version of MSE.fully single -bootstrap only!!! #############
# 
# EBP_fully <- function(est.par, datsamp,  datpop, nis, Nis, D, samparea, L= 200, hlist){
#   nh = length(hlist)
#   T.par = paste0("h",1:nh) 
#   
#   hat.alpha = est.par["alpha"]; hat.delta = est.par["delta"]; 
#   hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"];
#   
#   # Generate yij 
#   usD.b <- rgamma(D, shape = hat.delta, rate = hat.delta) %>% setNames(1:D)
#   betaij.b <- exp(hat.beta0 + hat.beta1*datpop$x)*usD.b[datpop$area]
#   yij.b <- rgamma(length(betaij.b), shape = hat.alpha, rate = betaij.b)  
#   
#   # Bootstrap ver. Parameters:
#   datpop.b = datpop; datpop.b$yij = yij.b;
#   TR.b = sapply(1:nh, function(i){ 
#     tapply(datpop.b$yij, datpop.b$area, hlist[[i]])}) %>%
#     as.data.frame(.) %>% setNames(T.par) 
#   rm(datpop, datsamp); gc();
#   
#   # Bootstrap ver. predictors:
#   datsamp.b = datpop.b %>% filter(Iunits==1)
#   ini.values.b = initials.propsed(datsamp.b)
#   est.par.b = optim(ini.values.b, fn = loglik.pars, datsamp = datsamp.b)$par
#   est.par.b["alpha"] = abs(est.par.b["alpha"]); est.par.b["delta"] = abs(est.par.b["delta"]) ;
#   
#   EB.b = ebs.pred(est.par.b, datsamp.b, datpop.b, nis, Nis, D, L , hlist, nh)$estpssampa %>% .[,T.par]
#   
#   MSE_bt = ((EB.b - TR.b)^2) 
#   
#   rm(datpop.b, datsamp.b);gc();
#   
#   return(MSE_bt)
# }
# 
# 
# MSE_Fully <- function(est.par, datsamp,  datpop, nis, Nis, D, samparea, L= 200, hlist, B){
#   
#   MSE_ys <- lapply(1:B, function(i){
#     EBP_fully(est.par, datsamp,  datpop, nis, Nis, D, samparea, L, hlist)
#   }) 
#   
#   numer = MSE_ys %>% Reduce("+",.)
#   denom = length(MSE_ys)
#   
#   # numer/denom
#   
#   # MSE_ys %>% Reduce("+",.) %>% ./length(MSE_ys)
#   rm(datsamp,datpop); gc()
#   return(numer/denom)
# }

