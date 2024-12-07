
## Generate data sets 

gen.dat = function(D, alpha, mean.u, var.u, Ni, beta0, beta1, xij){
  
  Nis <- rep(Ni, D) %>% setNames(1:D)
  areafacpop <- rep(1:D, Nis)
  
  # Generate stratindpop. stratum1 =n1 stratume2 = n2
  N_Stratum = round(D/2)
  stratindarea <- c( rep(1, N_Stratum),
                     rep(2 , D-N_Stratum)) %>% setNames(1:D)   # names(stratindarea) <- 1:D
  stratindpop <- stratindarea[areafacpop]
  
  # Generate area random effect :usD
  delta = mean.u/var.u # changed 20220208
  # delta = mean.u*tau # E(u_i)= 0.5
  usD <- rgamma(D, shape = delta, rate = delta) %>% setNames(1:D)
  
  # Generate yij 
  # xij <-  runif(sum(Nis), 0, 2) ## muted 20220222
  betaij <- exp(beta0  + beta1*xij)*usD[areafacpop]
  yij <- rgamma(length(betaij), shape = alpha, rate = betaij)
  
  # n1i = 0.10*Ni; n2i = 0.20*Ni;
  n1i = 25; n2i = 50 # 08122024
  samp.df = data.frame(area = areafacpop,strat = stratindpop, ind = 1:length(areafacpop))
  
  sampind = lapply(1:D, function(d){
    
      area.d = samp.df[samp.df$area==d, ]
      
      if(unique(area.d$strat)==1){
       sampind.d = sample(samp.df[samp.df$area==d, "ind"], n1i) %>% sort(.)
       
      }else{
       sampind.d = sample(samp.df[samp.df$area==d, "ind"], n2i) %>% sort(.)
      }
    }) %>% Reduce("append",.)
  Iunits = samp.df$ind %in% sampind

  datpop <- data.frame(area = areafacpop, strat = stratindpop, x = xij, 
                       yij = yij, Iunits = Iunits )
  
  datsamp <- datpop %>% filter(Iunits==1)
  
  true.par = c(alpha, beta0, beta1, delta)
  names(true.par) = c("alpha", "beta0", "beta1", "delta")
  
  
  D = length(unique(datpop$area))
  Nis = tapply(datpop$yij, datpop$area, length)
  nis = tapply(datsamp$yij, datsamp$area, length)
  samparea = unique(datsamp$area) ## added 02142022
  
  list(datpop = datpop, datsamp = datsamp, true.par = true.par, D = D, Nis = Nis, nis = nis, samparea = samparea)
  
}


## logloklihood for model-parameters

initials.propsed = function(datsamp){
  glme.ini = glmer(yij ~ x + (1|area), family = Gamma(link="log"),  data = datsamp)
  ini.alpha = 1/sigma(glme.ini)^2 # Moment Estimator for alpha with yij|j \in s : (mean(datsamp$yij))^2/var(datsamp$yij) # changed 20220209
  ini.beta0 = log(ini.alpha)-fixef(glme.ini)[1] 
  ini.beta1 = -fixef(glme.ini)[2] 
  ini.est.ranef = abs(unlist(ranef(glme.ini)))
  ini.delta = mean(ini.est.ranef)^2/var(ini.est.ranef)
  # ini.delta =   mean(ini.est.ranef)/var(ini.est.ranef)
  ini.values = c(ini.alpha,ini.beta0,ini.beta1, ini.delta)
  names(ini.values) = c("alpha","beta0","beta1","delta")
  return(ini.values)
}

loglik.pars = function(par, datsamp){
  hat.alpha = abs(par["alpha"])
  hat.beta0 = par["beta0"]; hat.beta1 = par["beta1"]
  hat.delta = abs(par["delta"]); 
  # hat.tau = par["tau"] # changed 20220209
  hat.betas = c(hat.beta0, hat.beta1)
  
  nis = tapply(datsamp$yij, datsamp$area, length)
  
  
  sum.lik= sapply(1:length(unique(datsamp$area)), function(d){
    ni.d = nis[d]
    datsamp.d = datsamp[datsamp$area==d,]
    sum.xij.d = cbind(1, datsamp.d$x) %>% apply(.,2,sum)
    
    sum.lik.i = hat.delta*log(hat.delta) + (hat.alpha-1)*sum(log(datsamp.d$yij)) + hat.alpha*sum(sum.xij.d*hat.betas)+lgamma((ni.d*hat.alpha)+hat.delta) -  ## changed 07112024 from log(gamma) to lgamma for numerical issue.
      lgamma(hat.delta) - (ni.d * lgamma(hat.alpha)) - ((ni.d*hat.alpha)+hat.delta)*log(sum(datsamp.d$yij*exp(hat.beta0+(hat.beta1*datsamp.d$x)))+hat.delta)
    return(sum.lik.i)
  }) %>% sum(.)
  
  return(-sum.lik)
  
}

# Prediction : Empirical Bayes Prediction

onecycle.ebs.pred = function(est.par, datsamp, datpop, nis, Nis, D){
  
  hat.alpha = est.par["alpha"]; hat.delta = est.par["delta"]; 
  hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"];

  # hat.tau = est.par["tau"]
  
  # Generate conditional Ui
  est.rate = datsamp %>% mutate(est.rate =  yij*exp(hat.beta0+hat.beta1*x)) %>% group_by(area) %>%
            summarise_at(vars(est.rate), sum) %>% dplyr::select(est.rate) %>% unlist(.)
  # est.rate2 = tapply((datsamp$yij)*(exp(hat.beta0+hat.beta1*datsamp$x)), datsamp$area,sum)
  
  hat.u.shape = ((nis*hat.alpha)+hat.delta); hat.u.rate = hat.delta + est.rate;
  cond.ui = rgamma(length(nis), shape = hat.u.shape, rate = hat.u.rate)  %>% setNames(1:D)
  
  # Generate pop yij
  condui.pop = cond.ui[as.character(datpop$area)]
  datpop$condui = condui.pop
  hat.y.rate = datpop %>% mutate(est.rat = exp(hat.beta0+hat.beta1*x)*condui) %>% dplyr::select(est.rat) %>% unlist(.)
  # hat.y.rate2 = exp(hat.beta0+hat.beta1*datpop$x)*datpop$condui
  pred.yij = rgamma(length(hat.y.rate), shape = hat.alpha, rate = hat.y.rate)
  
  pred.yij[datpop$Iunits==TRUE] = datsamp$yij
  
  return(pred.yij)
}



ebs.pred = function(est.par, datsamp, datpop, nis, Nis, D, L = 300, hlist, nh){
  
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
  
  return(list(estpssampa = estpssampa, M1hat = M1hat))
}




###########HZ method ################################


HZ.est.par = function(glme.ini.log, datsamp){
  
  HZ.pars = c(fixef(glme.ini.log), sqrt(unlist(VarCorr(glme.ini.log))), 1/sigma(glme.ini.log)^2)
  names(HZ.pars) = c("beta0", "beta1","phi.u","phi")
  
  return(HZ.pars)
}

one.cycle.hz.pred = function(HZ.est.pars, datsamp, datpop, nis, Nis, D, L2,hlist, nh){
  beta0.hz = HZ.est.pars["beta0"];  beta1.hz =  HZ.est.pars["beta1"];
  phi.u.hz = HZ.est.pars["phi.u"];  phi.hz =  HZ.est.pars["phi"];
  
  #one.cycle is to generate u_d ~ normal(0, .hat{phi.u})
  uD = rnorm(D, mean = 0, sd = phi.u.hz ) %>% setNames(1:D)
  new.datpop = datpop %>% mutate(uD = uD[area], mu.hat.dj = exp(beta0.hz+beta1.hz*x+uD), samp.unit.pyij = dgamma(yij, shape = phi.hz, rate =phi.hz/mu.hat.dj))## use loglink
  
  #Denominator: 
  denom = new.datpop %>% filter(Iunits==1) %>% group_by(area) %>% summarise_at(vars(samp.unit.pyij), prod) #should be prod
  
  uD.yij = replicate(L2, rgamma(nrow(new.datpop), shape = phi.hz, rate = phi.hz/(new.datpop$mu.hat.dj)))
  uD.yij[datpop$Iunits==1, ] = datsamp$yij
  
  T.par = paste0("h",1:nh)
  samparea = unique(datsamp$area)
  
  theta_rep <- sapply(1:nh, function(i) {                                     # changed 08192021 for gini coefficient  abs(y) and exp(y) in gini
    apply(uD.yij, 2 , function(y){tapply((y), datpop$area, hlist[[i]])})}) %>%  as_tibble() %>%
    setNames(paste0("h", 1:nh)) %>% 
    cbind(area = rep(samparea, times = L2)) %>% group_by(area) %>% summarise_all(mean)
  
  #numerator:
  numer = sapply(1:nh, function(.){
    Temp.par = paste0("h", .)
    unlist(theta_rep[, Temp.par])*unlist(denom[,"samp.unit.pyij"])})
  colnames(numer) = T.par
  denom = denom %>% dplyr::select(samp.unit.pyij) %>% unlist(.)
  
  rm(theta_rep, uD.yij, datpop, datsamp, new.datpop);gc();
  
  return(list(denom = denom, numer = numer))
}

HZ.EBP = function(HZ.est.pars, datsamp, datpop, nis, Nis, D, L1=100, L2=200,hlist, nh){
  
  beta0.hz = HZ.est.pars["beta0"];  beta1.hz =  HZ.est.pars["beta1"];
  phi.u.hz = HZ.est.pars["phi.u"];  phi.hz =  HZ.est.pars["phi"];
  
  HZ.SIR = lapply(1:L1, function(.){one.cycle.hz.pred(HZ.est.pars, datsamp, datpop, nis, Nis, D, L2,hlist, nh)})
  
  denom.d = lapply(1:length(HZ.SIR), function(.){HZ.SIR[[.]][["denom"]]}) %>% Reduce("+",.)/length(HZ.SIR)
  nunom.d = lapply(1:length(HZ.SIR), function(.){HZ.SIR[[.]][["numer"]]}) %>% Reduce("+",.)/length(HZ.SIR)
  
  rm(HZ.SIR, datpop, datsamp);gc();
  return(nunom.d/denom.d)
  
}

HZ.PLUG = function(glme.ini.log, HZ.est.pars, datsamp, datpop, nis, Nis, D,hlist, nh){
  est.area.ranef = ranef(glme.ini.log) %>% unlist(.) %>% setNames(1:D)
  beta0.hz = HZ.est.pars["beta0"];  beta1.hz =  HZ.est.pars["beta1"];
  phi.u.hz = HZ.est.pars["phi.u"];  phi.hz =  HZ.est.pars["phi"];
  new.datpop = datpop %>% mutate(est.ranef = est.area.ranef[area], mu.tilde.dj = exp(beta0.hz+beta1.hz*x+est.ranef))
  new.datpop$mu.tilde.dj[new.datpop$Iunits==1] = datsamp$yij
  
  HZ.plug = sapply(1:nh, function(i){ 
    tapply(new.datpop$mu.tilde.dj, new.datpop$area, hlist[[i]])}) %>%
    as.data.frame(.) %>% setNames(paste0("h",1:nh))
  rm(new.datpop,datpop, datsamp); gc();
  
  return(HZ.plug)
  
}

one.cycle.HZ.marginal = function(glme.ini.log, HZ.est.pars, datsamp, datpop, nis, Nis, D,hlist, nh){
  est.area.ranef = ranef(glme.ini.log) %>% unlist(.) %>% setNames(1:D)
  beta0.hz = HZ.est.pars["beta0"];  beta1.hz =  HZ.est.pars["beta1"];
  phi.u.hz = HZ.est.pars["phi.u"];  phi.hz =  HZ.est.pars["phi"];
  
  new.datpop = datpop %>% mutate(est.ranef = est.area.ranef[area], 
                                 mu.tilde.dj = exp(beta0.hz+beta1.hz*x+est.ranef))
  one.y.marginal = rgamma(nrow(datpop), shape = phi.hz, rate = phi.hz/(new.datpop$mu.tilde.dj))
  
  one.y.marginal[new.datpop$Iunits==1] = datsamp$yij
  
  rm(new.datpop, datpop, datsamp); gc();
  return(one.y.marginal)
}

HZ.Marginal = function(glme.ini.log, HZ.est.pars, datsamp, datpop, nis, Nis, D,hlist, nh, L=200){
  
  y.marginals = replicate(L, one.cycle.HZ.marginal(glme.ini.log, HZ.est.pars, datsamp, datpop, nis, Nis, D,hlist, nh))
  
  T.par = paste0("h",1:nh)
  samparea = unique(datsamp$area)
  
  theta_rep <- sapply(1:nh, function(i) {                                     # changed 08192021 for gini coefficient  abs(y) and exp(y) in gini
    apply(y.marginals, 2 , function(y){tapply((y), datpop$area, hlist[[i]])})}) %>%  as_tibble() %>%
    setNames(paste0("h", 1:nh)) %>% 
    cbind(area = rep(samparea, times = L))
  
  HZ.marginal <- theta_rep %>% group_by(area) %>% summarise_at(vars(h1,h2,h3,h4),mean) %>% dplyr::select(h1,h2,h3,h4) %>% as.data.frame(.)
  
  rm(datpop, datsamp, theta_rep, y.marginals); gc();
  
  return(HZ.marginal)  
  
}

######### MSE Calculation ################

Tmse = function(estpssampa,popsampa,nh){
  Tmse.res = sapply(1:nh, function(i){
    colname=paste0("h",i)
    ((estpssampa[,colname])-popsampa[,colname])^2
  })
  colnames(Tmse.res) = paste0("h",1:nh)
  
  Tmse.res  
  
}

Summary_MSE.dfs = function(named_vec){
  # named_vec = avg.vec
  T.name = names(named_vec)
  val = named_vec
  
  MSE_T = grep("MSE",T.name)
  method = gsub("(.*)_h.*","\\1",T.name)
  h_id = gsub(".*_(h.*)", "\\1",T.name)
  
  long.form = data.frame(h_id=h_id, method= method,val = val) %>% 
    mutate(row_ID = 1:length(named_vec), MSE_T = ifelse(row_ID %in% MSE_T, 1, 0))  %>%
    dplyr::select(MSE_T, h_id, method, val)
  
  Methods = c("h_id","MSE.True","MSE.TRcond","MSE.noBC","MSE.add","MSE.mult","MSE.hm","MSE.comp",
              # "MSE.noBC200","MSE.add200","MSE.mult200","MSE.hm200","MSE.comp200",
              "MSE.single","MSE.double","MSE.bc_double") # changed 08192021
  MSE.df = filter(long.form, MSE_T ==1)[,-1] %>% spread(.,method,val) %>% dplyr::select(one_of(Methods)) 
  
  MSE.df.Ucond = MSE.df;  MSE.df.cond = MSE.df;
  MSE.df.Ucond[,-c(1,2,3)] = ((MSE.df.Ucond[,-c(1,2,3)]- MSE.df.Ucond[,"MSE.True"])/MSE.df.Ucond[,"MSE.True"])*100
  MSE.df.cond[,-c(1,2,3)] = ((MSE.df.cond[,-c(1,2,3)]- MSE.df.cond[,"MSE.TRcond"])/MSE.df.Ucond[,"MSE.TRcond"])*100
  
  MSE.df.Ucond$BEST = colnames(MSE.df.Ucond[,Methods[-c(1,2,3)]])[apply(abs(MSE.df.Ucond[,Methods[-c(1,2,3)]]),1,which.min)]
  MSE.df.cond$BEST = colnames(MSE.df.cond[,Methods[-c(1,2,3)]])[apply(abs(MSE.df.cond[,Methods[-c(1,2,3)]]),1,which.min)]
  
  CI_T = grepl("Naive|Norm",T.name)
  CI_95 = grepl(".95",T.name)
  
  # long.form.CI = data.frame(h_id=h_id, method= method,val = val) %>% 
  #   mutate(row_ID = 1:length(named_vec), CI_T = ifelse(row_ID %in% CI_T, 1, 0))  %>%
  #   dplyr::select(CI_T, h_id, method, val)
  
  long.form.CI = data.frame(h_id=h_id, method= method,val = val) %>% 
    mutate(row_ID = 1:length(named_vec), CI95_T = CI_T*CI_95)  %>%
    dplyr::select(CI95_T, h_id, method, val)
  
  Methods.CI = c("h_id","Norm.noBC","Norm.add","Norm.mult","Norm.comp","Norm.hm",
                 "Norm.noBC200","Norm.add200","Norm.mult200","Norm.comp200","Norm.hm200","Norm.single","Naive")  
  CI.df = filter(long.form.CI, CI95_T ==1) %>% spread(., method, val)
  
  # EB.bias = filter(long.form, MSE_T!=1) %>% spread(., method, val) %>% dplyr::select(one_of(Methods.EB))
  # EB.bias[,-c(1,2)] = ((EB.bias[,-c(1,2)]-EB.bias[,"TR"])/EB.bias[,"TR"])*100
  # EB.bias$BEST = colnames(EB.bias[,Methods.EB[-c(1,2)]])[apply(abs(EB.bias[,Methods.EB[-c(1,2)]]),1,which.min)]
  # 
  return(list( MSE.df.Ucond = MSE.df.Ucond, MSE.df.cond = MSE.df.cond, CI.df = CI.df))
  
}

Test_statistics = function(variable_names, res_list){
  res_dfs = sapply(1:length(res_list), function(.){
            temp.name = variable_names[grep("lead.bias",variable_names)]
            temp.vec = res_list[[.]][grepl("lead.bias",variable_names),] %>% apply(.,1,mean)
            names(temp.vec) = temp.name
            temp.vec
          })
  
 Test_stat = apply( res_dfs,1,mean)/sqrt(apply( res_dfs,1,var)/length(res_list))
 Test_stat_mat = matrix(Test_stat, ncol=2, byrow = T)
 rownames(Test_stat_mat)   = paste0("h",1:4)
 colnames(Test_stat_mat)   = c("True.bias","Est.bias")
 return(Test_stat_mat) 
}






