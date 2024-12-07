
## Generate data sets 

gen.dat = function(D, alpha, mean.u, var.u, Ni, beta0, beta1, xij,a,b, phi, pop.model, wt.model){
  
  Nis <- rep(Ni, D) %>% setNames(1:D)
  areafacpop <- rep(1:D, Nis)
  
  # Generate stratindpop. stratum1 =n1 stratume2 = n2
  N_Stratum = round(D/2)
  stratindarea <- c( rep(1, N_Stratum),
                     rep(2 , D-N_Stratum)) %>% setNames(1:D)   # names(stratindarea) <- 1:D
  stratindpop <- stratindarea[areafacpop]
  n1i = 0.025*Ni; n2i = 0.05*Ni; ## Changed 07262024
  
  ## Population Model
  delta = mean.u/var.u # changed 20220208
  if(pop.model == "Gam"){
    
    # Generate area random effect :usD
   
    # delta = mean.u*tau # E(u_i)= 0.5
    usD <- rgamma(D, shape = delta, rate = delta) %>% setNames(1:D)
    
    # Generate yij 
    # xij <-  runif(sum(Nis), 0, 2) #should be fixed before running simulation loop!!!  ## muted 20220222
    betaij <- exp(beta0  + beta1*xij)*usD[areafacpop]
    yij <- rgamma(length(betaij), shape = alpha, rate = betaij)
    
  }else {
    
    nuD <- rnorm(D, sd = phi) %>% setNames(1:D)
    # Generate yij 
    # xij <-  runif(sum(Nis), 0, 2) #should be fixed before running simulation loop!!!  ## muted 20220222
    muij <- exp(0.5+ 0.05*xij + nuD[areafacpop])
    yij <- rgamma(length(muij), shape = alpha, rate = alpha/muij)
    
  }
  nis <- rep(n1i, length(stratindpop))
  nis[stratindpop == 2] <- n2i
  
  ## Weight Model
  
  if(wt.model == "Exp" ){
    
    # Generate pij
    deltaij <- rgamma(sum(Nis), delta, delta)
    zij <- exp((a*xij+b*yij) + deltaij/30)
    zijUs <-  tapply(zij, areafacpop , sum)
    
    # ; nis[stratindpop == 3] <- 15
    probij <- nis*zij/(zijUs[areafacpop])
    
  }else {
    
    # Generate pij
    deltaij <- rnorm(sum(Nis), 0, 1)
    b_new = ifelse(b ==0, 0, exp(b))
    probij_unscaled <- 1 + xij + b_new*yij + deltaij/20 ## changed for linear model weighting 07112024 
    probijUs <- tapply(probij_unscaled, areafacpop, sum)
    
    probij <- nis*probij_unscaled/(probijUs[areafacpop])
    
  }
  
  
  #Generate unit indicator
  Iunits <- unlist(lapply(1:D, function(i){UPsystematic(probij[areafacpop == i],  eps = 1e-200)} )) # VVVVV
  # Iunits[!(areafacpop %in% samparea)] <- 0
  # tapply(Iunits, areafacpop, sum)
  
  datpop <- data.frame(area = areafacpop, strat = stratindpop, x = xij, 
                       yij = yij, Iunits = Iunits, probij = probij )
  
  datsamp <- datpop %>% filter(Iunits==1)
  
  true.par = c(alpha, beta0, beta1, delta, a, b)
  names(true.par) = c("alpha", "beta0", "beta1", "delta", "a","b")
  
  
  D = length(unique(datpop$area))
  Nis = tapply(datpop$yij, datpop$area, length)
  nis = tapply(datsamp$yij, datsamp$area, length)
    
  list(datpop = datpop, datsamp = datsamp, true.par = true.par, D = D, Nis = Nis, nis = nis)
  
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
  
  hat.alpha = est.par["alpha"]
  hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"]
  hat.delta = est.par["delta"]; 
  # hat.tau = est.par["tau"]
  
  # Generate conditional Ui
  est.rate = datsamp %>% mutate(est.rate =  yij*exp(hat.beta0+hat.beta1*x)) %>% group_by(area) %>% summarise_at(vars(est.rate), sum) %>% dplyr::select(est.rate) %>% unlist(.)
  hat.u.shape = ((nis*hat.alpha)+hat.delta); hat.u.rate = hat.delta + est.rate;
  cond.ui = rgamma(length(nis), shape = hat.u.shape, rate = hat.u.rate)  %>% setNames(1:D)
  
  # Generate pop yij
  condui.pop = cond.ui[as.character(datpop$area)]
  datpop$condui = condui.pop
  hat.y.rate = datpop %>% mutate(est.rat = exp(hat.beta0+hat.beta1*x)*condui) %>% dplyr::select(est.rat) %>% unlist(.)
  
  pred.yij = rgamma(length(hat.y.rate), shape = hat.alpha, rate = hat.y.rate)
  
  pred.yij[datpop$Iunits==TRUE] = datsamp$yij
  
  return(pred.yij)
}



ebs.pred = function(est.par, datsamp, datpop, nis, Nis, D, L  , hlist, nh){
  
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

#######################################################################################################################
#### ADDED for Informative Sampling:
#######################################################################################################################
sse.ab <- function(pars, D, datsamp){
  hat.a <- pars[1]; hat.b <- (pars[2]); hat.kappai<-pars[(1:D)+2]
  names(hat.kappai) = 1:D
  wij <- 1/datsamp$prob
  sum((wij - hat.kappai[datsamp$area]*exp(hat.a*datsamp$x-hat.b*datsamp$yij))^2) # changed -hat.a to hat.a 11202024
}

wt_par_ab = function (datsamp, D){
  D = datsamp$area %>% unique(.) %>% length(.)
  ini.val.ab = lm(log(1/probij)~0+x+yij+as.factor(area), data = datsamp)
  coef.ab  = coef(ini.val.ab)
  coef.ab[(1:D)+2] <- exp(coef(ini.val.ab)[(1:D)+2])
  coef_ab = optim(coef.ab, fn = sse.ab, D = D, datsamp = datsamp, control = list(maxit = 10000))$par
  coef_ab
}

onecycle.ebs.pred.INFOR = function(est.par, datsamp, datpop, nis, Nis, D, a, b){
  
  hat.alpha = est.par["alpha"]
  hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"]
  hat.delta = est.par["delta"]; 
  hat.a = a;
  hat.b = b;
  # hat.tau = est.par["tau"]
  
  # Generate conditional Ui
  est.rate = datsamp %>% mutate(est.rate =  yij*exp(hat.beta0+hat.beta1*x)) %>% group_by(area) %>% summarise_at(vars(est.rate), sum) %>% dplyr::select(est.rate) %>% unlist(.)
  hat.u.shape = ((nis*hat.alpha)+hat.delta); hat.u.rate = hat.delta + est.rate;
  cond.ui = rgamma(length(nis), shape = hat.u.shape, rate = hat.u.rate)  %>% setNames(1:D)
  
  # Generate pop yij
  condui.pop = cond.ui[as.character(datpop$area)]
  datpop$condui = condui.pop
  hat.y.rate = datpop %>% mutate(est.rat = exp(hat.beta0+hat.beta1*x)*condui) %>% dplyr::select(est.rat) %>% unlist(.)
  
  pred.yij = rgamma(length(hat.y.rate), shape = hat.alpha, rate = hat.y.rate +hat.b)
  
  pred.yij[datpop$Iunits==TRUE] = datsamp$yij
  
  return(pred.yij)
}


 ## Used f_{p}(yij|.)
ebs.pred.INFOR = function(est.par, datsamp, datpop, nis, Nis, D, L  , hlist, nh, a,b){
  
  yNs = replicate(L, onecycle.ebs.pred.INFOR(est.par, datsamp, datpop, nis, Nis, D,a,b))  
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


onecycle.ebs.pred.INFOR.fciyij = function(est.par, datsamp, datpop, nis, Nis, D, a, b, kappa.i){
  
  hat.alpha = est.par["alpha"]
  hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"]
  hat.delta = est.par["delta"]; 
  hat.a = a;
  hat.b = b;
  # hat.tau = est.par["tau"]
  
  # Generate conditional Ui
  est.rate = datsamp %>% mutate(est.rate =  yij*exp(hat.beta0+hat.beta1*x)) %>% group_by(area) %>% summarise_at(vars(est.rate), sum) %>% dplyr::select(est.rate) %>% unlist(.)
  hat.u.shape = ((nis*hat.alpha)+hat.delta); hat.u.rate = hat.delta + est.rate;
  cond.ui = rgamma(length(nis), shape = hat.u.shape, rate = hat.u.rate)  %>% setNames(1:D)
  
  # Generate pop yij
  u.runif = runif(nrow(datpop))
  Fciyij = function( y_ij, EsiWij, hat.alpha, hat.b, betaij){
    (EsiWij*pgamma(y_ij, hat.alpha,rate = hat.b+betaij)-pgamma(y_ij, hat.alpha,rate = betaij))/(EsiWij-1)
  }
  Fciyij <- Vectorize(Fciyij)
  Fciyij.inv <- function(x, EsiWij.i, hat.alpha, hat.b, betaij.i, u.runif.i){Fciyij(x, EsiWij.i, hat.alpha, hat.b, betaij.i)-u.runif.i}
  Fciyij.inv <- Vectorize(Fciyij.inv)
  area.df = data.frame(area=unique(datsamp$area), kappa_i = kappa.i, nis = nis, Nis = Nis, cond.ui = cond.ui)
  new.df = datpop %>% left_join(., area.df) %>% 
    mutate(betaij = exp(hat.beta0+hat.beta1*x)*cond.ui
           , EsiWij = kappa_i*exp(hat.a*x)*(1+(hat.b/betaij))^(-hat.alpha))
  
  EsiWij = new.df$EsiWij ; betaij = new.df$betaij
  pred.yij = lapply(1:nrow(new.df), function(i){
    uniroot(f = function(x){Fciyij(x, EsiWij[i], hat.alpha, hat.b, betaij[i])-u.runif[i]},interval = c(0,100))$root
  }) %>% unlist(.)
  
  pred.yij[datpop$Iunits==TRUE] = datsamp$yij
  
  return(pred.yij)
}

## Used f_{ci}(yij|xij, ui,Ii=1)
ebs.pred.INFOR.fciyij = function(est.par, datsamp, datpop, nis, Nis, D, L  , hlist, nh, a,b, kappa.i){
  
  yNs = replicate(L, onecycle.ebs.pred.INFOR.fciyij(est.par, datsamp, datpop, nis, Nis, D, a, b, kappa.i))  
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


onecycle.ebs.pred.INFOR.fp = function(est.par, datsamp, datpop, nis, Nis, D, a, b){
  
  hat.alpha = est.par["alpha"]
  hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"]
  hat.delta = est.par["delta"]; 
  hat.a = a;
  hat.b = b;
  # hat.tau = est.par["tau"]
  
  # Generate conditional Ui
  est.rate = datsamp %>% mutate(est.rate =  yij*exp(hat.beta0+hat.beta1*x)) %>% group_by(area) %>% summarise_at(vars(est.rate), sum) %>% dplyr::select(est.rate) %>% unlist(.)
  hat.u.shape = ((nis*hat.alpha)+hat.delta); hat.u.rate = hat.delta + est.rate;
  cond.ui = rgamma(length(nis), shape = hat.u.shape, rate = hat.u.rate)  %>% setNames(1:D)
  
  # Generate pop yij
  condui.pop = cond.ui[as.character(datpop$area)]
  datpop$condui = condui.pop
  hat.y.rate = datpop %>% mutate(est.rat = exp(hat.beta0+hat.beta1*x)*condui) %>% dplyr::select(est.rat) %>% unlist(.)
  
  pred.yij = rgamma(length(hat.y.rate), shape = hat.alpha, rate = hat.y.rate +hat.b)
  
  # pred.yij[datpop$Iunits==TRUE] = datsamp$yij
  
  return(pred.yij)
}


## Used f_{p}(yij|.)
ebs.pred.INFOR.fp = function(est.par, datsamp, datpop, nis, Nis, D, L, hlist, nh, a,b){
  
  yNs = replicate(L, onecycle.ebs.pred.INFOR.fp(est.par, datsamp, datpop, nis, Nis, D,a,b))  
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




#########################################################################################################################################################
#########################################################################################################################################################


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

HZ.EBP = function(HZ.est.pars, datsamp, datpop, nis, Nis, D, L1, L2,hlist, nh){
  
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

HZ.Marginal = function(glme.ini.log, HZ.est.pars, datsamp, datpop, nis, Nis, D,hlist, nh, L){
  
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

Tbias = function(estpssampa,popsampa,nh){
  Tbias.res = sapply(1:nh, function(i){
    colname=paste0("h",i)
    ((estpssampa[,colname])-popsampa[,colname])
  })
  colnames(Tbias.res) = paste0("h",1:nh)
  
  Tbias.res  
  
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
  
  # Methods = c("h_id","TMSE_B","TMSE_B.clsd","TMSE_EB","TMSE_EB.clsd", "TMSE_EB_infor" ,"TMSE_EB_HZ", "TMSE_EB_HZmarginal", "TMSE_EB_HZplug","TMSE_EB_dir") # changed 08192021
  # Methods = c("h_id","TMSE_B","TMSE_EB", "TMSE_EB_infor","TMSE_EB_infor_fp","TMSE_EB_infor_fci" ,"TMSE_EB_HZ", "TMSE_EB_HZmarginal", "TMSE_EB_HZplug","TMSE_EB_dir") # changed 07082022
  Methods = c("h_id","TMSE_B","TMSE_EB", "TMSE_EB_infor","TMSE_EB_infor_fci","TMSE_EB_HZ", "TMSE_EB_HZmarginal", "TMSE_EB_HZplug","TMSE_EB_dir") # changed 07212024
  MSE.df = filter(long.form, MSE_T ==1)[,-1] %>% spread(.,method,val) %>% dplyr::select(one_of(Methods)) 
  MSE.df$BEST = colnames(MSE.df[,Methods[-c(1:2)]])[apply(abs(MSE.df[,Methods[-c(1:2)]]),1,which.min)]
  
  # Methods.EB = c("h_id","TBIAS_EB","TBIAS_EB_infor","TBIAS_EB_infor_fp","TBIAS_EB_infor_fci","TBIAS_EB_infor","TBIAS_EB_HZ","TBIAS_EB_HZmarginal","TBIAS_EB_HZplug","TBIAS_EB_dir") 
  Methods.EB = c("h_id","TBIAS_EB","TBIAS_EB_infor","TBIAS_EB_infor_fci","TBIAS_EB_infor","TBIAS_EB_HZ","TBIAS_EB_HZmarginal","TBIAS_EB_HZplug","TBIAS_EB_dir") # changed 07212024
  EB.bias = filter(long.form, MSE_T!=1) %>% spread(., method, val) %>% dplyr::select(one_of(Methods.EB))
  # EB.bias[,-c(1,2)] = ((EB.bias[,-c(1,2)]-EB.bias[,"TR"])/EB.bias[,"TR"])*100
  EB.bias$BEST = colnames(EB.bias[,Methods.EB[-c(1)]])[apply(abs(EB.bias[,Methods.EB[-c(1)]]),1,which.min)]
  
  return(list(EB.bias = EB.bias, MSE.df = MSE.df))
  
}

##########################################################
#### 05162022 Closed Formula for Mean predictor added#####
##########################################################
ebs_clsd.pred= function(est.par, datsamp, datpop, nis, Nis, D,  hlist, nh,TR,hat.b){
  
  hat.alpha = est.par["alpha"]; hat.delta = est.par["delta"]; 
  hat.beta0 = est.par["beta0"]; hat.beta1 = est.par["beta1"];
  
  
  Eui.inv.cond = datsamp %>% mutate(cond.ui = yij*exp(hat.beta0+hat.beta1*x)) %>% 
    group_by(area) %>% summarise_at(vars(cond.ui),sum) %>%
    left_join(., as_tibble(data.frame(area = 1:D, ni = nis))) %>%
    mutate(Eui.cond = (cond.ui + hat.delta)/(ni*hat.alpha+hat.delta-1))
  
  clsd_pred = datpop %>% left_join(., Eui.inv.cond, by = "area") %>%
    mutate(margE.yij = hat.alpha * exp(-hat.beta0 - hat.beta1*x)*Eui.cond,
           margE.yij = ifelse(Iunits==TRUE, yij, margE.yij)) %>% 
    group_by(area) %>% summarise_at(vars(margE.yij), mean) %>%data.frame(.) %>%
    setNames(c("area", "h1"))
  
  
  return(clsd_pred = clsd_pred)
}

#############################################################################
