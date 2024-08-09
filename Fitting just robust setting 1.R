

library(dplyr)
library(sf)
library(sp)

library(TMB)
library(starve)
compile("nngp_updated_model.cpp")#,"&> C:/Users/mcdonaldra/Documents/errors.txt")
dyn.load(dynlib("nngp_updated_model"))
compile("nngp_updated_model_nonstat.cpp")#,"&> C:/Users/mcdonaldra/Documents/errors.txt")
dyn.load(dynlib("nngp_updated_model_nonstat"))
compile("nngp_updated_subsim_nonstat.cpp")
dyn.load(dynlib("nngp_updated_subsim_nonstat"))
compile("nngp_updated_robust.cpp")
dyn.load(dynlib("nngp_updated_robust"))

library(ggplot2)
library(optimr)

# Logit function
logitp=function(p){log(p/(1-p))}
# Inverse logist function
logitpi=function(t){exp(t)/(1+exp(t))}

#Create a 100X100km square
x_coord<-c(175,225,225,175)
y_coord<-c(175,175,225,225)
xy<-data.frame(x=x_coord,y=y_coord)
sf_poly<-st_as_sf(xy,coords=c("x","y")) %>% st_bbox() %>%
  st_as_sfc() %>% st_as_sf()

sf_grid<-st_make_grid(sf_poly,cellsize=1)
sf_centroids<-st_as_sf(st_centroid(sf_grid))
sf_centroids$ID<-1:length(sf_centroids$x)
sf_grid2<-st_as_sf(sf_grid)

sf_grid2$ID<-1:nrow(sf_grid2)

bot_left_ids<-rep(1:20,20)
for (row in 2:20){
  bot_left_ids[((row-1)*20+1):(row*20)]<-bot_left_ids[((row-1)*20+1):(row*20)]+((row-1)*50)
}

n_vess<-40
n_t<-2
strat_n_t<-5
big_n_t<-20
n_k<-3
n_k2<-4
n_sets<-70

props_alloc_sets<-round(c(5/(5+7+10),7/(5+7+10),10/(5+7+10))*n_sets)

nrep<-10
rob_list<-list()
for (sub in 1:20){
  load(paste0("sub_data_list_set_1_",sub,"RData"))
  for (n in 1:nrep){
    fit_data_rob<-sub_data_list[[1]][[n]]
    fit_data_rob$strat_d<-rep(c(rep(0,16),rep(1,22),rep(2,32)),big_n_t)
    fit_data_rob$robcode<-2
    fit_data_rob$tc<-c(50,8,4)
    
    par_list_rob<-list(theta=logitp(0.5),
                       vess_eff_t=rep(0,length(unique(fit_data_rob$vess_id))),
                       vess_eff_nt=rep(0,length(unique(fit_data_rob$vess_id))),
                       pre_pars_rands=array(data=c(-10,logitp((0.99+1)/2),-1,-8,logitp((0.99+1)/2),-1),dim=c(3,2)),
                       pre_pars_omegas=array(data=c(log(0.1),log(5),0,log(0.1),log(5),0),dim=c(3,2)),
                       log_sigma_vess_nt=-0.5,
                       log_sigma_vess_t=-0.5,
                       rands=matrix(rep(0,big_n_t*2),ncol=2),
                       pg_re=array(rep(0,(length(fit_data_rob$pg_edges)+2)*big_n_t*2),dim=c((length(fit_data_rob$pg_edges)+2),big_n_t,2)),
                       tg_re=matrix(rep(0,length(fit_data_rob$tg_edges)*2),ncol=2))
    random_rob<-c("vess_eff_t","vess_eff_nt","rands","pg_re","tg_re")
    
    map<-list(pre_pars_omegas=as.factor(array(c(1,2,NA,3,4,NA),dim=c(3,2))),
              pre_pars_rands=as.factor(array(c(5,NA,6,7,NA,8),dim=c(3,2))))
    
    obj_rob<-MakeADFun(fit_data_rob,par_list_rob,random=random_rob,DLL="nngp_updated_robust",map=map,silent=F)
    Opt_rob<-optimx::optimr(obj_rob$par,obj_rob$fn,obj_rob$gr,control=list(maxit=100000),method="nlminb")
    while (Opt_rob$message=="iteration limit reached without convergence (10)"){
      obj_rob$par<-obj_rob$env$last.par.best[-which(names(obj_rob$env$last.par.best)%in% random_rob)]
      Opt_rob<-optimx::optimr(obj_rob$par,obj_rob$fn,obj_rob$gr,control=list(maxit=100000),method="nlminb")
    }
    rep_rob<-sdreport(obj_rob)
    
    rob_list[[n]]<-list()
    rob_list[[n]]$obj<-obj_rob
    rob_list[[n]]$Opt<-Opt_rob
    rob_list[[n]]$rep<-rep_rob
    
  }
  save(rob_list,file=paste0("rob_fit_set_1_fixed_",sub,".RData"))
}

nrep<-10
rob_list<-list()
for (sub in 1:20){
  load(paste0("sub_data_list_set_1_",sub,"RData"))
  for (n in 1:nrep){
    fit_data_rob<-sub_data_list[[2]][[n]]
    fit_data_rob$strat_d<-rep(c(rep(0,16),rep(1,22),rep(2,32)),big_n_t)
    fit_data_rob$robcode<-2
    fit_data_rob$tc<-c(50,8,4)
    
    par_list_rob<-list(theta=logitp(0.5),
                       vess_eff_t=rep(0,length(unique(fit_data_rob$vess_id))),
                       vess_eff_nt=rep(0,length(unique(fit_data_rob$vess_id))),
                       pre_pars_rands=array(data=c(-10,logitp((0.99+1)/2),-1,-8,logitp((0.99+1)/2),-1),dim=c(3,2)),
                       pre_pars_omegas=array(data=c(log(0.1),log(5),0,log(0.1),log(5),0),dim=c(3,2)),
                       log_sigma_vess_nt=-0.5,
                       log_sigma_vess_t=-0.5,
                       rands=matrix(rep(0,big_n_t*2),ncol=2),
                       pg_re=array(rep(0,(length(fit_data_rob$pg_edges)+2)*big_n_t*2),dim=c((length(fit_data_rob$pg_edges)+2),big_n_t,2)),
                       tg_re=matrix(rep(0,length(fit_data_rob$tg_edges)*2),ncol=2))
    random_rob<-c("vess_eff_t","vess_eff_nt","rands","pg_re","tg_re")
    
    map<-list(pre_pars_omegas=as.factor(array(c(1,2,NA,3,4,NA),dim=c(3,2))),
              pre_pars_rands=as.factor(array(c(5,NA,6,7,NA,8),dim=c(3,2))))
    
    obj_rob<-MakeADFun(fit_data_rob,par_list_rob,random=random_rob,DLL="nngp_updated_robust",map=map,silent=F)
    Opt_rob<-optimx::optimr(obj_rob$par,obj_rob$fn,obj_rob$gr,control=list(maxit=100000),method="nlminb")
    while (Opt_rob$message=="iteration limit reached without convergence (10)"){
      obj_rob$par<-obj_rob$env$last.par.best[-which(names(obj_rob$env$last.par.best)%in% random_rob)]
      Opt_rob<-optimx::optimr(obj_rob$par,obj_rob$fn,obj_rob$gr,control=list(maxit=100000),method="nlminb")
    }
    rep_rob<-sdreport(obj_rob)
    
    rob_list[[n]]<-list()
    rob_list[[n]]$obj<-obj_rob
    rob_list[[n]]$Opt<-Opt_rob
    rob_list[[n]]$rep<-rep_rob
    
  }
  save(rob_list,file=paste0("rob_fit_set_1_strat_",sub,".RData"))
}

