
library(dplyr)
library(sf)
library(sp)

library(TMB)
library(starve)
compile("nngp_updated_model.cpp")#,"&> C:/Users/mcdonaldra/Documents/errors.txt")
dyn.load(dynlib("nngp_updated_model"))
compile("nngp_updated_subsim.cpp")
dyn.load(dynlib("nngp_updated_subsim"))

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

#Setting up
n_vess<-40
n_t<-2
strat_n_t<-5
big_n_t<-20
n_k<-3
n_k2<-4

nrep<-10
for (sub in 1:20){
  strat_list<-list()
  fit_sim_stuff<-list()
  rand_list<-list()
  fixed_prop_list<-list()
  fixed_hal_alloc_list<-list()
  strat_prop_list<-list()
  strat_hal_alloc_list<-list()
  for (n in 1:nrep){
    n_sets<-70
    
    while(TRUE){
      props_alloc_sets<-round(c(5/(5+7+10),7/(5+7+10),10/(5+7+10))*n_sets)
      
      #Simulating data:
      n_strat<-n_t*length(sf_grid)
      n_strats_sep<-rep(length(sf_grid),n_t)
      n_fixed<-0
      #Setting up H_i and H_j
      H_j<-rep(NA,n_strat+n_fixed)
      H_i<-rep(NA,n_strat)
      
      main_n2<-rep(NA,n_t)
      for (t in 1:n_t){
        main_n2[t]<-n_strats_sep[t]
      }
      
      main_n<-rep(NA,n_t)
      for (t in 1:n_t){
        main_n[t]<-(n_fixed/n_t)+main_n2[t]
      }
      
      soaky<-rep(1,sum(main_n))
      
      big_A<-matrix(rep(0.3,sum(main_n)*n_k),nrow=sum(main_n),ncol=n_k)
      big_H<-matrix(rep(0.3,sum(main_n2)*n_k2),nrow=sum(main_n2),ncol=n_k2)
      
      H_j<-rep(1,n_strat+n_fixed)
      H_i<-which(H_j==1)-1
      
      vess_id<-sample(c(1:n_vess),sum(main_n),replace=T)
      vess_id<-vess_id-1
      
      sub_starv<-data.frame(obs=big_A[,1],
                            greg=big_A[,2],
                            sim_obs=rpois(nrow(big_A),2),
                            sim_obs2=rpois(nrow(big_A),5),
                            year=rep(1:n_t,each=length(sf_grid)),
                            x=st_coordinates(sf_centroids)[,1],
                            y=st_coordinates(sf_centroids)[,2]) %>% st_as_sf(coords=c("x","y"))
      sub_starv$year<-sub_starv$year-min(sub_starv$year)
      
      pre_set_nodes<-sample(1:nrow(sub_starv),50,replace=F)
      
      strv_stuff<-strv_prepare(cbind(obs,greg)~time(year,type="independent")+space("matern",nu=1),
                               data=sub_starv,n_neighbours = 3,
                               nodes=sub_starv[pre_set_nodes,4])
      
      data<-list(A=big_A,H=big_H,H_i=H_i,H_j=H_j,
                 vess_id=vess_id,s=soaky,
                 n_t=n_t,n=main_n,n_k=3,n2=main_n2,n_k2=4,
                 pg_edges=starve:::convert_idxR_to_C(strv_stuff@process@persistent_graph)@edges,
                 pg_dists=strv_stuff@process@persistent_graph@distances,
                 tg_t=strv_stuff@process@tg_re@locations$year,
                 tg_edges=starve:::convert_idxR_to_C(strv_stuff@process@transient_graph)@edges,
                 tg_dists=strv_stuff@process@transient_graph@distances,
                 cv_code=c(2,2),
                 idx=array(data=c(strv_stuff@observations@data_predictions@locations$graph_idx,
                                  strv_stuff@observations@data_predictions@locations$year),
                           dim=c(length(big_A[,1]),2)),
                 hook_num=c(700,300))
      
      data$idx[,1]<-data$idx[,1]-1
      
      par_list<-list()
      
      par_list$theta<-0.83
      
      par_list$vess_eff_t<-rep(0,length(unique(vess_id)))
      par_list$vess_eff_nt<-rep(0,length(unique(vess_id)))
      
      par_list$pre_pars_rands<-array(data=c(-8,logitp((0.99+1)/2),log(0.3),-6,logitp((0.99+1)/2),log(0.3)),dim=c(3,2))
      par_list$pre_pars_omegas<-array(data=c(log(0.01),log(20),0,log(0.02),log(10),0),dim=c(3,2))
      
      par_list$log_sigma_vess_nt<--1.5
      par_list$log_sigma_vess_t<--1
      
      par_list$rands<-strv_stuff@process@time_effects$w
      par_list$pg_re<-strv_stuff@process@pg_re$w
      par_list$tg_re<-strv_stuff@process@tg_re@values$w
      
      
      map<-list(pre_pars_omegas=as.factor(array(c(1,2,NA,3,4,NA),dim=c(3,2))),
                pre_pars_rands=as.factor(array(c(5,NA,6,7,NA,8),dim=c(3,2))))
      
      random<-c("vess_eff_t","vess_eff_nt","rands","pg_re","tg_re")
      
      strat_sim_data<-list(
        A=matrix(nrow=length(sf_grid)*strat_n_t,ncol=n_k),
        H=matrix(nrow=length(sf_grid)*strat_n_t,ncol=n_k2),
        H_j=rep(1,length(sf_grid)*strat_n_t),
        H_i=c(1:(length(sf_grid)*strat_n_t))-1,
        vess_id=rep(1,length(sf_grid)*strat_n_t),
        s=rep(1,length(sf_grid)*strat_n_t),
        n_t=n_t,n=main_n,n_k=3,n2=main_n2,n_k2=4,
        pg_edges=starve:::convert_idxR_to_C(strv_stuff@process@persistent_graph)@edges,
        pg_dists=strv_stuff@process@persistent_graph@distances,
        tg_t=strv_stuff@process@tg_re@locations$year,
        tg_edges=starve:::convert_idxR_to_C(strv_stuff@process@transient_graph)@edges,
        tg_dists=strv_stuff@process@transient_graph@distances,
        cv_code=c(2,2),
        idx=array(data=c(strv_stuff@observations@data_predictions@locations$graph_idx,
                         strv_stuff@observations@data_predictions@locations$year),
                  dim=c(length(sf_grid)*strat_n_t,2)),
        hook_num=c(700,300)
      )
      
      strat_sim_re<-list(
        vess_eff_t=rep(0,n_vess),
        vess_eff_nt=rep(0,n_vess),
        rands=matrix(nrow=strat_n_t,ncol=2),
        pg_re=array(dim=c(dim(strv_stuff@process@pg_re$w)[1],strat_n_t,dim(strv_stuff@process@pg_re$w)[3])),
        tg_re=matrix(nrow=nrow(strv_stuff@process@tg_re@values$w)*strat_n_t,ncol=2),
        ldat=rep(0,strat_n_t*length(sf_grid)),
        ldant=rep(0,strat_n_t*length(sf_grid))
      )
      
      #First 2 years
      sim_obj<-MakeADFun(data,par_list,random=random,DLL="nngp_updated_model",map=map,silent=F)
      prev_sim<-sim_obj$simulate(complete=T)
      
      strat_sim_data$A[1:(length(sf_grid)*2),]<-prev_sim$A
      strat_sim_data$H[1:(length(sf_grid)*2),]<-prev_sim$H
      strat_sim_data$s[1:(length(sf_grid)*2)]<-prev_sim$s
      strat_sim_re$vess_eff_t<-prev_sim$vess_eff_t
      strat_sim_re$vess_eff_nt<-prev_sim$vess_eff_nt
      strat_sim_re$rands[1:2,]<-prev_sim$rands
      strat_sim_re$pg_re[,1:2,]<-prev_sim$pg_re[,1:2,]
      strat_sim_re$tg_re[1:(length(sf_grid)*2-2*(dim(strat_sim_re$pg_re)[1])),]<-prev_sim$tg_re
      strat_sim_data$idx[1:(length(sf_grid)*2),]<-data$idx
      strat_sim_re$ldat[1:(length(sf_grid)*2)]<-prev_sim$ldat
      strat_sim_re$ldant[1:(length(sf_grid)*2)]<-prev_sim$ldant
      
      for (t in 3:strat_n_t){
        new_soaky<-c(prev_sim$s[(main_n[1]+1):sum(main_n)],rep(1,main_n[2]))
        new_big_A<-rbind(prev_sim$A[(main_n[1]+1):sum(main_n),],
                         matrix(rep(0.3,sum(main_n)/2*n_k),nrow=sum(main_n)/2,ncol=n_k))
        new_big_H<-rbind(prev_sim$H[(main_n[1]+1):sum(main_n),],
                         matrix(rep(0.3,sum(main_n2)/2*n_k2),nrow=sum(main_n2)/2,ncol=n_k2))
        new_vess_id<-c(vess_id[(main_n[1]+1):sum(main_n)],sample(c(1:n_vess),sum(main_n)/2,replace=T)-1)
        
        new_data<-list(A=new_big_A,H=new_big_H,H_i=H_i,H_j=H_j,
                       vess_id=new_vess_id,s=new_soaky,
                       n_t=n_t,n=main_n,n_k=3,n2=main_n2,n_k2=4,
                       pg_edges=starve:::convert_idxR_to_C(strv_stuff@process@persistent_graph)@edges,
                       pg_dists=strv_stuff@process@persistent_graph@distances,
                       tg_t=strv_stuff@process@tg_re@locations$year,
                       tg_edges=starve:::convert_idxR_to_C(strv_stuff@process@transient_graph)@edges,
                       tg_dists=strv_stuff@process@transient_graph@distances,
                       cv_code=c(2,2),
                       idx=array(data=c(strv_stuff@observations@data_predictions@locations$graph_idx,
                                        strv_stuff@observations@data_predictions@locations$year),
                                 dim=c(length(big_A[,1]),2)),
                       hook_num=c(700,300))
        new_data$idx[,1]<-new_data$idx[,1]-1
        
        new_par_list<-par_list
        new_par_list$vess_eff_t<-prev_sim$vess_eff_t
        new_par_list$vess_eff_nt<-prev_sim$vess_eff_nt
        
        new_par_list$rands<-rbind(prev_sim$rands[2,],strv_stuff@process@time_effects$w[2,])
        new_par_list$pg_re[,1,]<-prev_sim$pg_re[,2,]
        new_par_list$tg_re<-rbind(prev_sim$tg_re[(main_n[1]-dim(par_list$pg_re)[1]+1):(sum(main_n)-2*dim(par_list$pg_re)[1]),],strv_stuff@process@tg_re@values$w[1:(main_n[1]-dim(par_list$pg_re)[1]),])
        
        new_obj<-MakeADFun(new_data,new_par_list,random=random,DLL="nngp_updated_subsim",map=map,silent=F)
        prev_sim<-new_obj$simulate(complete=T)
        
        strat_sim_data$A[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),]<-prev_sim$A[(length(sf_grid)+1):(2*length(sf_grid)),]
        strat_sim_data$H[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),]<-prev_sim$H[(length(sf_grid)+1):(2*length(sf_grid)),]
        strat_sim_data$s[((t-1)*length(sf_grid)+1):(t*length(sf_grid))]<-prev_sim$s[(length(sf_grid)+1):(2*length(sf_grid))]
        strat_sim_re$rands[t,]<-prev_sim$rands[2,]
        strat_sim_re$pg_re[,t,]<-prev_sim$pg_re[,2,]
        strat_sim_re$tg_re[((t-1)*length(sf_grid)+1-dim(par_list$pg_re)[1]):((t)*length(sf_grid)-2*dim(par_list$pg_re)[1]),]<-prev_sim$tg_re[(length(sf_grid)+1-dim(par_list$pg_re)[1]):((2)*length(sf_grid)-2*dim(par_list$pg_re)[1]),]
        strat_sim_data$idx[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),]<-new_data$idx[(length(sf_grid)+1):(2*length(sf_grid)),]
        strat_sim_re$ldat[((t-1)*length(sf_grid)+1):(t*length(sf_grid))]<-prev_sim$ldat[(length(sf_grid)+1):(2*length(sf_grid))]
        strat_sim_re$ldant[((t-1)*length(sf_grid)+1):(t*length(sf_grid))]<-prev_sim$ldant[(length(sf_grid)+1):(2*length(sf_grid))]
        
      }
      strat_list[[n]]<-list(strat_sim_data=strat_sim_data,
                            strat_sim_re=strat_sim_re)
      
      #Stratify first one as an example
      sf_grid3<-rbind(sf_grid2[strat_sim_data$idx[1:(main_n[1]),1]+1,],sf_grid2[strat_sim_data$idx[1:(main_n[1]),1]+1,],
                      sf_grid2[strat_sim_data$idx[1:(main_n[1]),1]+1,],sf_grid2[strat_sim_data$idx[1:(main_n[1]),1]+1,])
      plot_strats_frame<-cbind(sf_grid3[,1],strat_list[[n]]$strat_sim_data$A[-(1:main_n[1]),1]+strat_list[[n]]$strat_sim_data$H[-(1:main_n[1]),2],
                               strat_list[[n]]$strat_sim_data$s[-(1:main_n[1])],year=rep(1:4,each=length(sf_grid)))
      colnames(plot_strats_frame)<-c("Catch","Soak_time","Year","geometry")
      st_geometry(plot_strats_frame)<-"geometry"
      plot_strats_frame$catch_rate<-plot_strats_frame$Catch/plot_strats_frame$Soak_time/1000
      
      # strat_plot<-ggplot()+geom_sf(data=plot_strats_frame,aes(fill=catch_rate))+
      #   facet_wrap(~Year)+scale_fill_viridis_c()
      
      mean_vec<-data.frame(plot_strats_frame$catch_rate[1:length(sf_grid)],
                           plot_strats_frame$catch_rate[(length(sf_grid)+1):(2*length(sf_grid))],
                           plot_strats_frame$catch_rate[(2*length(sf_grid)+1):(3*length(sf_grid))],
                           plot_strats_frame$catch_rate[(3*length(sf_grid)+1):(4*length(sf_grid))]) %>% rowMeans()
      sf_mean_frame<-data.frame(geometry=sf_grid2$x,means=mean_vec) %>% st_as_sf()
      
      #To get a spatially-weighed average, only take neighboring cells, so distance =1
      dist_mat<-as.matrix(stats::dist(st_coordinates(sf_centroids)))
      #Multiply by the ones that are distance 1, then everything else =0
      #Add em up, divide by number of 1s
      spat_moving_avg<-rep(0,length(mean_vec))
      for (i in 1:length(spat_moving_avg)){
        indices<-c(i,which(dist_mat[i,]<2 & dist_mat[i,]>0))
        n_mean<-length(indices)
        spat_moving_avg[i]<-sum(mean_vec[indices])/n_mean
      }
      
      sf_spat_avg<-data.frame(geometry=sf_grid2$x,means=spat_moving_avg) %>% st_as_sf()
      
      #They are ordered when put in integer, 1 is lowest
      sf_spat_avg$bins<-as.integer(cut(sf_spat_avg$means,breaks=4))
      sf_spat_avg$bins[which(sf_spat_avg$bins==4)]<-3
      sf_spat_avg$id<-1:nrow(sf_spat_avg)
      
      props_sets<-round((aggregate(means~bins,data=sf_spat_avg,FUN=length)$means/nrow(sf_spat_avg))*n_sets)
      check_rows<-c(nrow(sf_spat_avg[which(sf_spat_avg$bins==1),]),
                    nrow(sf_spat_avg[which(sf_spat_avg$bins==2),]),
                    nrow(sf_spat_avg[which(sf_spat_avg$bins==3),]))
      if(!any(props_sets<5)) break
      if(!any(check_rows>(props_alloc_sets*3))) break
    }
    #Now, we have our 3 strata, so can select locations
    #Let's set the # of data points at 70 to reflect cusk, and simulate 20 years
    n_t_forward<-20
    n_sets<-70
    
    #Random sampling;
    rand_sets<-c()
    for (i in 1:n_t_forward){
      rand_sets<-c(rand_sets,sample(sf_spat_avg$id,n_sets,replace=F))
    }
    
    #Randomly selected fixed stations:
    fixed_sets<-rep(sample(sf_spat_avg$id,n_sets,replace=F),n_t_forward)
    
    #Fixed stations with proportional allocation
    if (sum(props_sets)>n_sets) {
      rem_ind<-sample(1:3,1)
      props_sets[rem_ind]<-props_sets[rem_ind]-1
    } else if (sum(props_sets)<n_sets){
      add_ind<-sample(1:3,1)
      props_sets[add_ind]<-props_sets[add_ind]+1
    }
    fixed_prop_sets<-c()
    for (i in 1:3){
      fixed_prop_sets<-c(fixed_prop_sets,sample(subset(sf_spat_avg,bins==i)$id,props_sets[i],replace=F))
    }
    fixed_prop_sets<-rep(fixed_prop_sets,n_t_forward)
    
    #Fixed Stations following part of halibut design:
    #It is 5:7:10 for low catch:med catch:high catch
    fixed_hal_alloc_sets<-c()
    for (i in 1:3){
      fixed_hal_alloc_sets<-c(fixed_hal_alloc_sets,sample(subset(sf_spat_avg,bins==i)$id,props_alloc_sets[i],replace=F))
    }
    fixed_hal_alloc_sets<-rep(fixed_hal_alloc_sets,n_t_forward)
    
    #Stratified random design
    #Proportional allocations
    strat_prop_sets<-c()
    for (t in 1:n_t_forward){
      for (i in 1:3){
        strat_prop_sets<-c(strat_prop_sets,sample(subset(sf_spat_avg,bins==i)$id,props_sets[i],replace=F))
      }
    }
    
    #Same allocation as halibut one
    strat_hal_alloc_sets<-c()
    for (t in 1:n_t_forward){
      for (i in 1:3){
        strat_hal_alloc_sets<-c(strat_hal_alloc_sets,sample(subset(sf_spat_avg,bins==i)$id,props_alloc_sets[i],replace=F))
      }
    }
    
    fit_sim_data<-list(
      A=matrix(nrow=length(sf_grid)*n_t_forward,ncol=n_k),
      H=matrix(nrow=length(sf_grid)*n_t_forward,ncol=n_k2),
      H_j=rep(1,length(sf_grid)*n_t_forward),
      H_i=c(1:(length(sf_grid)*n_t_forward))-1,
      vess_id=rep(1,length(sf_grid)*n_t_forward),
      s=rep(1,length(sf_grid)*n_t_forward),
      n_t=n_t,n=main_n,n_k=3,n2=main_n2,n_k2=4,
      pg_edges=starve:::convert_idxR_to_C(strv_stuff@process@persistent_graph)@edges,
      pg_dists=strv_stuff@process@persistent_graph@distances,
      tg_t=strv_stuff@process@tg_re@locations$year,
      tg_edges=starve:::convert_idxR_to_C(strv_stuff@process@transient_graph)@edges,
      tg_dists=strv_stuff@process@transient_graph@distances,
      cv_code=c(2,2),
      idx=array(data=c(strv_stuff@observations@data_predictions@locations$graph_idx,
                       strv_stuff@observations@data_predictions@locations$year),
                dim=c(length(sf_grid)*n_t_forward,2)),
      hook_num=c(700,300)
    )
    
    fit_sim_re<-list(
      vess_eff_t=rep(0,n_vess),
      vess_eff_nt=rep(0,n_vess),
      rands=matrix(nrow=n_t_forward,ncol=2),
      pg_re=array(dim=c(dim(strv_stuff@process@pg_re$w)[1],n_t_forward,dim(strv_stuff@process@pg_re$w)[3])),
      tg_re=matrix(nrow=nrow(strv_stuff@process@tg_re@values$w)*n_t_forward,ncol=2),
      ldat=rep(0,n_t_forward*length(sf_grid)),
      ldant=rep(0,n_t_forward*length(sf_grid))
    )
    
    for (t in 1:n_t_forward){
      new_soaky<-c(prev_sim$s[(main_n[1]+1):sum(main_n)],rep(1,main_n[2]))
      new_big_A<-rbind(prev_sim$A[(main_n[1]+1):sum(main_n),],
                       matrix(rep(0.3,sum(main_n)/2*n_k),nrow=sum(main_n)/2,ncol=n_k))
      new_big_H<-rbind(prev_sim$H[(main_n[1]+1):sum(main_n),],
                       matrix(rep(0.3,sum(main_n2)/2*n_k2),nrow=sum(main_n2)/2,ncol=n_k2))
      new_vess_id<-c(vess_id[(main_n[1]+1):sum(main_n)],sample(c(1:n_vess),sum(main_n)/2,replace=T)-1)
      
      new_data<-list(A=new_big_A,H=new_big_H,H_i=H_i,H_j=H_j,
                     vess_id=new_vess_id,s=new_soaky,
                     n_t=n_t,n=main_n,n_k=3,n2=main_n2,n_k2=4,
                     pg_edges=starve:::convert_idxR_to_C(strv_stuff@process@persistent_graph)@edges,
                     pg_dists=strv_stuff@process@persistent_graph@distances,
                     tg_t=strv_stuff@process@tg_re@locations$year,
                     tg_edges=starve:::convert_idxR_to_C(strv_stuff@process@transient_graph)@edges,
                     tg_dists=strv_stuff@process@transient_graph@distances,
                     cv_code=c(2,2),
                     idx=array(data=c(strv_stuff@observations@data_predictions@locations$graph_idx,
                                      strv_stuff@observations@data_predictions@locations$year),
                               dim=c(length(big_A[,1]),2)),
                     hook_num=c(700,300))
      new_data$idx[,1]<-new_data$idx[,1]-1
      
      new_par_list<-par_list
      new_par_list$vess_eff_t<-prev_sim$vess_eff_t
      new_par_list$vess_eff_nt<-prev_sim$vess_eff_nt
      
      new_par_list$rands<-rbind(prev_sim$rands[2,],strv_stuff@process@time_effects$w[2,])
      new_par_list$pg_re[,1,]<-prev_sim$pg_re[,2,]
      new_par_list$tg_re<-rbind(prev_sim$tg_re[(main_n[1]-dim(par_list$pg_re)[1]+1):(sum(main_n)-2*dim(par_list$pg_re)[1]),],strv_stuff@process@tg_re@values$w[1:(main_n[1]-dim(par_list$pg_re)[1]),])
      
      new_obj<-MakeADFun(new_data,new_par_list,random=random,DLL="nngp_updated_subsim",map=map,silent=F)
      prev_sim<-new_obj$simulate(complete=T)
      
      fit_sim_data$A[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),]<-prev_sim$A[(length(sf_grid)+1):(2*length(sf_grid)),]
      fit_sim_data$H[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),]<-prev_sim$H[(length(sf_grid)+1):(2*length(sf_grid)),]
      fit_sim_data$s[((t-1)*length(sf_grid)+1):(t*length(sf_grid))]<-prev_sim$s[(length(sf_grid)+1):(2*length(sf_grid))]
      fit_sim_data$vess_id[((t-1)*length(sf_grid)+1):(t*length(sf_grid))]<-prev_sim$vess_id[(length(sf_grid)+1):(2*length(sf_grid))]
      fit_sim_re$rands[t,]<-prev_sim$rands[2,]
      fit_sim_re$pg_re[,t,]<-prev_sim$pg_re[,2,]
      if (t == 1) {
        fit_sim_re$tg_re[1:((t)*length(sf_grid)-dim(par_list$pg_re)[1]),]<-prev_sim$tg_re[(length(sf_grid)+1-dim(par_list$pg_re)[1]):((2)*length(sf_grid)-2*dim(par_list$pg_re)[1]),]
      } else fit_sim_re$tg_re[((t-1)*length(sf_grid)+1-dim(par_list$pg_re)[1]):((t)*length(sf_grid)-2*dim(par_list$pg_re)[1]),]<-prev_sim$tg_re[(length(sf_grid)+1-dim(par_list$pg_re)[1]):((2)*length(sf_grid)-2*dim(par_list$pg_re)[1]),]
      fit_sim_data$idx[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),]<-new_data$idx[(length(sf_grid)+1):(2*length(sf_grid)),]
      fit_sim_re$ldat[((t-1)*length(sf_grid)+1):(t*length(sf_grid))]<-prev_sim$ldat[(length(sf_grid)+1):(2*length(sf_grid))]
      fit_sim_re$ldant[((t-1)*length(sf_grid)+1):(t*length(sf_grid))]<-prev_sim$ldant[(length(sf_grid)+1):(2*length(sf_grid))]
    }
    
    fit_sim_stuff[[n]]<-list(data=fit_sim_data,re=fit_sim_re,rand_sets=rand_sets,
                             fixed_prop_sets=fixed_prop_sets,
                             fixed_hal_alloc_sets=fixed_hal_alloc_sets,
                             strat_prop_sets=strat_prop_sets,
                             strat_hal_alloc_sets=strat_hal_alloc_sets)
    
    #Can fit model now
    vess_id_rand<-c()
    vess_id_fixed_prop<-c()
    vess_id_fixed_hal_alloc<-c()
    vess_id_strat_prop<-c()
    vess_id_strat_hal_alloc<-c()
    s_rand<-c()
    s_fixed_prop<-c()
    s_fixed_hal_alloc<-c()
    s_strat_prop<-c()
    s_strat_hal_alloc<-c()
    A_rand<-matrix(nrow=n_sets*n_t_forward,ncol=n_k)
    A_fixed_prop<-matrix(nrow=n_sets*n_t_forward,ncol=n_k)
    A_fixed_hal_alloc<-matrix(nrow=n_sets*n_t_forward,ncol=n_k)
    A_strat_prop<-matrix(nrow=n_sets*n_t_forward,ncol=n_k)
    A_strat_hal_alloc<-matrix(nrow=n_sets*n_t_forward,ncol=n_k)
    H_rand<-matrix(nrow=n_sets*n_t_forward,ncol=n_k2)
    H_fixed_prop<-matrix(nrow=n_sets*n_t_forward,ncol=n_k2)
    H_fixed_hal_alloc<-matrix(nrow=n_sets*n_t_forward,ncol=n_k2)
    H_strat_prop<-matrix(nrow=n_sets*n_t_forward,ncol=n_k2)
    H_strat_hal_alloc<-matrix(nrow=n_sets*n_t_forward,ncol=n_k2)
    sub_vess<-data.frame(vess=fit_sim_data$vess_id[1:length(sf_grid)],
                         id=fit_sim_data$idx[1:length(sf_grid),1]+1)
    sub_s<-data.frame(s=fit_sim_data$s[1:length(sf_grid)],
                      id=fit_sim_data$idx[1:length(sf_grid),1]+1)
    sub_A<-data.frame(A1=fit_sim_data$A[1:length(sf_grid),1],
                      A2=fit_sim_data$A[1:length(sf_grid),2],
                      A3=fit_sim_data$A[1:length(sf_grid),3],
                      id=fit_sim_data$idx[1:length(sf_grid),1]+1)
    sub_H<-data.frame(H1=fit_sim_data$H[1:length(sf_grid),1],
                      H2=fit_sim_data$H[1:length(sf_grid),2],
                      H3=fit_sim_data$H[1:length(sf_grid),3],
                      H4=fit_sim_data$H[1:length(sf_grid),4],
                      id=fit_sim_data$idx[1:length(sf_grid),1]+1)
    for (j in 1:n_sets){
      vess_id_rand[j]<-sub_vess$vess[which(sub_vess$id==rand_sets[j])]
      vess_id_fixed_prop[j]<-sub_vess$vess[which(sub_vess$id==fixed_prop_sets[j])]
      vess_id_fixed_hal_alloc[j]<-sub_vess$vess[which(sub_vess$id==fixed_hal_alloc_sets[j])]
      vess_id_strat_prop[j]<-sub_vess$vess[which(sub_vess$id==strat_prop_sets[j])]
      vess_id_strat_hal_alloc[j]<-sub_vess$vess[which(sub_vess$id==strat_hal_alloc_sets[j])]
      s_rand[j]<-sub_s$s[which(sub_s$id==rand_sets[j])]
      s_fixed_prop[j]<-sub_s$s[which(sub_s$id==fixed_prop_sets[j])]
      s_fixed_hal_alloc[j]<-sub_s$s[which(sub_s$id==fixed_hal_alloc_sets[j])]
      s_strat_prop[j]<-sub_s$s[which(sub_s$id==strat_prop_sets[j])]
      s_strat_hal_alloc[j]<-sub_s$s[which(sub_s$id==strat_hal_alloc_sets[j])]
      A_rand[j,]<-as.matrix(sub_A[which(sub_A$id==rand_sets[j]),1:3])
      A_fixed_prop[j,]<-as.matrix(sub_A[which(sub_A$id==fixed_prop_sets[j]),1:3])
      A_fixed_hal_alloc[j,]<-as.matrix(sub_A[which(sub_A$id==fixed_hal_alloc_sets[j]),1:3])
      A_strat_prop[j,]<-as.matrix(sub_A[which(sub_A$id==strat_prop_sets[j]),1:3])
      A_strat_hal_alloc[j,]<-as.matrix(sub_A[which(sub_A$id==strat_hal_alloc_sets[j]),1:3])
      H_rand[j,]<-as.matrix(sub_H[which(sub_H$id==rand_sets[j]),1:4])
      H_fixed_prop[j,]<-as.matrix(sub_H[which(sub_H$id==fixed_prop_sets[j]),1:4])
      H_fixed_hal_alloc[j,]<-as.matrix(sub_H[which(sub_H$id==fixed_hal_alloc_sets[j]),1:4])
      H_strat_prop[j,]<-as.matrix(sub_H[which(sub_H$id==strat_prop_sets[j]),1:4])
      H_strat_hal_alloc[j,]<-as.matrix(sub_H[which(sub_H$id==strat_hal_alloc_sets[j]),1:4])
    }
    for (t in 2:n_t_forward){
      sub_vess<-data.frame(vess=fit_sim_data$vess_id[((t-1)*length(sf_grid)+1):(t*length(sf_grid))],
                           id=fit_sim_data$idx[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),1]+1)
      sub_s<-data.frame(s=fit_sim_data$s[((t-1)*length(sf_grid)+1):(t*length(sf_grid))],
                        id=fit_sim_data$idx[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),1]+1)
      sub_A<-data.frame(A1=fit_sim_data$A[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),1],
                        A2=fit_sim_data$A[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),2],
                        A3=fit_sim_data$A[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),3],
                        id=fit_sim_data$idx[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),1]+1)
      sub_H<-data.frame(H1=fit_sim_data$H[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),1],
                        H2=fit_sim_data$H[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),2],
                        H3=fit_sim_data$H[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),3],
                        H4=fit_sim_data$H[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),4],
                        id=fit_sim_data$idx[((t-1)*length(sf_grid)+1):(t*length(sf_grid)),1]+1)
      for (j in 1:n_sets){
        vess_id_rand[(t-1)*n_sets+j]<-sub_vess$vess[which(sub_vess$id==rand_sets[j])]
        vess_id_fixed_prop[(t-1)*n_sets+j]<-sub_vess$vess[which(sub_vess$id==fixed_prop_sets[j])]
        vess_id_fixed_hal_alloc[(t-1)*n_sets+j]<-sub_vess$vess[which(sub_vess$id==fixed_hal_alloc_sets[j])]
        vess_id_strat_prop[(t-1)*n_sets+j]<-sub_vess$vess[which(sub_vess$id==strat_prop_sets[j])]
        vess_id_strat_hal_alloc[(t-1)*n_sets+j]<-sub_vess$vess[which(sub_vess$id==strat_hal_alloc_sets[j])]
        s_rand[(t-1)*n_sets+j]<-sub_s$s[which(sub_s$id==rand_sets[j])]
        s_fixed_prop[(t-1)*n_sets+j]<-sub_s$s[which(sub_s$id==fixed_prop_sets[j])]
        s_fixed_hal_alloc[(t-1)*n_sets+j]<-sub_s$s[which(sub_s$id==fixed_hal_alloc_sets[j])]
        s_strat_prop[(t-1)*n_sets+j]<-sub_s$s[which(sub_s$id==strat_prop_sets[j])]
        s_strat_hal_alloc[(t-1)*n_sets+j]<-sub_s$s[which(sub_s$id==strat_hal_alloc_sets[j])]
        A_rand[(t-1)*n_sets+j,]<-as.matrix(sub_A[which(sub_A$id==rand_sets[j]),1:3])
        A_fixed_prop[(t-1)*n_sets+j,]<-as.matrix(sub_A[which(sub_A$id==fixed_prop_sets[j]),1:3])
        A_fixed_hal_alloc[(t-1)*n_sets+j,]<-as.matrix(sub_A[which(sub_A$id==fixed_hal_alloc_sets[j]),1:3])
        A_strat_prop[(t-1)*n_sets+j,]<-as.matrix(sub_A[which(sub_A$id==strat_prop_sets[j]),1:3])
        A_strat_hal_alloc[(t-1)*n_sets+j,]<-as.matrix(sub_A[which(sub_A$id==strat_hal_alloc_sets[j]),1:3])
        H_rand[(t-1)*n_sets+j,]<-as.matrix(sub_H[which(sub_H$id==rand_sets[j]),1:4])
        H_fixed_prop[(t-1)*n_sets+j,]<-as.matrix(sub_H[which(sub_H$id==fixed_prop_sets[j]),1:4])
        H_fixed_hal_alloc[(t-1)*n_sets+j,]<-as.matrix(sub_H[which(sub_H$id==fixed_hal_alloc_sets[j]),1:4])
        H_strat_prop[(t-1)*n_sets+j,]<-as.matrix(sub_H[which(sub_H$id==strat_prop_sets[j]),1:4])
        H_strat_hal_alloc[(t-1)*n_sets+j,]<-as.matrix(sub_H[which(sub_H$id==strat_hal_alloc_sets[j]),1:4])
      }
    }
    
    locs_rand<-st_coordinates(sf_centroids[rand_sets,])
    locs_fixed_prop<-st_coordinates(sf_centroids[fixed_prop_sets,])
    locs_fixed_hal_alloc<-st_coordinates(sf_centroids[fixed_hal_alloc_sets,])
    locs_strat_prop<-st_coordinates(sf_centroids[strat_prop_sets,])
    locs_strat_hal_alloc<-st_coordinates(sf_centroids[strat_hal_alloc_sets,])
    
    starve_rand<-data.frame(obs=rep(1,n_sets*n_t_forward),
                            greg=rep(1,n_sets*n_t_forward),
                            year=rep(1:n_t_forward,each=n_sets)-1,
                            x=locs_rand[,1],y=locs_rand[,2]) %>% 
      st_as_sf(coords=c("x","y"))
    starve_fixed_prop<-data.frame(obs=rep(1,n_sets*n_t_forward),
                                  greg=rep(1,n_sets*n_t_forward),
                                  year=rep(1:n_t_forward,each=n_sets)-1,
                                  x=locs_fixed_prop[,1],y=locs_fixed_prop[,2]) %>% 
      st_as_sf(coords=c("x","y"))
    starve_fixed_hal_alloc<-data.frame(obs=rep(1,n_sets*n_t_forward),
                                       greg=rep(1,n_sets*n_t_forward),
                                       year=rep(1:n_t_forward,each=n_sets)-1,
                                       x=locs_fixed_hal_alloc[,1],y=locs_fixed_hal_alloc[,2]) %>% 
      st_as_sf(coords=c("x","y"))
    starve_strat_prop<-data.frame(obs=rep(1,n_sets*n_t_forward),
                                  greg=rep(1,n_sets*n_t_forward),
                                  year=rep(1:n_t_forward,each=n_sets)-1,
                                  x=locs_strat_prop[,1],y=locs_strat_prop[,2]) %>% 
      st_as_sf(coords=c("x","y"))
    starve_strat_hal_alloc<-data.frame(obs=rep(1,n_sets*n_t_forward),
                                       greg=rep(1,n_sets*n_t_forward),
                                       year=rep(1:n_t_forward,each=n_sets)-1,
                                       x=locs_strat_hal_alloc[,1],y=locs_strat_hal_alloc[,2]) %>% 
      st_as_sf(coords=c("x","y"))
    
    strv_rand<-strv_prepare(cbind(obs,greg)~time(year,type="independent")+space("matern",nu=1),
                            data=starve_rand,n_neighbours = 3,
                            nodes=starve_rand[sample(1:nrow(starve_rand),50,replace=F),4])
    strv_fixed_prop<-strv_prepare(cbind(obs,greg)~time(year,type="independent")+space("matern",nu=1),
                                  data=starve_fixed_prop,n_neighbours = 3,
                                  nodes=starve_fixed_prop[sample(1:nrow(starve_fixed_prop),50,replace=F),4])
    strv_fixed_hal_alloc<-strv_prepare(cbind(obs,greg)~time(year,type="independent")+space("matern",nu=1),
                                       data=starve_fixed_hal_alloc,n_neighbours = 3,
                                       nodes=starve_fixed_hal_alloc[sample(1:nrow(starve_fixed_hal_alloc),50,replace=F),4])
    strv_strat_prop<-strv_prepare(cbind(obs,greg)~time(year,type="independent")+space("matern",nu=1),
                                  data=starve_strat_prop,n_neighbours = 3,
                                  nodes=starve_strat_prop[sample(1:nrow(starve_strat_prop),50,replace=F),4])
    strv_strat_hal_alloc<-strv_prepare(cbind(obs,greg)~time(year,type="independent")+space("matern",nu=1),
                                       data=starve_strat_hal_alloc,n_neighbours = 3,
                                       nodes=starve_strat_hal_alloc[sample(1:nrow(starve_strat_hal_alloc),50,replace=F),4])
    
    #Gotta change the data, both Hs, s, and I think thats it?
    
    fit_data_rand<-list(
      A=A_rand,
      H=H_rand,
      H_j=rep(1,n_sets*n_t_forward),
      H_i=c(1:(n_sets*n_t_forward))-1,
      vess_id=vess_id_rand,
      s=s_rand,
      n_t=n_t_forward,n=rep(n_sets,n_t_forward),n_k=3,n2=rep(n_sets,n_t_forward),n_k2=4,
      pg_edges=starve:::convert_idxR_to_C(strv_rand@process@persistent_graph)@edges,
      pg_dists=strv_rand@process@persistent_graph@distances,
      tg_t=strv_rand@process@tg_re@locations$year,
      tg_edges=starve:::convert_idxR_to_C(strv_rand@process@transient_graph)@edges,
      tg_dists=strv_rand@process@transient_graph@distances,
      cv_code=c(2,2),
      idx=array(data=c(strv_rand@observations@data_predictions@locations$graph_idx,
                       strv_rand@observations@data_predictions@locations$year),
                dim=c(n_sets*n_t_forward,2)),
      hook_num=c(700,300)
    )
    fit_data_rand$idx[,1]<-fit_data_rand$idx[,1]-1
    
    par_list_rand<-list(theta=logitp(0.5),
                        vess_eff_t=rep(0,length(unique(fit_data_rand$vess_id))),
                        vess_eff_nt=rep(0,length(unique(fit_data_rand$vess_id))),
                        pre_pars_rands=array(data=c(-10,logitp((0.99+1)/2),-1,-8,logitp((0.99+1)/2),-1),dim=c(3,2)),
                        pre_pars_omegas=array(data=c(log(0.1),log(5),0,log(0.1),log(5),0),dim=c(3,2)),
                        log_sigma_vess_nt=-0.5,
                        log_sigma_vess_t=-0.5,
                        rands=strv_rand@process@time_effects$w,
                        pg_re=strv_rand@process@pg_re$w,
                        tg_re=strv_rand@process@tg_re@values$w)
    random<-c("vess_eff_t","vess_eff_nt","rands","pg_re","tg_re")
    
    map<-list(pre_pars_omegas=as.factor(array(c(1,2,NA,3,4,NA),dim=c(3,2))),
              pre_pars_rands=as.factor(array(c(5,NA,6,7,NA,8),dim=c(3,2))))
    
    obj_rand<-MakeADFun(fit_data_rand,par_list_rand,random=random,DLL="nngp_updated_model",map=map,silent=F)
    Opt_rand<-optimx::optimr(obj_rand$par,obj_rand$fn,obj_rand$gr,control=list(maxit=100000),method="nlminb")
    rep_rand<-sdreport(obj_rand)
    
    rand_list[[n]]<-list()
    rand_list[[n]]$obj<-obj_rand
    rand_list[[n]]$Opt<-Opt_rand
    rand_list[[n]]$rep<-rep_rand
    
    
    fit_data_fixed_prop<-list(
      A=A_fixed_prop,
      H=H_fixed_prop,
      H_j=rep(1,n_sets*n_t_forward),
      H_i=c(1:(n_sets*n_t_forward))-1,
      vess_id=vess_id_fixed_prop,
      s=s_fixed_prop,
      n_t=n_t_forward,n=rep(n_sets,n_t_forward),n_k=3,n2=rep(n_sets,n_t_forward),n_k2=4,
      pg_edges=starve:::convert_idxR_to_C(strv_fixed_prop@process@persistent_graph)@edges,
      pg_dists=strv_fixed_prop@process@persistent_graph@distances,
      tg_t=strv_fixed_prop@process@tg_re@locations$year,
      tg_edges=starve:::convert_idxR_to_C(strv_fixed_prop@process@transient_graph)@edges,
      tg_dists=strv_fixed_prop@process@transient_graph@distances,
      cv_code=c(2,2),
      idx=array(data=c(strv_fixed_prop@observations@data_predictions@locations$graph_idx,
                       strv_fixed_prop@observations@data_predictions@locations$year),
                dim=c(n_sets*n_t_forward,2)),
      hook_num=c(700,300)
    )
    fit_data_fixed_prop$idx[,1]<-fit_data_fixed_prop$idx[,1]-1
    
    par_list_fixed_prop<-list(theta=logitp(0.5),
                              vess_eff_t=rep(0,length(unique(fit_data_fixed_prop$vess_id))),
                              vess_eff_nt=rep(0,length(unique(fit_data_fixed_prop$vess_id))),
                              pre_pars_rands=array(data=c(-10,logitp((0.99+1)/2),-1,-8,logitp((0.99+1)/2),-1),dim=c(3,2)),
                              pre_pars_omegas=array(data=c(log(0.1),log(5),0,log(0.1),log(5),0),dim=c(3,2)),
                              log_sigma_vess_nt=-0.5,
                              log_sigma_vess_t=-0.5,
                              rands=strv_fixed_prop@process@time_effects$w,
                              pg_re=strv_fixed_prop@process@pg_re$w,
                              tg_re=strv_fixed_prop@process@tg_re@values$w)
    
    obj_fixed_prop<-MakeADFun(fit_data_fixed_prop,par_list_fixed_prop,random=random,DLL="nngp_updated_model",map=map,silent=F)
    Opt_fixed_prop<-optimx::optimr(obj_fixed_prop$par,obj_fixed_prop$fn,obj_fixed_prop$gr,control=list(maxit=100000),method="nlminb")
    rep_fixed_prop<-sdreport(obj_fixed_prop)
    
    fixed_prop_list[[n]]<-list()
    fixed_prop_list[[n]]$obj<-obj_fixed_prop
    fixed_prop_list[[n]]$Opt<-Opt_fixed_prop
    fixed_prop_list[[n]]$rep<-rep_fixed_prop
    
    fit_data_fixed_hal_alloc<-list(
      A=A_fixed_hal_alloc,
      H=H_fixed_hal_alloc,
      H_j=rep(1,n_sets*n_t_forward),
      H_i=c(1:(n_sets*n_t_forward))-1,
      vess_id=vess_id_fixed_hal_alloc,
      s=s_fixed_hal_alloc,
      n_t=n_t_forward,n=rep(n_sets,n_t_forward),n_k=3,n2=rep(n_sets,n_t_forward),n_k2=4,
      pg_edges=starve:::convert_idxR_to_C(strv_fixed_hal_alloc@process@persistent_graph)@edges,
      pg_dists=strv_fixed_hal_alloc@process@persistent_graph@distances,
      tg_t=strv_fixed_hal_alloc@process@tg_re@locations$year,
      tg_edges=starve:::convert_idxR_to_C(strv_fixed_hal_alloc@process@transient_graph)@edges,
      tg_dists=strv_fixed_hal_alloc@process@transient_graph@distances,
      cv_code=c(2,2),
      idx=array(data=c(strv_fixed_hal_alloc@observations@data_predictions@locations$graph_idx,
                       strv_fixed_hal_alloc@observations@data_predictions@locations$year),
                dim=c(n_sets*n_t_forward,2)),
      hook_num=c(700,300)
    )
    fit_data_fixed_hal_alloc$idx[,1]<-fit_data_fixed_hal_alloc$idx[,1]-1
    
    par_list_fixed_hal_alloc<-list(theta=logitp(0.5),
                                   vess_eff_t=rep(0,length(unique(fit_data_fixed_hal_alloc$vess_id))),
                                   vess_eff_nt=rep(0,length(unique(fit_data_fixed_hal_alloc$vess_id))),
                                   pre_pars_rands=array(data=c(-10,logitp((0.99+1)/2),-1,-8,logitp((0.99+1)/2),-1),dim=c(3,2)),
                                   pre_pars_omegas=array(data=c(log(0.1),log(5),0,log(0.1),log(5),0),dim=c(3,2)),
                                   log_sigma_vess_nt=-0.5,
                                   log_sigma_vess_t=-0.5,
                                   rands=strv_fixed_hal_alloc@process@time_effects$w,
                                   pg_re=strv_fixed_hal_alloc@process@pg_re$w,
                                   tg_re=strv_fixed_hal_alloc@process@tg_re@values$w)
    
    obj_fixed_hal_alloc<-MakeADFun(fit_data_fixed_hal_alloc,par_list_fixed_hal_alloc,random=random,DLL="nngp_updated_model",map=map,silent=F)
    Opt_fixed_hal_alloc<-optimx::optimr(obj_fixed_hal_alloc$par,obj_fixed_hal_alloc$fn,obj_fixed_hal_alloc$gr,control=list(maxit=100000),method="nlminb")
    rep_fixed_hal_alloc<-sdreport(obj_fixed_hal_alloc)
    
    fixed_hal_alloc_list[[n]]<-list()
    fixed_hal_alloc_list[[n]]$obj<-obj_fixed_hal_alloc
    fixed_hal_alloc_list[[n]]$Opt<-Opt_fixed_hal_alloc
    fixed_hal_alloc_list[[n]]$rep<-rep_fixed_hal_alloc
    
    fit_data_strat_prop<-list(
      A=A_strat_prop,
      H=H_strat_prop,
      H_j=rep(1,n_sets*n_t_forward),
      H_i=c(1:(n_sets*n_t_forward))-1,
      vess_id=vess_id_strat_prop,
      s=s_strat_prop,
      n_t=n_t_forward,n=rep(n_sets,n_t_forward),n_k=3,n2=rep(n_sets,n_t_forward),n_k2=4,
      pg_edges=starve:::convert_idxR_to_C(strv_strat_prop@process@persistent_graph)@edges,
      pg_dists=strv_strat_prop@process@persistent_graph@distances,
      tg_t=strv_strat_prop@process@tg_re@locations$year,
      tg_edges=starve:::convert_idxR_to_C(strv_strat_prop@process@transient_graph)@edges,
      tg_dists=strv_strat_prop@process@transient_graph@distances,
      cv_code=c(2,2),
      idx=array(data=c(strv_strat_prop@observations@data_predictions@locations$graph_idx,
                       strv_strat_prop@observations@data_predictions@locations$year),
                dim=c(n_sets*n_t_forward,2)),
      hook_num=c(700,300)
    )
    fit_data_strat_prop$idx[,1]<-fit_data_strat_prop$idx[,1]-1
    
    par_list_strat_prop<-list(theta=logitp(0.5),
                              vess_eff_t=rep(0,length(unique(fit_data_strat_prop$vess_id))),
                              vess_eff_nt=rep(0,length(unique(fit_data_strat_prop$vess_id))),
                              pre_pars_rands=array(data=c(-10,logitp((0.99+1)/2),-1,-8,logitp((0.99+1)/2),-1),dim=c(3,2)),
                              pre_pars_omegas=array(data=c(log(0.1),log(5),0,log(0.1),log(5),0),dim=c(3,2)),
                              log_sigma_vess_nt=-0.5,
                              log_sigma_vess_t=-0.5,
                              rands=strv_strat_prop@process@time_effects$w,
                              pg_re=strv_strat_prop@process@pg_re$w,
                              tg_re=strv_strat_prop@process@tg_re@values$w)
    
    obj_strat_prop<-MakeADFun(fit_data_strat_prop,par_list_strat_prop,random=random,DLL="nngp_updated_model",map=map,silent=F)
    Opt_strat_prop<-optimx::optimr(obj_strat_prop$par,obj_strat_prop$fn,obj_strat_prop$gr,control=list(maxit=100000),method="nlminb")
    rep_strat_prop<-sdreport(obj_strat_prop)
    
    strat_prop_list[[n]]<-list()
    strat_prop_list[[n]]$obj<-obj_strat_prop
    strat_prop_list[[n]]$Opt<-Opt_strat_prop
    strat_prop_list[[n]]$rep<-rep_strat_prop
    
    fit_data_strat_hal_alloc<-list(
      A=A_strat_hal_alloc,
      H=H_strat_hal_alloc,
      H_j=rep(1,n_sets*n_t_forward),
      H_i=c(1:(n_sets*n_t_forward))-1,
      vess_id=vess_id_strat_hal_alloc,
      s=s_strat_hal_alloc,
      n_t=n_t_forward,n=rep(n_sets,n_t_forward),n_k=3,n2=rep(n_sets,n_t_forward),n_k2=4,
      pg_edges=starve:::convert_idxR_to_C(strv_strat_hal_alloc@process@persistent_graph)@edges,
      pg_dists=strv_strat_hal_alloc@process@persistent_graph@distances,
      tg_t=strv_strat_hal_alloc@process@tg_re@locations$year,
      tg_edges=starve:::convert_idxR_to_C(strv_strat_hal_alloc@process@transient_graph)@edges,
      tg_dists=strv_strat_hal_alloc@process@transient_graph@distances,
      cv_code=c(2,2),
      idx=array(data=c(strv_strat_hal_alloc@observations@data_predictions@locations$graph_idx,
                       strv_strat_hal_alloc@observations@data_predictions@locations$year),
                dim=c(n_sets*n_t_forward,2)),
      hook_num=c(700,300)
    )
    fit_data_strat_hal_alloc$idx[,1]<-fit_data_strat_hal_alloc$idx[,1]-1
    
    par_list_strat_hal_alloc<-list(theta=logitp(0.5),
                                   vess_eff_t=rep(0,length(unique(fit_data_strat_hal_alloc$vess_id))),
                                   vess_eff_nt=rep(0,length(unique(fit_data_strat_hal_alloc$vess_id))),
                                   pre_pars_rands=array(data=c(-10,logitp((0.99+1)/2),-1,-8,logitp((0.99+1)/2),-1),dim=c(3,2)),
                                   pre_pars_omegas=array(data=c(log(0.1),log(5),0,log(0.1),log(5),0),dim=c(3,2)),
                                   log_sigma_vess_nt=-0.5,
                                   log_sigma_vess_t=-0.5,
                                   rands=strv_strat_hal_alloc@process@time_effects$w,
                                   pg_re=strv_strat_hal_alloc@process@pg_re$w,
                                   tg_re=strv_strat_hal_alloc@process@tg_re@values$w)
    
    obj_strat_hal_alloc<-MakeADFun(fit_data_strat_hal_alloc,par_list_strat_hal_alloc,random=random,DLL="nngp_updated_model",map=map,silent=F)
    Opt_strat_hal_alloc<-optimx::optimr(obj_strat_hal_alloc$par,obj_strat_hal_alloc$fn,obj_strat_hal_alloc$gr,control=list(maxit=100000),method="nlminb")
    rep_strat_hal_alloc<-sdreport(obj_strat_hal_alloc)
    
    strat_hal_alloc_list[[n]]<-list()
    strat_hal_alloc_list[[n]]$obj<-obj_strat_hal_alloc
    strat_hal_alloc_list[[n]]$Opt<-Opt_strat_hal_alloc
    strat_hal_alloc_list[[n]]$rep<-rep_strat_hal_alloc
    
    sub_data_list<-list()
    sub_data_list[[1]]<-fit_data_fixed_hal_alloc
    sub_data_list[[2]]<-fit_data_strat_hal_alloc
  }
  save(strat_list,file=paste0("strat_list_set_1_",sub,".RData"))
  save(sub_data_list,file=paste0("sub_data_list_set_1_",sub,"RData"))
  save(fit_sim_stuff,file=paste0("fit_sim_stuff_set_1_",sub,".RData"))
  save(rand_list,file=paste0("rand_model_fit_set_1_",sub,".RData"))
  save(fixed_prop_list,file=paste0("fixed_prop_model_fit_set_1_",sub,".RData"))
  save(fixed_hal_alloc_list,file=paste0("fixed_hal_alloc_model_fit_set_1_",sub,".RData"))
  save(strat_prop_list,file=paste0("strat_prop_model_fit_set_1_",sub,".RData"))
  save(strat_hal_alloc_list,file=paste0("strat_hal_alloc_model_fit_set_1_",sub,".RData"))
}



#Now setting 2, will have to very suddenly drop the mean levels or vary where they are
