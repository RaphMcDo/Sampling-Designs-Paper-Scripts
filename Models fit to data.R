library(TMB)
library(sf)
library(starve)
library(optimx)

compile("nngp_updated_model.cpp")
dyn.load(dynlib("nngp_updated_model"))

load("preset_data_ST_MEM.RData")

sys1<-Sys.time()
obj<-MakeADFun(data,par_list,random=random,DLL="nngp_updated_model",map=map,silent=F)
sys2<-Sys.time()
Opt<-optimx::optimr(obj$par,obj$fn,obj$gr,control=list(maxit=100000),method="nlminb")
rep<-sdreport(obj)
Report<-obj$report()
sys3<-Sys.time()

compile("nngp_updated_robust.cpp")
dyn.load("nngp_updated_robust")

load("preset_data_rob.RData")

sys1<-Sys.time()
obj_rob<-MakeADFun(data,par_list,random=random,DLL="nngp_updated_robust",map=map,silent=F)
sys2<-Sys.time()
Opt_rob<-optimx::optimr(obj_rob$par,obj_rob$fn,obj_rob$gr,control=list(maxit=100000),method="nlminb")
rep_rob<-sdreport(obj_rob)
Report_rob<-obj_rob$report()
sys3<-Sys.time()

models<-c("ST-MEM","Rob-ST-MEM")

vess_effect_t<-data.frame(val=rep(NA,length(unique(vess_id))*2),
                          sd=rep(NA,length(unique(vess_id))*2),
                          model=rep(models,each=length(unique(vess_id))))
vess_effect_nt<-data.frame(val=rep(NA,length(unique(vess_id))*2),
                           sd=rep(NA,length(unique(vess_id))*2),
                           model=rep(models,each=length(unique(vess_id))))
for (mod in models){
  for (i in 1:length(unique(vess_id))){
    if (mod == models[[1]]){
      vess_effect_t[which(vess_effect_t$model == mod),1:2][i,]<-c(rep$value[which(names(rep$value)=="vess_eff_t")][i],rep$sd[which(names(rep$value)=="vess_eff_t")][i])
      vess_effect_nt[which(vess_effect_nt$model == mod),1:2][i,]<-c(rep$value[which(names(rep$value)=="vess_eff_nt")][i],rep$sd[which(names(rep$value)=="vess_eff_nt")][i]) 
    } else if (mod == models[[2]]){
      vess_effect_t[which(vess_effect_t$model == mod),1:2][i,]<-c(rep_rob$value[which(names(rep_rob$value)=="vess_eff_t")][i],rep_rob$sd[which(names(rep_rob$value)=="vess_eff_t")][i])
      vess_effect_nt[which(vess_effect_nt$model == mod),1:2][i,]<-c(rep_rob$value[which(names(rep_rob$value)=="vess_eff_nt")][i],rep_rob$sd[which(names(rep_rob$value)=="vess_eff_nt")][i])
    } 
  } 
}

rand_t<-data.frame(val=c(rep$value[which(names(rep$value)=="rands")][1:n_t],rep_rob$value[which(names(rep_rob$value)=="rands")][1:n_t]),
                   sd=c(rep$sd[which(names(rep$value)=="rands")][1:n_t],rep_rob$sd[which(names(rep_rob$value)=="rands")][1:n_t]),
                   model=rep(factor(models),each=n_t),
                   spec="Cusk",
                   year=rep(2001:2021,length(models)))
rand_nt<-data.frame(val=c(rep$value[which(names(rep$value)=="rands")][(n_t+1):(2*n_t)],rep_rob$value[which(names(rep_rob$value)=="rands")][(n_t+1):(2*n_t)]),
                    sd=c(rep$sd[which(names(rep$value)=="rands")][(n_t+1):(2*n_t)],rep_rob$sd[which(names(rep_rob$value)=="rands")][(n_t+1):(2*n_t)]),
                    model=rep(factor(models),each=n_t),
                    spec="Non-Target",
                    year=rep(2001:2021,length(models)))
rands_all<-rbind(rand_t,rand_nt)
library(forcats)
rands_all<-rands_all %>% mutate(model=fct_relevel(model,models))

rand_plot<-ggplot(data=rands_all)+
  geom_line(aes(x=year,y=val,col=model))+
  geom_point(aes(x=year,y=val,col=model))+
  geom_ribbon(aes(x=year,ymin=val-1.96*sd,ymax=val+1.96*sd,col=model,fill=model),alpha=0.2)+
  scale_color_viridis_d(name="Model")+
  scale_fill_viridis_d(name="Model")+
  theme_bw()+
  xlab("Year")+ylab("Random Intercept")+
  facet_wrap(~spec)

bid = read.csv("blockIDkey.csv", header = T)
loc_id = data.frame(lon=bid$lon.DecDeg, lat=bid$lat.DecDeg)
sf_loc_id<-st_as_sf(loc_id,coords=c("lon","lat"))
st_crs(sf_loc_id)<-4326
sf_sub_loc_id<-st_intersection(sf_loc_id,divs_4X)
sub_loc_id<-st_coordinates(sf_sub_loc_id)

library(gissr)
library(deldir)
library(spatstat)
library(alphahull)

bound_a = ashape(sub_loc_id, alpha = 0.084)
bound_a_index = bound_a$alpha.extremes
# Boundary points and survey stations
bound_a_pos = sub_loc_id[bound_a_index, ]
bound_a_pos = data.frame(long=bound_a_pos[,1], lat=bound_a_pos[,2])
bound_a_pos$num = 1:length(bound_a_pos$long)
#Get rid of dumb points
bound_a_pos<-bound_a_pos[-c(18,19),]
check_bound<-bound_a_pos[which(bound_a_pos$long>-64.5 & bound_a_pos$long<64.1 & bound_a_pos$lat>44.2 & bound_a_pos$lat<44.5),]

tf_bound_a_pos = SpatialPoints(bound_a_pos,proj4string = prj4s)
tf_bound_a_pos = spTransform(tf_bound_a_pos,utm.prj4s)
tf_bound_a_pos = data.frame(long = tf_bound_a_pos@coords[,1], lat = tf_bound_a_pos@coords[,2])
plot(tf_bound_a_pos)

tf_bound_a_pos2 = sort_points(tf_bound_a_pos, "lat", "long", clockwise = FALSE)
tf_bound_a_pos2 = data.frame(long = tf_bound_a_pos2$long, lat = tf_bound_a_pos2$lat, num=1:length(bound_a_pos$long))

tf_bound_a_pos22<-tf_bound_a_pos2[c(1,2,4,5,6,7,9,10,11,12,13,14,15,17,
                                    18,19,20,21,23,24,25,27,28,29,31,32,
                                    33,35,36,38,39,40,41,42,44,45,46,47,
                                    48,49,50,52,53,55,56,57,59,60,61,63,
                                    64,65,67,68,69,71,72,74,75,76,77,79,
                                    80,81,82,84,85,86,88,89,90,92,93,94,
                                    95,96,98,99,101,102,104,105,106,107,
                                    108,109,110,111,113,114,115,116,117,
                                    118,119,121,123,124,125,127,128,129,
                                    131,132,133,134,136,138,140,141,143,
                                    145,147,148,150,153,155,157,158,165,
                                    162,160,163,166,168,171,176,175,174,
                                    173,172,170,169,167,164,161,159,156,
                                    154,152,151,149,146,144,142,139,137,
                                    135,130,126,122,120,112,103,100,97,
                                    91,87,83,78,73,70,66,62,58,54,51,43,
                                    37,34,30,26,22,16,8,3,279,277,275,
                                    271,269,268,264,262,260,256,248,244,
                                    242,239,236,234,231,228,224,222,219,
                                    217,214,211,209,206,202,198,194,192,
                                    190,188,186,182,183,178,179,180,181,
                                    184,185,187,189,191,193,196,204,203,
                                    201,199,195,197,200,205,207,208,210,
                                    212,213,215,216,218,220,221,223,225,
                                    226,227,229,230,232,233,235,237,238,
                                    240,241,243,245,246,247,249,251,254,
                                    257,255,258,259,250,252,253,261,263,
                                    265,266,267,270,272,273,274,276,278,
                                    280,284,281,285,282,286,283,287,288),]

pp_list<-list()
diri_list<-list()
for (i in c(2001:2021)){
  temp<-subset(sf_rem_fix_data,YEAR==i)
  temp_dist<-st_coordinates(temp)
  
  pp_list[[i-2000]]=ppp(temp_dist[,1], temp_dist[,2], window=owin(poly=list(x=tf_bound_a_pos22$long, y=tf_bound_a_pos22$lat)))
  
  diri_list[[i-2000]] = dirichlet(pp_list[[i-2000]])
}

library(raster)
sf_loc_id<-st_as_sf(as.data.frame(sub_loc_id),coords=c("X","Y"),crs=prj4s) %>%
  st_transform(utm.prj4s)
station_IDs<-list()
for (i in c(1:length(diri_list))){
  temp<-subset(sf_rem_fix_data,YEAR==(i+2000))
  temp_dist<-st_coordinates(temp)
  
  temp_distance<-pointDistance(sf_loc_id,temp_dist,lonlat=F)
  
  station_ID<-c()
  
  for (j in 1:length(st_coordinates(sf_loc_id)[,1])) {
    station_ID<-c(station_ID,which(temp_distance[j,]==min(temp_distance[j,])))
  }
  
  station_IDs[[i]]<-station_ID
  
}

area_list<-list(list(),list(),list(),list())
straight_means<-list()
for (mod in 1:length(models)){
  if (mod==1){
    temp_Report<-Report
    temp_rep<-rep
  } else if (mod==2){
    temp_Report<-Report_rob
    temp_rep<-rep_rob
  } 
  for (i in c(1:length(diri_list))){
    wts = c()
    for (j in 1:length(temp_Report$ldat[Yearly_indices[i,]$ind_1:Yearly_indices[i,]$ind_2])){
      temp_area<-length(which(station_IDs[[i]]==j))*4
      wts[j]<-temp_area
    }
    mean_wts<-rep(1,length(wts))
    # Target
    dweight_ave_t = sum((temp_Report$ldat[Yearly_indices[i,]$ind_1:Yearly_indices[i,]$ind_2])*wts)/sum(wts)
    weighted_t<-log(dweight_ave_t)
    str_mean_t<-sum((temp_Report$ldat[Yearly_indices[i,]$ind_1:Yearly_indices[i,]$ind_2])*mean_wts)/sum(mean_wts)
    # Standard error calculation
    wts = as.matrix(wts)
    mean_wts<-as.matrix(mean_wts)
    # wi/sum(wi)
    wts_2 = wts/(sum(wts))
    mean_wts_2<-mean_wts/sum(mean_wts)
    # covariance matrix for estimated lambda t
    cov = temp_rep$cov
    # cov_t = cov[(Yearly_indices[i,]$ind_1+109):(Yearly_indices[i,]$ind_2+109), (Yearly_indices[i,]$ind_1+109):(Yearly_indices[i,]$ind_2+109)]
    cov_t = cov[which(names(temp_rep$value)=="log_ldat")[Yearly_indices[i,]$ind_1:Yearly_indices[i,]$ind_2],which(names(temp_rep$value)=="log_ldat")[Yearly_indices[i,]$ind_1:Yearly_indices[i,]$ind_2]]
    
    # Standard error for Dirichlet method
    se2_dir_t = t(wts_2)%*%cov_t%*%wts_2
    se_t<-se2_dir_t
    mean_se2_t<-t(mean_wts_2)%*%cov_t%*%mean_wts_2
    # Non-target
    dweight_ave_nt = sum((temp_Report$ldant[Yearly_indices[i,]$ind_1:Yearly_indices[i,]$ind_2])*wts)/sum(wts)
    weighted_nt<-log(dweight_ave_nt)
    str_mean_nt<-sum((temp_Report$ldant[Yearly_indices[i,]$ind_1:Yearly_indices[i,]$ind_2])*mean_wts)/sum(mean_wts)
    # covariance matrix for estimated lambda nt
    # cov_nt = cov[(tail(Yearly_indices$ind_2,n=1)+Yearly_indices[i,]$ind_1+109):(tail(Yearly_indices$ind_2,n=1)+Yearly_indices[i,]$ind_2+109), (tail(Yearly_indices$ind_2,n=1)+Yearly_indices[i,]$ind_1+109) :(tail(Yearly_indices$ind_2,n=1)+Yearly_indices[i,]$ind_2+109)]
    cov_nt = cov[which(names(temp_rep$value)=="log_ldant")[Yearly_indices[i,]$ind_1:Yearly_indices[i,]$ind_2],which(names(temp_rep$value)=="log_ldant")[Yearly_indices[i,]$ind_1:Yearly_indices[i,]$ind_2]]
    # Standard error for Dirichlet method
    se2_dir_nt = t(wts_2)%*%cov_nt%*%wts_2
    se_nt<-se2_dir_nt
    mean_se2_nt<-t(mean_wts_2)%*%cov_nt%*%mean_wts_2
    
    area_list[[mod]][[i]]<-data.frame(weighted_t,sqrt(se_t),weighted_nt,sqrt(se_nt))
    straight_means[[i]]<-data.frame(str_mean_t,sqrt(mean_se2_t),str_mean_nt,sqrt(mean_se2_t))
  }
}

weighted_avgs_frame_nngp<-data.frame(w_t=rep(NA,n_t*2),se_t=rep(NA,n_t*2),
                                     w_nt=rep(NA,n_t*2),se_nt=rep(NA,n_t*2),
                                     model=rep(factor(models),each=n_t),
                                     year=rep(2001:2021,2))
str_means<-data.frame(w_t=rep(NA,n_t),se_t=rep(NA,n_t),
                      w_nt=rep(NA,n_t),se_nt=rep(NA,n_t),
                      year=2001:2021)

for (mod in 1:2){
  for (i in c(1:n_t)){
    weighted_avgs_frame_nngp[which(weighted_avgs_frame_nngp$model==models[mod]),]$w_t[i]<-area_list[[mod]][[i]]$weighted_t
    weighted_avgs_frame_nngp[which(weighted_avgs_frame_nngp$model==models[mod]),]$se_t[i]<-area_list[[mod]][[i]]$sqrt.se_t.
    weighted_avgs_frame_nngp[which(weighted_avgs_frame_nngp$model==models[mod]),]$w_nt[i]<-area_list[[mod]][[i]]$weighted_nt
    weighted_avgs_frame_nngp[which(weighted_avgs_frame_nngp$model==models[mod]),]$se_nt[i]<-area_list[[mod]][[i]]$sqrt.se_nt.
  }
}

all_weighted_avgs<-weighted_avgs_frame_nngp
all_weighted_avgs<-all_weighted_avgs %>% mutate(model=fct_relevel(model,c(models)))

comp_ldat<-ggplot(data=all_weighted_avgs)+
  geom_point(aes(x=year,y=exp(w_t),col=model,shape=model))+
  geom_line(aes(x=year,y=exp(w_t),col=model,lty=model))+
  geom_ribbon(aes(x=year,ymax=exp(w_t+1.96*se_t),ymin=exp(w_t-1.96*se_t),col=model,fill=model),alpha=0.2)+
  theme_bw()+
  ylab("Mean Catch Rate (# of fish per hook per hour)")+
  xlab("Year")+
  scale_color_viridis_d(name="Model")+
  scale_fill_viridis_d(name="Model")+
  scale_shape(name="Model")+
  scale_linetype(name="Model")

comp_ldant<-ggplot(data=all_weighted_avgs)+
  geom_point(aes(x=year,y=exp(w_nt),col=model,shape=model))+
  geom_line(aes(x=year,y=exp(w_nt),col=model,lty=model))+
  geom_ribbon(aes(x=year,ymax=exp(w_nt+1.96*se_nt),ymin=exp(w_nt-1.96*se_nt),col=model,fill=model),alpha=0.2)+
  theme_bw()+
  ylab("Mean Catch Rate (# of fish per hook per hour)")+
  xlab("Year")+
  scale_color_viridis_d(name="Model")+
  scale_fill_viridis_d(name="Model")+
  scale_shape(name="Model")+
  scale_linetype(name="Model")

