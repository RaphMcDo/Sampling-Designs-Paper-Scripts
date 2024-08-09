#include <TMB.hpp> // Links in the TMB libraries

//A lot of these things were taken (with permission) from Ethan Lawler's staRVe package
#include "include_subsim/sim_multinom.hpp"
#include "include_subsim/utils.hpp"
#include "include_subsim/conditional_normal.hpp"
#include "include_subsim/covariance.hpp"
#include "include_subsim/time_series2.hpp"
#include "include_subsim/dag.hpp"
#include "include_subsim/persistent_graph.hpp"
#include "include_subsim/pg_cache.hpp"
#include "include_subsim/transient_graph.hpp"
#include "include_subsim/tg_cache.hpp"
#include "include_subsim/nngp.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  using namespace density;
  
  DATA_MATRIX(H); //Halibut data with fine_scale info on hook condition
  DATA_MATRIX(A); //Halibut data with non-fine scale info on hook condition
  PARAMETER(theta); // Parameters for pnt (probability of caught non-target fish escaping)
  DATA_VECTOR(s); // Soak time
  
  DATA_INTEGER(n_t); //Number of years
  DATA_IVECTOR(n); //For A (vector of number of observations in each year for fixed station dataset)
  DATA_INTEGER(n_k); //For A (number of columns)
  DATA_IVECTOR(n2); //For H (vector of number of observations in each year for stratified dataset)
  DATA_INTEGER(n_k2); //For H (number of columns)
  DATA_IVECTOR(H_i); //Indicator if subset of data present for full MEM to identify which ldat/ldant to be related to
  DATA_IVECTOR(H_j); //Indicator if subset of data present for multinomial part
  DATA_IVECTOR(vess_id);
  
  PARAMETER_VECTOR(vess_eff_t); //Vessel effects for target species
  PARAMETER_VECTOR(vess_eff_nt); //Vessel effects for non-target species
  
  PARAMETER_ARRAY(pre_pars_rands); //NNGP temporal component, all fixed so they don't do anything
  PARAMETER_ARRAY(pre_pars_omegas); //NNGP spatial parameters
  PARAMETER(log_sigma_vess_t); //Variance for vessel effects (target species)
  PARAMETER(log_sigma_vess_nt); //Variance for vessel effects (non-target species)
  
  //Parameter modifications
  array<Type> pars_rands = pre_pars_rands;
  for(int v=0; v<pre_pars_rands.dim(1); v++) {
    pars_rands(0,v) = pre_pars_rands(0,v); // mean, no transformation
    pars_rands(1,v) = 2*invlogit(pre_pars_rands(1,v))-1; // ar1, --> (-1,1)
    pars_rands(2,v) = exp(pre_pars_rands(2,v)); // marginal sd, --> (0,Inf)
  }
  array<Type> pars_omegas = pre_pars_omegas;
  for(int v=0; v<pre_pars_omegas.dim(1); v++) {
    pars_omegas(0,v) = exp(pre_pars_omegas(0,v));
    pars_omegas(1,v) = exp(pre_pars_omegas(1,v));
    pars_omegas(2,v) = exp(pre_pars_omegas(2,v));
    //Matern [sd, range, nu] --> [(0,Inf), (0,Inf), (0,Inf)]
  }
  Type sigma_vess_t = exp(log_sigma_vess_t);
  Type sigma_vess_nt = exp(log_sigma_vess_nt);
  
  //Set up likelihood components
  vector <Type> nll_joint(7); nll_joint.setZero();
  
  PARAMETER_ARRAY(rands); //Random intercept in each year for target and non-target (fixed to zero, don't do anything here)
  
  //Setting up the temporal structure (random intercept)
  time_series<Type> rand_ints {rands, pars_rands};
  nll_joint(0) -= rand_ints.loglikelihood();
  SIMULATE{
    rands = rand_ints.simulate().get_re();
    REPORT(rands);
  }
  
  //Setting up spatio-temporal component
  //  Set up persistent graph
  PARAMETER_ARRAY(pg_re); // [space,time,var]
  DATA_STRUCT(pg_edges,directed_graph);
  DATA_STRUCT(pg_dists,dag_dists);
  dag<Type> pg_g {pg_edges.dag, pg_dists.dag_dist};
  persistent_graph<Type> pg {pg_re, pg_re, pg_g};
  
  // Set up transient graph
  PARAMETER_ARRAY(tg_re); // [idx,var]
  DATA_IVECTOR(tg_t);
  DATA_STRUCT(tg_edges,directed_graph);
  DATA_STRUCT(tg_dists,dag_dists);
  dag<Type> tg_g {tg_edges.dag, tg_dists.dag_dist};
  transient_graph<Type> tg {tg_re, tg_re, tg_g, tg_t, pg.dim_t()};
  
  //Set up covariance
  DATA_IVECTOR(cv_code);
  vector<covariance<Type> > cv(cv_code.size());
  for(int v=0; v<cv_code.size(); v++) {
    cv(v) = covariance<Type> {vector<Type>(pars_omegas.col(v)), cv_code(v)};
  }
  
  //Finish up spatio-temporal component
  nngp<Type> process {pg, tg, cv};
  nll_joint(1) -= process.loglikelihood(rand_ints);
  SIMULATE{
    process.simulate(rand_ints);
    pg_re = process.get_pg_re();
    tg_re = process.get_tg_re();
    REPORT(pg_re);
    REPORT(tg_re);
  }
  
  //Simulate soaking time, values based on observed soak times
  SIMULATE{
    int county = 0;
    for (int t=1; t < n_t; t++){
      if (t>0)  county = county + n(t-1);
      for (int i = 0; i < n(t); i++){
        s((county+i)) = exp(rnorm(log(Type(250)),Type(0.2)))/60;
      }
    }
    REPORT(s);
  }
  
  //Vessel effects
  Type n_vess = vess_eff_t.size();
  for (int i=0; i < n_vess; i++){
    nll_joint(3) -= dnorm(vess_eff_t(i),Type(0.0),sigma_vess_t,true);
    nll_joint(4) -= dnorm(vess_eff_nt(i),Type(0.0),sigma_vess_nt,true);
  }
  
  // SIMULATE{
  //   for (int i=0; i < n_vess; i++){
  //     vess_eff_t(i) = rnorm(Type(0.0),sigma_vess_t);
  //     vess_eff_nt(i) = rnorm(Type(0.0),sigma_vess_nt);
  //   }
  //   REPORT(vess_eff_t);
  //   REPORT(vess_eff_nt);
  // }
  
  DATA_IARRAY(idx); //Necessary to identify proper graph + node + year
  DATA_MATRIX(non_stat) // To make non-stationary temporal mean, since it's all normal can just add/substract
  // ldat - the relative abundance indices for target species (halibut)
  vector <Type> ldat(A.col(0).size());
  vector<Type> ldant(A.col(0).size());
  //These are mostly so I can just extract them more easily to look at the values
  vector <Type> spat_temp_t(A.col(0).size());
  vector<Type> spat_temp_nt(A.col(0).size());
  int counter1 = 0;
  for (int t=0; t < n_t; t++){
    if (t > 0) counter1 = counter1+ n(t-1);
    for(int i=0; i<n(t); ++i){
      ldat((counter1+i))=exp(vess_eff_t(vess_id((counter1+i)))+process(idx((counter1+i),0),idx((counter1+i),1),0)+non_stat(idx((counter1+i),0),idx((counter1+i),1)));
      spat_temp_t((counter1+i))=process(idx((counter1+i),0),idx((counter1+i),1),0);
    }
    
    // ldant - the relative abundance indices for non-target species
    for(int i=0; i<n(t); ++i){
      // Adding omega (random field) and covariates to lambda nt
      ldant((counter1+i))=exp(vess_eff_nt(vess_id((counter1+i)))+process(idx((counter1+i),0),idx((counter1+i),1),1));
      spat_temp_nt((counter1+i))=process(idx((counter1+i),0),idx((counter1+i),1),1);
    }
  }
  REPORT(spat_temp_t);
  REPORT(spat_temp_nt);
  
  // pnt - the probability of caught non-target species escaping
  Type pnt;
  pnt=invlogit(theta);
  ADREPORT(pnt);
  
  //Corresponding probabilities for fixed stations dataset and 700 hooks portion of stratified
  // Reduced MEM
  
  int counter2 = 0;
  matrix<Type> q(A.col(0).size(), n_k);
  for (int t = 0; t < n_t; t++){
    if (t > 0) counter2 = counter2 + n(t-1);
    for(int i=0;i<n(t);i++){
      q((counter2+i),0)=(Type(1)-exp(-(ldat((counter2+i))+ldant((counter2+i)))*s((counter2+i))))*(ldat((counter2+i))/(ldat((counter2+i))+ldant((counter2+i))));
      q((counter2+i),1)=(Type(1)-exp(-(ldat((counter2+i))+ldant((counter2+i)))*s((counter2+i))))*(ldant((counter2+i))/(ldat((counter2+i))+ldant((counter2+i))))*(Type(1)-pnt);
      q((counter2+i),2)=1-q((counter2+i),0)-q((counter2+i),1);
    }
  }
  
  // Corresponding probabilities for N_b, N_t, N_nt and N_e for 300 hooks portion of stratified dataset
  // N_b is the nmuber of baited hooks
  // N_t is the number of target species caught
  // N_nt is the number of non-target species caught
  // N_e is the number of empty hooks
  //Fine-scale hook info (full MEM)
  
  int counter4 = 0;
  matrix<Type> p(H.col(0).size(), n_k2);
  for (int t = 0; t < n_t; t++){
    if (t > 0) counter4 = counter4 + n2(t-1);
    for(int i=0;i<n2(t);i++){
      p((counter4+i),0)=exp(-(ldat(H_i(counter4+i))+ldant(H_i(counter4+i)))*s(H_i(counter4+i)));
      p((counter4+i),1)=(Type(1)-exp(-(ldat(H_i(counter4+i))+ldant(H_i(counter4+i)))*s(H_i(counter4+i))))*(ldat(H_i(counter4+i))/(ldat(H_i(counter4+i))+ldant(H_i(counter4+i))));
      p((counter4+i),2)=(Type(1)-exp(-(ldat(H_i(counter4+i))+ldant(H_i(counter4+i)))*s(H_i(counter4+i))))*(ldant(H_i(counter4+i))/(ldat(H_i(counter4+i))+ldant(H_i(counter4+i))))*(Type(1)-pnt);
      p((counter4+i),3)=1-p((counter4+i),0)-p((counter4+i),1)-p((counter4+i),2);
    }
  }
  
  // The likelihood function for multinomial distribution
  int counter3 = 0;
  int counter_H = 0;
  for (int t = 0; t < n_t; t++){
    if (t > 0) counter3 = counter3+n(t-1);
    for(int i=0; i < n(t); i++){
      if (H_j(counter3+i)==0){
        vector<Type> q_row=q.row((counter3+i));
        vector<Type> A_row=A.row((counter3+i));
        nll_joint(2) -= dmultinom(A_row,q_row,true);
      } else if (H_j(counter3+i)==1){
        vector<Type> q_row=q.row((counter3+i));
        vector<Type> A_row=A.row((counter3+i));
        vector<Type> p_row=p.row((counter_H));
        vector<Type> H_row=H.row((counter_H));
        nll_joint(2) -= (dmultinom(A_row,q_row,true) + dmultinom(H_row,p_row,true));
        counter_H += 1;
      }
    }
  }
  
  SIMULATE{
    DATA_VECTOR(hook_num);
    int counter_sim3 = 0;
    int counter_simH = n(0);
    for (int t = 1; t < n_t; t++){
      if (t > 0) counter_sim3 = counter_sim3+n(t-1);
      for(int i=0; i < n(t); i++){
        if (H_j(counter_sim3+i)==0){
          vector<Type> q_row=q.row((counter_sim3+i));
          vector<Type> A_row=A.row((counter_sim3+i));
          A_row = sim_multinom(hook_num(0)+hook_num(1),q_row);
          A.row((counter_sim3+i))=A_row;
        } else if (H_j(counter_sim3+i)==1){
          vector<Type> q_row=q.row((counter_sim3+i));
          vector<Type> A_row=A.row((counter_sim3+i));
          vector<Type> p_row=p.row((counter_simH));
          vector<Type> H_row=H.row((counter_simH));
          A_row = sim_multinom(hook_num(0),q_row);
          H_row = sim_multinom(hook_num(1),p_row);
          A.row((counter_sim3+i))=A_row;
          H.row(counter_simH)=H_row;
          counter_simH += 1;
        }
      }
    }
    REPORT(H);
    REPORT(A);
  }
  
  REPORT(nll_joint);
  
  // REPORT(D);
  REPORT(ldat);
  REPORT(ldant);
  REPORT(pnt);
  REPORT(q);
  REPORT(p);
  REPORT(rands);
  ADREPORT(rands);
  REPORT(pg_re);
  REPORT(tg_re);
  REPORT(pre_pars_rands);
  ADREPORT(pre_pars_rands);
  REPORT(pre_pars_omegas);
  ADREPORT(pre_pars_omegas);
  ADREPORT(ldat);
  ADREPORT(ldant);
  REPORT(vess_eff_t);
  ADREPORT(vess_eff_t);
  REPORT(vess_eff_nt);
  ADREPORT(vess_eff_nt);
  ADREPORT(sigma_vess_t);
  ADREPORT(sigma_vess_nt);
  
  Type nll = nll_joint.sum();
  return nll;
}
