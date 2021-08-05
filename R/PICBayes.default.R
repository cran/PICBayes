
PICBayes.default <-
function(
L,
R,
y,
xcov,
IC,
model,
scale.designX,
scaled,
xtrt,
zcov,
area,
binary,
I,
C,
nn,
order=3,
knots,
grids,
a_eta=1,
b_eta=1,
a_ga=1,
b_ga=1,
a_lamb=1,
b_lamb=1,
a_tau=1,
b_tau=1,
a_tau_trt=1,
b_tau_trt=1,
a_alpha=1,
b_alpha=1,
H=5,
a_tau_star=1,
b_tau_star=1,
a_alpha_trt=1,
b_alpha_trt=1,
H_trt=5,
a_tau_trt_star=1,
b_tau_trt_star=1,
beta_iter=1001,
phi_iter=1001,
beta_cand,
phi_cand,
beta_sig0=10,
x_user=NULL,
total=6000,
burnin=1000,
thin=1,
conf.int=0.95,
seed=1,
...){

   if (model=='PIC'){
   est<-PIC(L,R,y,xcov,IC,scale.designX,scaled,binary,order,knots,grids,
   a_eta,b_eta,a_ga,b_ga,beta_iter,beta_cand,beta_sig0,x_user,
   total,burnin,thin,conf.int,seed)
   } else
   if (model=='spatialPIC'){
   est<-spatialPIC(L,R,y,xcov,IC,scale.designX,scaled,area,binary,I,
   C,nn,order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_lamb,b_lamb,beta_iter,
   phi_iter,beta_cand,beta_sig0,x_user,total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterPIC_int'){
   est<-clusterPIC_int(L,R,y,xcov,IC,scale.designX,scaled,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_tau,b_tau,beta_iter,phi_iter,
   beta_cand,phi_cand,beta_sig0,x_user,total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterPIC_int_DP'){
   est<-clusterPIC_int_DP(L,R,y,xcov,IC,scale.designX,scaled,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_alpha,b_alpha,H,a_tau_star,
   b_tau_star,beta_iter,phi_iter,beta_cand,phi_cand,beta_sig0,
   x_user,total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterPIC_trt'){
   est<-clusterPIC_trt(L,R,y,xcov,IC,scale.designX,scaled,xtrt,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_tau,b_tau,a_tau_trt,b_tau_trt,
   beta_iter,phi_iter,beta_cand,phi_cand,beta_sig0,x_user,
   total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterPIC_trt_DP'){
   est<-clusterPIC_trt_DP(L,R,y,xcov,IC,scale.designX,scaled,xtrt,area,binary,
   I,order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_alpha,b_alpha,H,a_tau_star,
   b_tau_star,a_alpha_trt,b_alpha_trt,H_trt,a_tau_trt_star,b_tau_trt_star,
   beta_iter,phi_iter,beta_cand,phi_cand,beta_sig0,x_user,
   total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterPIC_Z'){
   est<-clusterPIC_Z(L,R,y,xcov,IC,scale.designX,scaled,zcov,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_tau,b_tau,beta_iter,phi_iter,
   beta_cand,phi_cand,beta_sig0,x_user,total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterPIC_Z_DP'){
   est<-clusterPIC_Z_DP(L,R,y,xcov,IC,scale.designX,scaled,zcov,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_alpha,b_alpha,H,a_tau_star,
   b_tau_star,beta_iter,phi_iter,beta_cand,phi_cand,beta_sig0,
   x_user,total,burnin,thin,conf.int,seed)
   } else
   if (model=='IC'){
   est<-IC(L,R,y,xcov,IC,scale.designX,scaled,binary,order,knots,grids,
   a_eta,b_eta,a_ga,b_ga,beta_iter,beta_cand,beta_sig0,x_user,
   total,burnin,thin,conf.int,seed)
   } else
   if (model=='spatialIC'){
   est<-spatialIC(L,R,y,xcov,IC,scale.designX,scaled,area,binary,I,C,nn,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_lamb,b_lamb,beta_iter,phi_iter,
   beta_cand,beta_sig0,x_user,total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterIC_int'){
   est<-clusterIC_int(L,R,y,xcov,IC,scale.designX,scaled,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_tau,b_tau,beta_iter,phi_iter,
   beta_cand,phi_cand,beta_sig0,x_user,total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterIC_int_DP'){
   est<-clusterIC_int_DP(L,R,y,xcov,IC,scale.designX,scaled,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_alpha,b_alpha,H,a_tau_star,
   b_tau_star,beta_iter,phi_iter,beta_cand,phi_cand,beta_sig0,
   x_user,total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterIC_trt'){
   est<-clusterIC_trt(L,R,y,xcov,IC,xtrt,scale.designX,scaled,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_tau,b_tau,a_tau_trt,b_tau_trt,
   beta_iter,phi_iter,beta_cand,phi_cand,beta_sig0,x_user,
   total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterIC_trt_DP'){
   est<-clusterIC_trt_DP(L,R,y,xcov,IC,scale.designX,scaled,xtrt,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_alpha,b_alpha,H,a_tau_star,
   b_tau_star,a_alpha_trt,b_alpha_trt,H_trt,a_tau_trt_star,b_tau_trt_star,
   beta_iter,phi_iter,beta_cand,phi_cand,beta_sig0,x_user,
   total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterIC_Z'){
   est<-clusterIC_Z(L,R,y,xcov,IC,scale.designX,scaled,zcov,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_tau,b_tau,beta_iter,phi_iter,
   beta_cand,phi_cand,beta_sig0,x_user,total,burnin,thin,conf.int,seed)
   } else
   if (model=='clusterIC_Z_DP'){
   est<-clusterIC_Z_DP(L,R,y,xcov,IC,scale.designX,scaled,zcov,area,binary,I,
   order,knots,grids,a_eta,b_eta,a_ga,b_ga,a_alpha,b_alpha,H,a_tau_star,
   b_tau_star,beta_iter,phi_iter,beta_cand,phi_cand,beta_sig0,
   x_user,total,burnin,thin,conf.int,seed)
   }

  output<-est

  output$call<-match.call()
  class(output)<-"PICBayes"
  output
}

