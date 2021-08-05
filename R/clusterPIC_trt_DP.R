clusterPIC_trt_DP <-
function(
L, 
R, 
y, 
xcov,
IC,
scale.designX,
scaled,
xtrt,
area,
binary, 
I,
order, 
knots,
grids, 
a_eta,
b_eta,
a_ga,
b_ga,
a_alpha,
b_alpha,
H,
a_tau_star,
b_tau_star,
a_alpha_trt,
b_alpha_trt,
H_trt,
a_tau_trt_star,
b_tau_trt_star,
beta_iter,
phi_iter,
beta_cand,
phi_cand,
beta_sig0, 
x_user,
total,
burnin,
thin,
conf.int,
seed){


Ispline<-function(x,order,knots){
# M Spline function with order+1. or I spline with order
# x is a row vector
# knots are a sequence of increasing points

k=order+1
m=length(knots)
n=m-2+k # number of free parameters for M spline family
t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
for (l in k:n){
    yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
}

yytem1=yy1
for (ii in 1:order){
   yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
   for (i in (k-ii):n){
      yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
   }
   yytem1=yytem2
}

index=rep(0,length(x))
for (i in 1:length(x)){
    index[i]=sum(t<=x[i])
}

yy=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

if (order==1){
   for (i in 2:n){
      yy[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
   }
}else{
   for (j in 1:length(x)){
      for (i in 2:n){
         if (i<(index[j]-order+1)){
            yy[i-1,j]=1
         }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
            yy[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
         }else{
            yy[i-1,j]=0
         }
      }
   }
}
return(yy)
}


### get Mspline bases ###
Mspline<-function(x,order,knots){

k1=order
m=length(knots)
n1=m-2+k1 # number of parameters
t1=c(rep(1,k1)*knots[1], knots[2:(m-1)], rep(1,k1)*knots[m]) # new knots

tem1=array(rep(0,(n1+k1-1)*length(x)),dim=c(n1+k1-1, length(x)))
for (l in k1:n1){
    tem1[l,]=(x>=t1[l] & x<t1[l+1])/(t1[l+1]-t1[l])
}

if (order==1){
   mbases=tem1
}else{
   mbases=tem1
   for (ii in 1:(order-1)){
      tem=array(rep(0,(n1+k1-1-ii)*length(x)),dim=c(n1+k1-1-ii, length(x)))
      for (i in (k1-ii):n1){
         tem[i,]=(ii+1)*((x-t1[i])*mbases[i,]+(t1[i+ii+1]-x)*mbases[i+1,])/(t1[i+ii+1]-t1[i])/ii
      }
      mbases=tem
   }
}
return(mbases)
}


poissrndpositive<-function(lambda){
  q=200
  t=seq(0,q,1)
  p=dpois(t,lambda)
  pp=cumsum(p[2:(q+1)])/(1-p[1])
  u=runif(1)

  while(u>pp[q]){
    q=q+1
    pp[q]=pp[q-1]+dpois(q,lambda)/(1-p[1])
  }
  ll=sum(u>pp)+1
}


### main routine ###
   set.seed(seed)
   L=matrix(L,ncol=1)
   R=matrix(R,ncol=1)
   y=matrix(y,ncol=1)
   xcov=as.matrix(xcov)
   xtrt=matrix(xtrt,ncol=1)
   area=matrix(area,ncol=1)
   IC=matrix(IC,ncol=1)
   p=ncol(xcov)

   if (scale.designX==TRUE){
   mean_X<-apply(xcov,2,mean)
   sd_X<-apply(xcov,2,sd)
   for (r in 1:p){ 
   if (scaled[r]==1) xcov[,r]<-(xcov[,r]-mean_X[r])/sd_X[r]
   }
   }

   n1=sum(IC==0)
   n2=sum(IC==1)
   N=n1+n2

   t<-rep(0,N)
   for (i in 1:N) {t[i]=ifelse(IC[i]==0,L[i],0)}

   ## generate basis functions
   K=length(knots)-2+order 
   kgrids=length(grids)
   G<-length(x_user)/p     # number of survival curves

   bmsT=Mspline(t[1:n1],order,knots)    # K*n1
   bisT=Ispline(t[1:n1],order,knots)
   bisL=Ispline(L[(n1+1):N],order,knots) # K*n2
   bisR=Ispline(R[(n1+1):N],order,knots)
   bisg=Ispline(grids,order,knots)

   ## initial value
   eta=rgamma(1,a_eta,rate=b_eta)
   gamcoef=matrix(rgamma(K, 1, rate=1),ncol=K)

   phicoef<-matrix(rep(0,I),ncol=1)
   phicoef2<-matrix(rep(0,N),ncol=1)
   for (j in 1:N){
   phicoef2[j]<-phicoef[area[j]]  #a vector of phi values for all subjects
   }

   phicoef_trt<-matrix(rep(0,I),ncol=1)
   phicoef2_trt<-matrix(rep(NA,N),ncol=1)
   for (j in 1:N){
   phicoef2_trt[j]<-phicoef_trt[area[j]]  #a vector of phi values for all subjects
   }

   beta=matrix(rep(0,p),p,1)
   beta_original=matrix(rep(0,p),ncol=1)

   u=array(rep(0,n1*K),dim=c(n1,K))
   for (i in 1:n1){
   u[i,]=rmultinom(1,1,gamcoef*t(bmsT[,i]))
   }

   alpha=1
   tau_star<-rep(1,H)

   ns<-v<-rep(NA,H)   
   tau<-rep(NA,I)
   LK<-array(rep(NA,I*H),dim=c(I,H))

   pi<-rdirichlet(1,rep(alpha/H,H))   #H x 1
   c<-sample(seq(1:H),I,replace=T,pi)
   for (h in 1:H) ns[h]<-length(c[c==h])  
   for (i in 1:I) tau[i]<-tau_star[c[i]]

   alpha_trt=1
   tau_trt_star<-rep(1,H_trt)

   ns_trt<-v_trt<-rep(NA,H)   
   LK_trt<-array(rep(NA,I*H_trt),dim=c(I,H_trt))
   tau_trt<-rep(NA,I)

   pi_trt<-rdirichlet(1,rep(alpha_trt/H_trt,H_trt))  
   c_trt<-sample(seq(1:H_trt),I,replace=T,pi_trt)
   for (h in 1:H_trt) ns_trt[h]<-length(c_trt[c_trt==h])
   for (i in 1:I) tau_trt[i]<-tau_trt_star[c_trt[i]]

   lambdatT=t(gamcoef%*%bmsT)
   LambdatT=t(gamcoef%*%bisT)  # n1 x 1 # lambda0 at each t_i value
   LambdatL=t(gamcoef%*%bisL) # n2 x 1 # lambda0 at each u_i value
   LambdatR=t(gamcoef%*%bisR) # n2 x 1 # lambda0 at each v_i value
   Lambdatg=t(gamcoef%*%bisg)  # lambda0 at each grids_i value 

   parbeta=array(rep(0,total*p),dim=c(total,p))
   parbeta_original=array(rep(0,total*p),dim=c(total,p))
   pareta=array(rep(0,total),dim=c(total,1)) 
   paralpha=array(rep(0,total),dim=c(total,1)) 
   paralpha_trt=array(rep(0,total),dim=c(total,1)) 
   partau_star=array(rep(0,total*H),dim=c(total,H))
   partau_trt_star=array(rep(0,total*H_trt),dim=c(total,H_trt))
   parpi=array(rep(0,total*H),dim=c(total,H))
   parpi_trt=array(rep(0,total*H_trt),dim=c(total,H_trt))
   parphi=array(rep(0,total*I),dim=c(total,I))
   parphi_trt=array(rep(0,total*I),dim=c(total,I))
   pargam=array(rep(0,total*K),dim=c(total,K))
   parsurv0=array(rep(0,total*kgrids),dim=c(total,kgrids))
   parlambdatT=array(rep(0,total*n1),dim=c(total,n1))
   parLambdatT=array(rep(0,total*n1),dim=c(total,n1))
   parLambdatL=array(rep(0,total*n2),dim=c(total,n2))
   parLambdatR=array(rep(0,total*n2),dim=c(total,n2))
   pardev=array(rep(0,total),dim=c(total,1))
   parfinv_exact=array(rep(0,total*n1),dim=c(total,n1))
   parfinv_IC=array(rep(0,total*n2),dim=c(total,n2))

   if (is.null(x_user)){parsurv=parsurv0} else {
   G<-length(x_user)/p
   parsurv=array(rep(0,total*kgrids*G),dim=c(total,kgrids*G))}
   
   ## iteration
   iter=1
   while (iter<total+1)
   { 
      # sample z, zz, w and ww
      z=array(rep(0,n2),dim=c(n2,1)); w=z
      zz=array(rep(0,n2*K),dim=c(n2,K)); ww=zz

      for (j in 1:n2){
      if (y[n1+j]==0){
      templam1=LambdatR[j]*exp(xcov[(n1+j),]%*%beta+phicoef[area[n1+j]]+phicoef_trt[area[n1+j]]*xtrt[(n1+j)])
      z[j]=poissrndpositive(templam1)
      zz[j,]=rmultinom(1,z[j],gamcoef*t(bisR[,j]))
      } else if (y[n1+j]==1){
      templam1=(LambdatR[j]-LambdatL[j])*exp(xcov[(n1+j),]%*%beta+phicoef[area[n1+j]]+phicoef_trt[area[n1+j]]*xtrt[(n1+j)])
      w[j]=poissrndpositive(templam1)
      ww[j,]=rmultinom(1,w[j],gamcoef*t(bisR[,j]-bisL[,j]))
      }
      }

      # sample beta
      te1=z*(y[(n1+1):N]==0)+w*(y[(n1+1):N]==1)   #n2 x 1
      te2=(LambdatR*(y[(n1+1):N]==0)+LambdatR*(y[(n1+1):N]==1)+LambdatL*(y[(n1+1):N]==2))   #n2 x 1
      te3=LambdatT    #n1 x 1

# sample beta MH-sampler
for (r in 1:p){
if (binary[r]==0){
beta1<-beta2<-beta
if (iter<beta_iter) sd_cand<-beta_cand[r] else sd_cand<-sd(parbeta[1:(iter-1),r])   #specify sd for candidate distribution
xt<-beta[r] #current state
yt<-rnorm(1,xt,sd_cand) #candidate point
beta1[r]<-yt
beta2[r]<-xt
log_f1<-sum(yt*xcov[1:n1,r]-te3*exp(xcov[1:n1,]%*%beta1+phicoef2[1:n1]+phicoef2_trt[1:n1]*xtrt[1:n1]))+sum(yt*xcov[(n1+1):N,r]*te1)-sum(exp(xcov[(n1+1):N,]%*%beta1+phicoef2[(n1+1):N]+phicoef2_trt[(n1+1):N]*xtrt[(n1+1):N])*te2)-0.5*(yt^2)/(beta_sig0^2)
log_f2<-sum(xt*xcov[1:n1,r]-te3*exp(xcov[1:n1,]%*%beta2+phicoef2[1:n1]+phicoef2_trt[1:n1]*xtrt[1:n1]))+sum(xt*xcov[(n1+1):N,r]*te1)-sum(exp(xcov[(n1+1):N,]%*%beta2+phicoef2[(n1+1):N]+phicoef2_trt[(n1+1):N]*xtrt[(n1+1):N])*te2)-0.5*(xt^2)/(beta_sig0^2)
num<-log_f1 #proposal dist is symmetric
den<-log_f2
if (log(runif(1))<(num-den)) beta[r]<-yt else beta[r]<-xt
}

if (binary[r]==1 & p>1){
te4=sum(xcov[1:n1,r])+sum(xcov[(n1+1):N,r]*te1)
te5=sum(te3*exp(as.matrix(xcov[1:n1,-r])%*%as.matrix(beta[-r])+phicoef2[1:n1]+phicoef2_trt[1:n1]*xtrt[1:n1])*xcov[1:n1,r])+sum(te2*exp(as.matrix(xcov[(n1+1):N,-r])%*%as.matrix(beta[-r])+phicoef2[(n1+1):N]+phicoef2_trt[(n1+1):N]*xtrt[(n1+1):N])*xcov[(n1+1):N,r])
beta[r]<-log(rgamma(1,a_ga+te4,rate=b_ga+te5))
}

if (binary[r]==1 & p==1){
te4=sum(xcov[1:n1,r])+sum(xcov[(n1+1):N,r]*te1)
te5=sum(te3*exp(phicoef2[1:n1]+phicoef2_trt[1:n1]*xtrt[1:n1])*xcov[1:n1,r])+sum(te2*exp(phicoef2[(n1+1):N]+phicoef2_trt[(n1+1):N]*xtrt[(n1+1):N])*xcov[(n1+1):N,r])
beta[r]<-log(rgamma(1,a_ga+te4,rate=b_ga+te5))
}
}

# convert to beta_original
if (scale.designX==TRUE){
for (r in 1:p) beta_original[r]<-ifelse(scaled[r]==1,beta[r]/sd_X[r],beta[r])
} 
if (scale.designX==FALSE) beta_original<-beta
#print(beta_original)

      # sample gamcoef
      for (l in 1:K){
         tempa=1+sum(u[,l])+sum(zz[,l]*(y[(n1+1):N]==0)+ww[,l]*(y[(n1+1):N]==1))
         tempb=eta+sum(bisT[l,]*exp(xcov[1:n1,]%*%beta+phicoef2[1:n1]+phicoef2_trt[1:n1]*xtrt[1:n1]))+sum(((bisR[l,])*(y[(n1+1):N]==0)+(bisR[l,])*(y[(n1+1):N]==1)
+(bisL[l,])*(y[(n1+1):N]==2))*exp(xcov[(n1+1):N,]%*%beta+phicoef2[(n1+1):N]+phicoef2_trt[(n1+1):N]*xtrt[(n1+1):N]))
         gamcoef[l]=rgamma(1,tempa,rate=tempb)
      }

      lambdatT=t(gamcoef%*%bmsT)
LambdatT=t(gamcoef%*%bisT) # n1 x 1
      LambdatL=t(gamcoef%*%bisL) # n2 x 1
      LambdatR=t(gamcoef%*%bisR) # n2 x 1

# sample u
      u=array(rep(0,n1*K),dim=c(n1,K))
for (i in 1:n1){
u[i,]=rmultinom(1,1,gamcoef*t(bmsT[,i]))
}

      #sample eta
      eta=rgamma(1,a_eta+K, rate=b_eta+sum(gamcoef))

# sample phi_int MH-sampler
phi1<-array(rep(0,N),dim=c(N,1))
phi2<-array(rep(0,N),dim=c(N,1))
for (i in 1:I){
phi1<-phi2<-phicoef2
if (iter<phi_iter) sd_cand<-phi_cand else sd_cand<-sd(parphi[1:(iter-1),i]) #specify sd for candidate distribution
xt<-phicoef[i] #current state
yt<-rnorm(1,xt,sd_cand) # candidate point
phi1[area==i]<-yt
phi2[area==i]<-xt
log_f1<-sum(phi1[1:n1])-sum(te3*exp(xcov[1:n1,]%*%beta+phi1[1:n1]+phicoef2_trt[1:n1]*xtrt[1:n1]))+sum(phi1[(n1+1):N]*te1)-sum((exp(xcov[(n1+1):N,]%*%beta+phi1[(n1+1):N]+phicoef2_trt[(n1+1):N]*xtrt[(n1+1):N]))*te2)-0.5*tau[i]*(yt^2)
log_f2<-sum(phi2[1:n1])-sum(te3*exp(xcov[1:n1,]%*%beta+phi2[1:n1]+phicoef2_trt[1:n1]*xtrt[1:n1]))+sum(phi2[(n1+1):N]*te1)-sum((exp(xcov[(n1+1):N,]%*%beta+phi2[(n1+1):N]+phicoef2_trt[(n1+1):N]*xtrt[(n1+1):N]))*te2)-0.5*tau[i]*(xt^2)
num<-log_f1 #proposal dist is symmetric
den<-log_f2
if (log(runif(1))<(num-den)) phicoef[i]<-yt else phicoef[i]<-xt
phicoef2[area==i]<-phicoef[i]
}

# sample phi_trt MH-sampler
phi1_trt<-array(rep(0,N),dim=c(N,1))
phi2_trt<-array(rep(0,N),dim=c(N,1))
for (i in 1:I){
phi1_trt<-phi2_trt<-phicoef2_trt
if (iter<phi_iter) sd_cand<-phi_cand else sd_cand<-sd(parphi_trt[1:(iter-1),i]) #specify sd for candidate distribution
xt<-phicoef_trt[i] #current state
yt<-rnorm(1,xt,sd_cand) # candidate point
phi1_trt[area==i]<-yt
phi2_trt[area==i]<-xt
log_f1<-sum(phi1_trt[1:n1]*xtrt[1:n1])-sum(te3*exp(xcov[1:n1,]%*%beta+phicoef2[1:n1]+phi1_trt[1:n1]*xtrt[1:n1]))+sum(phi1_trt[(n1+1):N]*xtrt[(n1+1):N]*te1)-sum(te2*exp(xcov[(n1+1):N,]%*%beta+phicoef2[(n1+1):N]+phi1_trt[(n1+1):N]*xtrt[(n1+1):N]))-0.5*tau_trt[i]*(yt^2)
log_f2<-sum(phi2_trt[1:n1]*xtrt[1:n1])-sum(te3*exp(xcov[1:n1,]%*%beta+phicoef2[1:n1]+phi2_trt[1:n1]*xtrt[1:n1]))+sum(phi2_trt[(n1+1):N]*xtrt[(n1+1):N]*te1)-sum(te2*exp(xcov[(n1+1):N,]%*%beta+phicoef2[(n1+1):N]+phi2_trt[(n1+1):N]*xtrt[(n1+1):N]))-0.5*tau_trt[i]*(xt^2)
num<-log_f1 #proposal dist is symmetric
den<-log_f2
if (log(runif(1))<(num-den)) phicoef_trt[i]<-yt else phicoef_trt[i]<-xt
phicoef2_trt[area==i]<-phicoef_trt[i]
}

# update tau_star
for (h in 1:H){
tau_star[h]<-rgamma(1,a_tau_star+0.5*ns[h],b_tau_star+0.5*t(phicoef[c==h])%*%phicoef[c==h])
}

# update tau_trt_star
for (h in 1:H_trt){
tau_trt_star[h]<-rgamma(1,a_tau_trt_star+0.5*ns_trt[h],rate=b_tau_trt_star+0.5*t(phicoef_trt[c_trt==h])%*%phicoef_trt[c_trt==h])
}

# update c, component ID for each center
for (i in 1:I){
for (h in 1:H){
LK[i,h]<-dnorm(phicoef[i],0,1/sqrt(tau_star[h]))   #annealing parameter to improve clustering 
}
p_c<-(pi*LK[i,])/sum(pi*LK[i,])
c[i]<-sample(seq(1:H),1,replace=T,p_c)
}

# update c_trt, component ID for each center
for (i in 1:I){
for (h in 1:H_trt){
LK_trt[i,h]<-dnorm(phicoef_trt[i],0,1/sqrt(tau_trt_star[h]))   #annealing parameter to improve clustering 
}
p_c_trt<-(pi_trt*LK_trt[i,])/sum(pi_trt*LK_trt[i,])
c_trt[i]<-sample(seq(1:H_trt),1,replace=T,p_c_trt)
}

# update ns
for (h in 1:H) ns[h]<-length(c[c==h])   #must allow zeros for empty components

# update ns_trt
for (h in 1:H_trt) ns_trt[h]<-length(c_trt[c_trt==h])   #must allow zeros for empty components

# update v
for (h in 1:(H-1)) v[h]<-rbeta(1,1+ns[h],alpha+sum(ns[(h+1):H]))
v[H]=1
cumv<-cumprod(1-v)
pi[1]<-v[1]
for (h in 2:H) pi[h]<-v[h]*cumv[h-1]

# update v_trt
for (h in 1:(H_trt-1)) v_trt[h]<-rbeta(1,1+ns_trt[h],alpha_trt+sum(ns_trt[(h+1):H_trt]))
v_trt[H_trt]=1
cumv_trt<-cumprod(1-v_trt)
pi_trt[1]<-v_trt[1]
for (h in 2:H_trt) pi_trt[h]<-v_trt[h]*cumv_trt[h-1]

# update alpha
alpha<-rgamma(1,a_alpha+H-1,b_alpha-sum(log(1-v)[1:(H-1)]))   #b_alpha-log(cumv[H-1]) or b_alpha-log(pi[H])

# update alpha_trt
alpha_trt<-rgamma(1,a_alpha_trt+H_trt-1,b_alpha_trt-sum(log(1-v_trt)[1:(H_trt-1)]))   #b_alpha-log(cumv[H-1]) or b_alpha-log(pi[H])

# update tau
for (i in 1:I) tau[i]<-tau_star[c[i]]

# update tau_trt
for (i in 1:I) tau_trt[i]<-tau_trt_star[c_trt[i]]

# calculate -2logL (exact obs.)
f_iter_exact<-lambdatT*exp(xcov[1:n1,]%*%beta+phicoef2[1:n1]+phicoef2_trt[1:n1]*xtrt[1:n1])*exp(-LambdatT*exp(xcov[1:n1,]%*%beta+phicoef2[1:n1]+phicoef2_trt[1:n1]*xtrt[1:n1]))

# calculate -2logL (IC obs.)
FL<-1-exp(-LambdatL*exp(xcov[(n1+1):N,]%*%beta+phicoef2[(n1+1):N]+phicoef2_trt[(n1+1):N]*xtrt[(n1+1):N]))   # n2*1
FR<-1-exp(-LambdatR*exp(xcov[(n1+1):N,]%*%beta+phicoef2[(n1+1):N]+phicoef2_trt[(n1+1):N]*xtrt[(n1+1):N]))   # n2*1 
f_iter_IC<-(FR^(y[(n1+1):N]==0))*((FR-FL)^(y[(n1+1):N]==1))*((1-FL)^(y[(n1+1):N]==2)) # n2*1, individual likelihood for each iteration

finv_iter_exact<-1/f_iter_exact
finv_iter_IC<-1/f_iter_IC # n*1, inverse of individual likelihood for each iteration

loglike<-sum(log(f_iter_exact))+sum(log(FR^(y[(n1+1):N]==0))+log((FR-FL)^(y[(n1+1):N]==1))+log((1-FL)^(y[(n1+1):N]==2)))
dev<--2*loglike      # -2logL

      parbeta[iter,]=beta
      parbeta_original[iter,]=beta_original
      pareta[iter]=eta
partau_star[iter,]=tau_star
partau_trt_star[iter,]=tau_trt_star
paralpha[iter]=alpha
paralpha_trt[iter]=alpha_trt
parphi[iter,]=phicoef
parphi_trt[iter,]=phicoef_trt
parpi[iter,]=pi
parpi_trt[iter,]=pi_trt
      pargam[iter,]=gamcoef
      ttt=gamcoef%*%bisg
	if (scale.designX==FALSE) {parsurv0[iter,]<-exp(-ttt)}
      if (scale.designX==TRUE) {parsurv0[iter,]<-exp(-ttt*exp(-sum((beta*mean_X/sd_X)[scaled==1])))}
      parlambdatT[iter,]=lambdatT
      parLambdatT[iter,]=LambdatT
parLambdatL[iter,]=LambdatL
parLambdatR[iter,]=LambdatR
pardev[iter]=dev
parfinv_exact[iter,]=finv_iter_exact
parfinv_IC[iter,]=finv_iter_IC

      if (is.null(x_user)){parsurv[iter,]=parsurv0[iter,]} else {
      A<-matrix(x_user,byrow=TRUE,ncol=p)
      B<-exp(A%*%beta)
      for (g in 1:G){
      parsurv[iter,((g-1)*kgrids+1):(g*kgrids)]=exp(-ttt*B[g,1])}
      }

      iter=iter+1
if (iter%%100==0) print(iter)

   } # end iteration

   wbeta=as.matrix(parbeta_original[seq((burnin+thin),total,by=thin),],ncol=p) # thinned beta samples
   wparsurv0=as.matrix(parsurv0[seq((burnin+thin),total,by=thin),],ncol=kgrids)
   wparsurv=as.matrix(parsurv[seq((burnin+thin),total,by=thin),],ncol=kgrids*G)
   coef<-apply(wbeta,2,mean)
   coef_ssd<-apply(wbeta,2,sd)
   coef_ci<-array(rep(0,p*2),dim=c(p,2))
   S0_m<-apply(wparsurv0,2,mean)
   S_m<- apply(wparsurv,2,mean)
   colnames(coef_ci)<-c(paste(100*(1-conf.int)/2,"%CI"),paste(100*(0.5+conf.int/2),"%CI"))
   for (r in 1:p) coef_ci[r,]<-quantile(wbeta[,r],c((1-conf.int)/2,0.5+conf.int/2))

   CPO_exact=1/apply(parfinv_exact[seq((burnin+thin),total,by=thin),],2,mean)
   CPO_IC=1/apply(parfinv_IC[seq((burnin+thin),total,by=thin),],2,mean)
   NLLK_exact=-sum(log(CPO_exact))      # smaller is better
   NLLK_IC=-sum(log(CPO_IC))      # smaller is better
   NLLK=NLLK_exact+NLLK_IC

   #calculate D(theta_bar) for IC
   LambdatL_m<-apply(parLambdatL[seq((burnin+thin),total,by=thin),],2,mean)
   LambdatR_m<-apply(parLambdatR[seq((burnin+thin),total,by=thin),],2,mean)
   beta_m<-apply(as.matrix(parbeta[seq((burnin+thin),total,by=thin),]),2,mean)
   phicoef_m<-apply(parphi[seq((burnin+thin),total,by=thin),],2,mean)
   phicoef2_m<-array(rep(0,N),dim=c(N,1))
   for (j in 1:N){
   phicoef2_m[j]<-phicoef_m[area[j]]
   }
   phicoef_trt_m<-apply(parphi_trt[seq((burnin+thin),total,by=thin),],2,mean)
   phicoef2_trt_m<-array(rep(0,N),dim=c(N,1))
   for (j in 1:N){
   phicoef2_trt_m[j]<-phicoef_trt_m[area[j]]
   }
   FL_m<-1-exp(-LambdatL_m*exp(as.matrix(xcov[(n1+1):N,])%*%as.matrix(beta_m)+phicoef2_m[(n1+1):N]+phicoef2_trt_m[(n1+1):N]*xtrt[(n1+1):N]))
   FR_m<-1-exp(-LambdatR_m*exp(as.matrix(xcov[(n1+1):N,])%*%as.matrix(beta_m)+phicoef2_m[(n1+1):N]+phicoef2_trt_m[(n1+1):N]*xtrt[(n1+1):N]))
   loglike_m_IC<-sum(log(FR_m^(y[(n1+1):N]==0))+log((FR_m-FL_m)^(y[(n1+1):N]==1))+log((1-FL_m)^(y[(n1+1):N]==2)))
   D_thetabar_IC<--2*loglike_m_IC      # -2logL

   #calculate D(theta_bar) for exact
   lambdatT_m<-apply(parlambdatT[seq((burnin+thin),total,by=thin),],2,mean)
   LambdatT_m<-apply(parLambdatT[seq((burnin+thin),total,by=thin),],2,mean)   
   loglike_m_exact<-sum(log(lambdatT_m)+as.matrix(xcov[1:n1,])%*%as.matrix(beta_m)+phicoef2_m[1:n1]+phicoef2_trt_m[1:n1]*xtrt[1:n1]-LambdatT_m*exp(as.matrix(xcov[1:n1,])%*%as.matrix(beta_m)+phicoef2_m[1:n1]+phicoef2_trt_m[1:n1]*xtrt[1:n1]))
   D_thetabar_exact<--2*loglike_m_exact     # -2logL
 
   D_bar=mean(pardev[seq((burnin+thin),total,by=thin)]) 
   D_thetabar=D_thetabar_IC+D_thetabar_exact
   DIC=2*D_bar-D_thetabar

   est<-list(
   N=nrow(xcov),
   nameX=colnames(xcov),
   parbeta=parbeta_original,
   parsurv0=parsurv0,
   parsurv=parsurv,
   parphi=parphi,
   parphi_trt=parphi_trt,
   paralpha=paralpha,
   paralpha_trt=paralpha_trt,
   partau_star=partau_star,
   partau_trt_star=partau_trt_star,
   coef = coef,
   coef_ssd = coef_ssd,
   coef_ci = coef_ci,
   S0_m = S0_m,
   S_m = S_m,
   grids=grids,
   DIC = DIC,
   NLLK = NLLK
   )
   est

}
