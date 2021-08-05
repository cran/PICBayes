clusterIC_Z <-
function(
L, 
R, 
y, 
xcov,
IC,
scale.designX,
scaled,
zcov,
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
a_tau,
b_tau,
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
   zcov=as.matrix(zcov)
   area=matrix(area,ncol=1)
   IC=matrix(IC,ncol=1)
   p=ncol(xcov)
   q=ncol(zcov)
   N=nrow(L)

   if (scale.designX==TRUE){
   mean_X<-apply(xcov,2,mean)
   sd_X<-apply(xcov,2,sd)
   for (r in 1:p){ 
   if (scaled[r]==1) xcov[,r]<-(xcov[,r]-mean_X[r])/sd_X[r]
   }
   }

   ## generate basis functions
   K=length(knots)-2+order
   kgrids=length(grids)

   bisL=Ispline(L,order,knots) # K*n2
   bisR=Ispline(R,order,knots)
   bisg=Ispline(grids,order,knots)

   ## initial value
   eta=rgamma(1,a_eta,rate=b_eta)
   tau=rgamma(q,a_tau,rate=b_tau)
   gamcoef=matrix(rgamma(K, 1, rate=1),ncol=K)

   phicoef<-matrix(rep(0,I*q),ncol=q)
   phicoef2<-matrix(rep(0,N*q),ncol=q)
   for (j in 1:N){
   phicoef2[j,]<-phicoef[area[j],]  #a vector of phi values for all subjects
   }

   beta=matrix(rep(0,p),p,1)
   beta_original=matrix(rep(0,p),ncol=1)

   LambdatL=t(gamcoef%*%bisL) # n2 x 1 # lambda0 at each u_i value
   LambdatR=t(gamcoef%*%bisR) # n2 x 1 # lambda0 at each v_i value
   Lambdatg=t(gamcoef%*%bisg)  # lambda0 at each grids_i value 

   parbeta=array(rep(0,total*p),dim=c(total,p))
   parbeta_original=array(rep(0,total*p),dim=c(total,p))
   pareta=array(rep(0,total),dim=c(total,1)) 
   partau=array(rep(0,total*q),dim=c(total,q))
   parphi=array(rep(0,I*q*total),dim=c(I,q,total))
   pargam=array(rep(0,total*K),dim=c(total,K))
   parsurv0=array(rep(0,total*kgrids),dim=c(total,kgrids))
   parLambdatL=array(rep(0,total*N),dim=c(total,N))
   parLambdatR=array(rep(0,total*N),dim=c(total,N))
   pardev=array(rep(0,total),dim=c(total,1))
   parfinv=array(rep(0,total*N),dim=c(total,N))

   if (is.null(x_user)){parsurv=parsurv0} else {
   G<-length(x_user)/p
   parsurv=array(rep(0,total*kgrids*G),dim=c(total,kgrids*G))}
      
   ## iteration
   iter=1
   while (iter<total+1)
   { 
      # sample z, zz, w and ww
      z=array(rep(0,N),dim=c(N,1)); w=z
      zz=array(rep(0,N*K),dim=c(N,K)); ww=zz

      for (j in 1:N){
      if (y[j]==0){
      templam1=LambdatR[j]*exp(xcov[j,]%*%beta+zcov[j,]%*%phicoef[area[j],])
      z[j]=poissrndpositive(templam1)
      zz[j,]=rmultinom(1,z[j],gamcoef*t(bisR[,j]))
      } else if (y[j]==1){
      templam1=(LambdatR[j]-LambdatL[j])*exp(xcov[j,]%*%beta+zcov[j,]%*%phicoef[area[j],])
      w[j]=poissrndpositive(templam1)
      ww[j,]=rmultinom(1,w[j],gamcoef*t(bisR[,j]-bisL[,j]))
      }
      }

      # sample beta
      te1=z*(y==0)+w*(y==1)   #n2 x 1
      te2=(LambdatR*(y==0)+LambdatR*(y==1)+LambdatL*(y==2))   #n2 x 1

#sample beta MH-sampler
for (r in 1:p){
if (binary[r]==0){
beta1<-beta2<-beta
if (iter<beta_iter) sd_cand<-beta_cand[r] else sd_cand<-sd(parbeta[1:(iter-1),r])   #specify sd for candidate distribution
xt<-beta[r] #current state
yt<-rnorm(1,xt,sd_cand) #candidate point
beta1[r]<-yt
beta2[r]<-xt
log_f1<-sum(yt*xcov[,r]*te1)-sum(exp(xcov%*%beta1+apply(zcov*phicoef2,1,sum))*te2)-0.5*(yt^2)/(beta_sig0^2)
log_f2<-sum(xt*xcov[,r]*te1)-sum(exp(xcov%*%beta2+apply(zcov*phicoef2,1,sum))*te2)-0.5*(xt^2)/(beta_sig0^2)
num<-log_f1 #proposal dist is symmetric
den<-log_f2
if (log(runif(1))<(num-den)) beta[r]<-yt else beta[r]<-xt
}

if (binary[r]==1 & p>1){
te4=sum(xcov[,r]*te1)
te5=sum(te2*exp(as.matrix(xcov[,-r])%*%as.matrix(beta[-r])+apply(zcov*phicoef2,1,sum))*xcov[,r])
beta[r]<-log(rgamma(1,a_ga+te4,rate=b_ga+te5))
}

if (binary[r]==1 & p==1){
te4=sum(xcov[,r]*te1)
te5=sum(te2*exp(apply(zcov*phicoef2,1,sum))*xcov[,r])
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
         tempa=1+sum(zz[,l]*(y==0)+ww[,l]*(y==1))
         tempb=eta+sum(((bisR[l,])*(y==0)+(bisR[l,])*(y==1)
+(bisL[l,])*(y==2))*exp(xcov%*%beta+apply(zcov*phicoef2,1,sum)))
         gamcoef[l]=rgamma(1,tempa,rate=tempb)
      }

      LambdatL=t(gamcoef%*%bisL) # n2 x 1
      LambdatR=t(gamcoef%*%bisR) # n2 x 1

      #sample eta
      eta=rgamma(1,a_eta+K, rate=b_eta+sum(gamcoef))

#sample phi MH-sampler
for (r in 1:q){
phi1<-array(rep(0,N*q),dim=c(N,q))
phi2<-array(rep(0,N*q),dim=c(N,q))
for (i in 1:I){
phi1<-phi2<-phicoef2
if (iter<phi_iter) sd_cand<-phi_cand else sd_cand<-sd(parphi[i,r,1:(iter-1)]) #specify sd for candidate distribution
xt<-phicoef[i,r] #current state
yt<-rnorm(1,xt,sd_cand) # candidate point
phi1[area==i,r]<-yt
phi2[area==i,r]<-xt
log_f1<-sum(zcov[,r]*phi1[,r]*te1)-sum((exp(xcov%*%beta+apply(zcov*phi1,1,sum)))*te2)-0.5*tau[r]*(yt^2)
log_f2<-sum(zcov[,r]*phi2[,r]*te1)-sum((exp(xcov%*%beta+apply(zcov*phi2,1,sum)))*te2)-0.5*tau[r]*(xt^2)
num<-log_f1 #proposal dist is symmetric
den<-log_f2
if (log(runif(1))<(num-den)) phicoef[i,r]<-yt else phicoef[i,r]<-xt
phicoef2[area==i,r]<-phicoef[i,r]
}
}

#sample tau
for (r in 1:q){
tau[r]<-rgamma(1,0.5*I+a_tau,rate=0.5*t(phicoef[,r])%*%phicoef[,r]+b_tau) 
}

#calculate -2logL
FL<-1-exp(-LambdatL*exp(xcov%*%beta+apply(zcov*phicoef2,1,sum)))   # n2*1
FR<-1-exp(-LambdatR*exp(xcov%*%beta+apply(zcov*phicoef2,1,sum)))   # n2*1 
f_iter<-(FR^(y==0))*((FR-FL)^(y==1))*((1-FL)^(y==2)) # n2*1, individual likelihood for each iteration

finv_iter<-1/f_iter # n*1, inverse of individual likelihood for each iteration

loglike<-sum(log(FR^(y==0))+log((FR-FL)^(y==1))+log((1-FL)^(y==2)))
dev<--2*loglike      # -2logL

      parbeta[iter,]=beta
      parbeta_original[iter,]=beta_original
      pareta[iter]=eta
partau[iter,]=tau
parphi[,,iter]=phicoef
      pargam[iter,]=gamcoef
      ttt=gamcoef%*%bisg
	if (scale.designX==FALSE) {parsurv0[iter,]<-exp(-ttt)}
      if (scale.designX==TRUE) {parsurv0[iter,]<-exp(-ttt*exp(-sum((beta*mean_X/sd_X)[scaled==1])))}
parLambdatL[iter,]=LambdatL
parLambdatR[iter,]=LambdatR
pardev[iter]=dev
parfinv[iter,]=finv_iter

      if (is.null(x_user)){parsurv[iter,]=parsurv0[iter,]} else {
      A<-matrix(x_user,byrow=TRUE,ncol=p)
	if (scale.designX==TRUE){
	for (r in 1:p){ 
	if (scaled[r]==1) A[,r]<-(A[,r]-mean_X[r])/sd_X[r]}
	}
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

   CPO=1/apply(parfinv[seq((burnin+thin),total,by=thin),],2,mean)
   NLLK=-sum(log(CPO))      # smaller is better

   #calculate D(theta_bar) 
   LambdatL_m<-apply(parLambdatL[seq((burnin+thin),total,by=thin),],2,mean)
   LambdatR_m<-apply(parLambdatR[seq((burnin+thin),total,by=thin),],2,mean)
   beta_m<-apply(as.matrix(parbeta[seq((burnin+thin),total,by=thin),]),2,mean)
   phicoef_m<-apply(parphi[,,seq((burnin+thin),total,by=thin)],c(1,2),mean)
   phicoef2_m<-array(rep(0,N*q),dim=c(N,q))
   for (j in 1:N){
   phicoef2_m[j,]<-phicoef_m[area[j],]
   }
   FL_m<-1-exp(-LambdatL_m*exp(as.matrix(xcov)%*%as.matrix(beta_m)+apply(zcov*phicoef2_m,1,sum)))
   FR_m<-1-exp(-LambdatR_m*exp(as.matrix(xcov)%*%as.matrix(beta_m)+apply(zcov*phicoef2_m,1,sum)))
   loglike_m<-sum(log(FR_m^(y==0))+log((FR_m-FL_m)^(y==1))+log((1-FL_m)^(y==2)))
   D_thetabar<--2*loglike_m      # -2logL
 
   D_bar=mean(pardev[seq((burnin+thin),total,by=thin)]) 
   DIC=2*D_bar-D_thetabar

   est<-list(
   N=nrow(xcov),
   nameX=colnames(xcov),
   parbeta=parbeta_original,
   parsurv0=parsurv0,
   parsurv=parsurv,
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
