SurvtoLR<-function(x){
    ## change Surv object to data.frame with L and R columns
    type<-attr(x,"type")
    if (type=="interval"){
        L<-R<-x[,1]
        R[x[,3]==0]<-0   #NA
        L[x[,3]==2]<-0
        R[x[,3]==3]<-x[x[,3]==3,2]
    } else if (type=="right") {
        L<-R<-x[,1]
        R[x[,2]==0]<-0   #NA
    } else {
    stop(paste("Surv obj type='",type,"' unrecognized",sep=""))
    }
    if (type=="interval"){
        y<-x[,1]
        y[x[,3]==0]<-2   
        y[x[,3]==1]<-3
        y[x[,3]==2]<-0
        y[x[,3]==3]<-1
    } else if (type=="right") {
        y<-x[,1]
        y[x[,2]==0]<-2 
        y[x[,2]==1]<-3     
    }
    out<-cbind(L,R,y)
    return(out)
}
