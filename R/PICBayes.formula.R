PICBayes.formula<-function (formula, data, ...) 
{
    # 'data' must be data.frame
    mf<-model.frame(formula=formula,data=data)  
    response<-model.response(mf)
    LR<-matrix(SurvtoLR(response),ncol=3)
    x<-model.matrix(attr(mf,'terms'),data=mf)  # get design matrix
    xcov<-x[,-1]  # delete intercept column
    est<-PICBayes.default(L=LR[,1],R=LR[,2],y=LR[,3],xcov=xcov,...)
    est$call <- match.call()
    est$formula <- formula
    est
}
