summary.PICBayes <-
function(object, ...){
   coef.table<-cbind(object$coef, object$coef_ssd, object$coef_ci[,1], object$coef_ci[,2])
   dimnames(coef.table)<-list(object$nameX,c('Mean','Std. Dev.','CI_lower','CI_upper'))

   cat("\nPosterior inference of regression coefficients\n")
   print.default(coef.table)

   cat("\nDeviance information criterion: DIC=", object$DIC, sep="")
   cat("\nNegative log-likelihood: NLLK=", object$NLLK, sep="")
   cat("\nNumber of subjects: N=", object$N, "\n", sep="")
}
