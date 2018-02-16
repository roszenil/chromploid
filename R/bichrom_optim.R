#' Calculates the maximum likelihood estimate and the value of the likelihood at the MLE of a BiChroM and BiChroM reduced models
#'
#' @param phy.tree phylogenetic tree
#' @param tip.values A vector with chromosome numbers and/or binary trait states. See Zenil-Ferguson et al 2017 supplement for an explanation of BichroM model.
#' @param model options are "bichrom" or "reducedbichrom"
#' @param model.args List of arguments passed in each specified model other than parameters like size (default size is 50 chromosomes)
#' @param starting.val vector of parameters in log scale that help with the optimization of the likelihood
#' @param root.dist vector of probabilities for the states at the root 
#' @param optim.options See opts from nloptr package, optimization options for the likelihood function
#' @author Rosana Zenil-Ferguson
#' @importFrom nloptr nloptr 
#' @export

bichrom_optim<-function(phy.tree, tip.values, model="bichrom",model.args=list(size=50),starting.val, root.dist, optim.options=list("algorithm"= "NLOPT_LN_SBPLX","ftol_rel"=1e-08,"print_level"=1,"maxtime"=170000000, "maxeval"=1000)){
	if(model=="bichrom"){
	results<-rep(0,12)
	mle<-nloptr::nloptr(x0=starting.val,eval_f=chromploid_nllike,opts=optim.options, phy=phy.tree, 	tip.values=tip.values,pi.0=root.dist,Q.FUN=Q_bichrom, Q.ARGS = model.args)
	print(mle) 
results[1:10]<-mle$solution 
results[11]<-mle$objective 
results[12]<-mle$status
results<-as.dataframe(results)
names(results)<-c("lambda0","lambda1","mu0","mu1","rho0","rho1","q01","q10","epsilon0","epsilon1","nloglike","convergencestatus")}else{
		 results<-rep(0,11)
		mle<-nloptr::nloptr(x0=starting.val,eval_f=chromploid_nllike,opts=optim.options, phy=phy.tree, tip.values=tip.values, pi.0=root.dist,Q.FUN=Q_reducedbichrom, Q.ARGS = model.args)
		print.mle
		results[1:9]<-mle$solution
		results[10]<-mle$objective
		results[11]<-mle$status}
        results<-as.dataframe(results)
        names(results)<-c("lambda0","lambda1","mu0","mu1","rho","q01","q10","epsilon0","epsilon1","nloglike","convergencestatus")
return(results)
}


