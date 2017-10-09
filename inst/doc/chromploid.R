## ----error=FALSE, message=FALSE, warning=FALSE---------------------------
library("devtools")
devtools::install_github("roszenil/chromploid")
library("chromploid")
library("geiger") 
library("nloptr")
#chromploid depends on geiger and nloptr so don't forget to load the packages

## ---- message=FALSE, warning=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=60)----
log.params<-log(c(0.12, 0.001, 0.25, 0.002,0.036, 0.006, 0.04,0.02, 1.792317852, 1.57E-14))
N<- 50
mymatrix<-Q_bichrom(log.theta=log.params, size=N)

## ------------------------------------------------------------------------
mytree<- sim.bdtree(b=0.055, stop="taxa", n=500, seed=1) 
# Seed was fixed to obtain always the same tree for this tutorial
plot(mytree,type="fan",cex=0.2, no.margin=TRUE)
set.seed(1) #Seed was fixed to obtain always the same sample for this tutorial
mysample<- chromploid_sim.char(phy=mytree, par=mymatrix,model="discrete", root=56)

## ---- tidy=TRUE, tidy.opts=list(width.cutoff=60)-------------------------
nllike<-chromploid_nllike(log.par=log.params, phy=mytree,   tip.values=mysample, pi.0=NULL, Q.FUN=Q_bichrom, Q.ARGS = list(size=50))
nllike

## ---- eval=FALSE,tidy=TRUE, tidy.opts=list(width.cutoff=60)--------------
#  # Error in tip.values[i, charnum + 1] : incorrect number of dimensions

## ---- eval=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=60)-------------
#  mysample<-data.frame(taxa=rownames(mysample),sample=as.numeric(mysample))
#  head(mysample)

## ---- eval=FALSE,tidy=TRUE, tidy.opts=list(width.cutoff=60)--------------
#  # Optimizations take hours to run, I suggest you run this routine in a cluster. For a 4700 tip tree starting with a x.0 close to the maximum it can take up to 96 hrs.
#  
#  library(nloptr)
#  x.0<- log(c(0.12, 0.001, 0.25, 0.002,0.036, 0.006, 0.04,0.02, 1.792317852, 1.57e-14))
#  #value where the optimization algorithm is going to start. A vector of 10 values for the parameters in log scale
#  
#  model.args<-list(size=50)
#  optimization.bichrom<-chromploid_optim(phy.tree=mytree,tip.values=mysample,model="bichrom", model.args=model.args, starting.val=x.0, root.dist=NULL)
#  
#  print(optimization.bichrom)

## ---- eval=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=60)-------------
#  x.0<- log(c(0.12, 0.001, 0.25, 0.002,0.010, 0.04,0.02, 1.792317852, 1.57e-14))
#  # Value where the optimization algorithm is going to start. A vector of 9 values for the parameters in log scale. The value of the hypothesis here is rho=0.010
#  
#  model.args=list(size=50, equal.param=c("rho0","rho1"),location.par=5)
#  #Q.ARGS is a list that has all the arguments included in Q_reducedbichrom that are not the parameters
#  #size= maximum number of haploid chromosome numbers to consider in your sample,
#  #equal.params=which parameters are equal based on the hypothesis H0,
#  #location.par is the position in the vector where rho value appears
#  
#  optimization.reducedbichrom<-chromploid_optim(phy.tree=mytree, tip.values=mysample, model="reducedbichrom",model.args=model.args, starting.val=x.0,root.dist=NULL)

## ---- eval=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=60)-------------
#  neglog.red<-optimization.reducedbichrom[10]
#  neglog.full<-optimization.bichrom[11]
#  alpha<-0.05
#  D =2*(neglog.red-neglog.full)
#  p.value<-pchisq(D, lower.tail=FALSE, df=1)
#  	if(p.value>0.05){
#  		reject<-0
#  	}else{reject<-1}

## ---- eval=FALSE, tidy=TRUE, tidy.opts=list(width.cutoff=60)-------------
#  logrho0.values<-log(seq(0.02,0.04,0.005))
#  # Calculate the likelihood values in a grid 0.02, 0.025, ..., 0.035, 0.04 (in log scale). True value is around 0.36 so values surrounding the true value are interesting
#  long<-length(logrho0.values)
#  x.0<- log(c(0.12, 0.001, 0.25, 0.002, 0.006, 0.04,0.02, 1.792317852, 1.57e-14))
#  # value where the optimization algorithm is going to start. A vector of 9 values for the nuisance parameters in log scale (everything except rho0)
#  profile.values<-rep(0,long) # Vector useful to save optimization results
#  my.options<-list("algorithm"= "NLOPT_LN_SBPLX","ftol_rel"=1e-08,"print_level"=1,"maxtime"=170000000, "maxeval"=1000)
#  # Options used in nloptr, for more information ?nloptr. This is a subplex algorithm, with high tolerance to fine a very precise optimum.
#  for(i in 1:long){
#    mle.profilerho0<-nloptr(x0=x.0,eval_f=chromploid_nllike,opts=my.options, phy=mytree, tip.values=mysample,pi.0=NULL,Q.FUN=Q_unibichrom, Q.ARGS = list(log.theta0=logrho0.values[i],size=50, param.name="rho0"))
#  #Q.ARGS is a list that has all the arguments included in Q_unibichrom
#  #Type ?Q_unibichrom to see a longer explanation of each argument
#    print(mle.profilerho0)
#    profile.values[i]<-mle.profilerho0$objective
#  }
#  
#  plot(exp(logrho0.values),profile.values, xlab=expression(rho_0),ylab="Negative log-likelihood", type="l", lwd=2)

