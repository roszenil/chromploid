#' Calculates a Q matrix for reduced BiChroM model
#' @details This function assumes that two parameter values are equal. This is useful for calculating likelihood ratio tests in BiChrom. By default the two parameters assumed equal are "rho0"="rho1" but user can specify which two parameters are the same
#'
#' @param log.pars a vector of size 9 with  ln values of parameters (lambda0, lambda1, mu0, mu1, rho, q01, q10, e0, e1 )
#' @param size Maximum number of chromosomes in the sample (recommended no more than 50, states larger than that should be coded as 51)
#' @param equal.param by default "rho0"="rho1" but user can specify the parameters forced to be equal
#' @param location position in log.pars where the value of the equal.param is located, by default location=5
#' @export
Q_reducedbichrom<-function(log.pars,size, equal.param=c("rho0","rho1"), location.par=5){
#Parameters
if(length(equal.param)!=2){
	stop("reduced model requires specifying which two parameters are equal")
}else{
	if (length(log.pars)!=9){
		stop("log.pars should be a vector of size nine")
	}else{
	index=c("lambda0","lambda1","mu0","mu1","rho0","rho1", "q01","q10","e0","e1")
	log.pars<-log.pars[-location.par]
	aux1<-which(index==equal.param[1])-1
	all.parms<-append(log.pars,log.pars[location.par], after=aux1)
	aux2<-which(index==equal.param[2])-1
	all.parms<-append(all.parms,log.pars[location.par],after=aux2)
theta<-exp(all.parms)
l.0<-theta[1]
l.1<-theta[2]
m.0<-theta[3]
m.1<-theta[4]
r.0<-theta[5]
r.1<-theta[6]
prob.01<-theta[7]
prob.10<-theta[8]
e.0<-theta[9]
e.1<-theta[10]
C<-size+1

#Building the sparse matrix
Q<-rep(0,4*C*C)
Q<- matrix(Q,ncol=2*C)
aux1<- floor(size/2)
aux2<- 2*C-1 ###############????????????

Q[1,1]<- -(l.0+r.0+ prob.01)
Q[1,2]<- l.0+r.0
Q[1,(C+1)]<-prob.01

for(i in 2:aux1){
	Q[i,i]<- -(m.0+l.0+r.0+prob.01)
	Q[i,(i-1)]<-m.0
	Q[i,(i+1)]<-l.0
	Q[i,(2*i)]<-r.0
	Q[i,(C+i)]<- prob.01
	}


for (i in (aux1+1):(C-1)){
	Q[i,i]<- -(m.0+l.0+r.0+prob.01)
	Q[i,(i-1)]<- m.0
	Q[i,(i+1)]<- l.0
	Q[i,C]<- r.0+Q[i,C]
	Q[i,(C+i)]<- prob.01
	}


Q[C,C]<- -(e.0)
Q[C,(2*C)]<-e.0

###################################################

Q[(C+1),(C+1)]<- -(l.1+r.1+prob.10)
Q[(C+1),(C+2)]<- l.1+r.1
Q[(C+1),1]<-prob.10

for(i in (C+2):(C+aux1)){
	Q[i,i]<- -(m.1+l.1+r.1+prob.10)
	Q[i,(i-1)]<-m.1
	Q[i,(i+1)]<-l.1
	Q[i,(2*i-C)]<-r.1
	Q[i,(i-C)]<- prob.10
	}


for(i in (C+aux1+1):aux2){
	Q[i,i]<- -(m.1+r.1+l.1+prob.10)
	Q[i,(i-1)]<-m.1
	Q[i,(i+1)]<-l.1
	Q[i,2*C]<- r.1+Q[i,2*C]
	Q[i,(i-C)]<- prob.10
	}

Q[(2*C),(2*C)]<- -(e.1)
Q[(2*C),C]<-e.1
return(Q)
}
}
}
