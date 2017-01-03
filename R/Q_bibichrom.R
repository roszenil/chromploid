#' Calculates Q-matrix for bivariate profile likelihoods for BiChroM model
#'
#' @param log.pars vector of size 8 indicating values in ln of nuisance parameters (excluding the ones of interest for the profile likelihood)
#' @param log.theta0 vector of size 2 indicating values parameters of interest in BiChroM. Options are (lambda0, lambda1, mu0, mu1, rho0, rho1, q01,q10, e01, e10)
#' @param size Maximum number of chromosomes in the sample (recommended no more than 50, states larger than that should be coded as 51)
#' @param param.name vector of size 2 indicating the names of parameters of interest. E.g. c("rho0","rho1")
#' @return Q a sparse matrix of size 2*(size+1)
#' @export

Q_bibichrom<-function(log.pars,log.theta0,size, param.name){
#Parameters
if(length(log.theta0)!=2){
	stop("log.theta0 should be a vector of dimension 2")
}else{
	index=c("lambda0","lambda1","mu0","mu1","rho0","rho1", "q01","q10","e0","e1")
	aux1<-which(index==param.name[1])-1
	all.parms<-append(log.pars,log.theta0[1], after=aux1)
	aux2<-which(index==param.name[2])-1
	all.parms<-append(all.parms,log.theta0[2],after=aux2)
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
