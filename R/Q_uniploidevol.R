#' Calculates Q-matrix for univariate profile likelihoods for Ploidevol model
#'
#' @param log.pars vector of size 5 indicating values in ln of nuisance parameters (excluding one of interest for the profile likelihood)
#' @param log.theta0 vector of size 1 indicating values parameters of interest in BiChroM. Options are (alpha, beta,delta, epsilon, rho, omega)
#' @param param.name vector of size 1 indicating the names of parameters of interest. E.g. "epsilon"
#' @return Q a sparse matrix of size 11x11
#' @export

Q_uniploidevol<-function(log.pars, log.theta0, param.name){
		if(length(log.theta0)>1){
			stop("log.theta0 should be a single value")
		}else{
		index=c("alpha","beta","delta","epsilon","rho","omega")
		aux<-which(index==param.name)-1
		all.parms<-append(log.pars,log.theta0, after=aux)
	pos.guess <- exp(all.parms)
	alpha     <- pos.guess[1]
	beta      <- pos.guess[2]
	delta     <- pos.guess[3]
	epsilon   <- pos.guess[4]
	rho       <- pos.guess[5]
	omega<- pos.guess[6]
	beta.subQ <- diag(rep(c(0,beta), 5));
	betamat <- matrix(0, nrow=11,ncol=11);
	betamat[2:11, 1:10] <- beta.subQ;
	omega.subQ<-diag(rep(c(omega, 0),5));
	betamat[1:10, 2:11]<-betamat[1:10, 2:11]+omega.subQ;

	alpha.subQ <- diag(rep(c(0,alpha),5));
	eps.subsub <- diag(rep(c(0,epsilon), 4));
	zero10 <- matrix(0,nrow=10,ncol=10);
	zero10[2:9,3:10] <- eps.subsub;
	alpha.subQ <- alpha.subQ + zero10;

	alphamat <- rbind(alpha.subQ,rep(0,10))
	alphamat<- cbind(rep(0,11),alphamat)

	Qmat <- betamat + alphamat;
	Qmat[1,3] <- epsilon + rho;
	Qmat[c(3,5,7,9,11),1] <- rep(delta,5);
	Qmat[3,7] <- rho;
	Qmat[5,11]<-rho;
#Qmat[5,11] <- rho; this were the transitions from 6x to 12x
#Qmat<-rbind(rep(0,11),Qmat) this were adding extinction to odd-ploids
#Qmat<-cbind((c(0,rep(c(0,mu),5),0)),Qmat)
	diag.Q <- apply(Qmat,1,sum)
	diag(Qmat) <- -diag.Q
	return(Qmat)
}
}
