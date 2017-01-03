#' Calculates Q-matrix for Ploidevol  model
#' @details This function specifies a Q-matrix of size 11x11 where states are ploidy values from 2x to 12x
#'
#' @param log.theta vector of size 6 with ln values of parameters (alpha, beta, delta, epsilon, rho, omega)
#' @return Q a sparse matrix of size 11x11
#' @export
Q_ploidevol <- function(log.theta){
	# all.parms <- [alpha, beta, delta, epsilon, rho, omega] in log space

	pos.guess <- exp(log.theta)
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
