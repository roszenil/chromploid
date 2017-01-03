#' Calculates Q-matrix for BiChroM model
#' @details Q_bichromlinear determines the Q-matrix for a model of chromosome number change associated with a binary state with traits coded as 0 or 1. It is linear because binary trait change rate is a linear function of chromosome number.
#' @details Note: This is a beta version, the linear function for binary trait change rate from state 0 to 1 is  kappa0+ (chromosome number)*q01
#' @param log.theta vector of size 12 indicating parameters in ln scale for BiChoM the order of the parameters is (lambda0, lambda1, mu0, mu1, rho0, rho1, kappa0, q01, kappa1, q10, e0, e1)
#' @param size Maximum number of chromosomes in the sample (recommended no more than 50, states larger than that should be coded as 51)
#' @return Q a sparse matrix of size 2*(size+1)
#' @export
Q_bichromlinear<-function(log.theta,size){
  #Parameters
  theta<-exp(log.theta)
  l.0<-theta[1]
  l.1<-theta[2]
  m.0<-theta[3]
  m.1<-theta[4]
  r.0<-theta[5]
  r.1<-theta[6]
  k.0<-theta[7]
  k.1<-theta[8]
  prob.01<-theta[9]
  prob.10<-theta[10]
  e.0<-theta[11]
  e.1<-theta[12]
  C<-size+1


  #Building the sparse matrix
  Q<-rep(0,4*C*C)
  Q<- matrix(Q,ncol=2*C)
  aux1<- floor(size/2)
  aux2<- 2*C-1 ###############????????????

  Q[1,1]<- -(l.0+r.0+ k.0+prob.01)
  Q[1,2]<- l.0+r.0
  Q[1,(C+1)]<-k.0+prob.01

  for(i in 2:aux1){
    Q[i,i]<- -(m.0+l.0+r.0+k.0+i*prob.01)
    Q[i,(i-1)]<-m.0
    Q[i,(i+1)]<-l.0
    Q[i,(2*i)]<-r.0
    Q[i,(C+i)]<- k.0+i*prob.01
  }


  for (i in (aux1+1):(C-1)){
    Q[i,i]<- -(m.0+l.0+r.0+k.0+i*prob.01)
    Q[i,(i-1)]<- m.0
    Q[i,(i+1)]<- l.0
    Q[i,C]<- r.0+Q[i,C]
    Q[i,(C+i)]<- k.0+i*prob.01
  }


  Q[C,C]<- -(e.0)
  Q[C,(2*C)]<-e.0

  ###################################################

  Q[(C+1),(C+1)]<- -(l.1+r.1+k.1+prob.10)
  Q[(C+1),(C+2)]<- l.1+r.1
  Q[(C+1),1]<-k.1+prob.10

  for(i in (C+2):(C+aux1)){
    Q[i,i]<- -(m.1+l.1+r.1+k.1+(i-C)*prob.10)
    Q[i,(i-1)]<-m.1
    Q[i,(i+1)]<-l.1
    Q[i,(2*i-C)]<-r.1
    Q[i,(i-C)]<- k.1+(i-C)*prob.10
  }


  for(i in (C+aux1+1):aux2){
    Q[i,i]<- -(m.1+l.1+r.1+k.1+(i-C)*prob.10)
    Q[i,(i-1)]<-m.1
    Q[i,(i+1)]<-l.1
    Q[i,2*C]<- r.1+Q[i,2*C]
    Q[i,(i-C)]<- k.1+(i-C)*prob.10
  }

  Q[(2*C),(2*C)]<- -(e.1)
  Q[(2*C),C]<-e.1
  return(Q)
}
