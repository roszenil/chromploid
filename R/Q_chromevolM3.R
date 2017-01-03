#' Calculates Q-matrix for chromevol M3 model (beta version)
#' @details Q_chomevolM3 determines the Q-matrix for a model of chromosome number change with duplications, dysploidy, and  demiploidy.
#'
#' @param log.theta vector of size 4 indicating parameters in ln scale for chromevol M3 the order of the parameters is (lambda, mu, rho, delta)
#' @param size Maximum number of chromosomes in the sample (recommended no more than 40)
#' @return Q a sparse matrix of size (size+10)
#' @export
Q_chromevolM3<-function(log.theta,size){
  #Parameters
  theta<-exp(log.theta)
  l.0<-theta[1]
  m.0<-theta[2]
  r.0<-theta[3]
  d.0<-theta[4]
  C<-size+10


  #Building the sparse matrix
  Q<-rep(0,C*C)
  Q<- matrix(Q,ncol=C)

  Q[1,1]<- -(l.0+r.0)
  Q[1,2]<- l.0+r.0
  aux1<-floor((C-1)/2)
  for(i in 2:aux1){
    Q[i,i]<- -(l.0+r.0+m.0+d.0)
    Q[i,(i-1)]<-d.0
    Q[i,(i+1)]<-l.0
    Q[i,(2*i)]<-r.0
    aux2<- ((i*1.5)%%1)
    aux3<-floor(i*1.5)
    if(aux2==0){
      Q[i, aux3]<-Q[i,aux3]+m.0
    }else{
      Q[i,aux3]<-Q[i,aux3]+m.0/2
      Q[i,(aux3+1)]<-Q[i,(aux3+1)]+m.0/2
  }}

  ##################################
  for(i in (aux1+1):(C-2)){
    Q[i,i]<- -(l.0+r.0+m.0+d.0)
    Q[i,(i-1)]<-d.0
    Q[i,(i+1)]<-l.0
    Q[i,C]<-r.0
    aux2<- ((i*1.5)%%1)
    aux3<-floor(i*1.5)
    if(aux3<C){
    if(aux2==0){
      Q[i, aux3]<-Q[i,aux3]+m.0
    }else{
      Q[i,aux3]<-Q[i,aux3]+m.0/2
      Q[i,(aux3+1)]<-Q[i,(aux3+1)]+m.0/2
    }
    }else{
    Q[i,C]<-Q[i,C]+m.0
  }}
  ###################################
  Q[C-1,C-2]<-d.0
  Q[C-1,C-1]<--(l.0+m.0+r.0+d.0)
  Q[C-1,C]<-l.0+m.0+r.0
  Q[C,(C-1)]<-d.0
  Q[C,C]<--d.0
  return(Q)
}
