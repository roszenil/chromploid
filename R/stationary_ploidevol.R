#' Calculates stationary distribution for root under Ploidevol model
#' @author Rosana Zenil-Ferguson
#'
#' @param log.parameters a vector of size 6 containing ln values for (alpha, beta, delta, epsilon, rho, omega)
#' @return A probability vector of size 11 indicating probabilities of the root being ploidy 2x, 3x,..., 12x
#' @seealso
#' \code{\link{negloglikelihood}}, \code{\link{Q_ploidevol}}
#' @export
stationary_ploidevol<-function(log.parameters){
	Q_aa<- Q_ploidevol(log.parameters)
	p.stationary <- Null(Q_aa)
	p.stationary <- p.stationary/sum(p.stationary)
	return(p.stationary)
}
