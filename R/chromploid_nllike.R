#' Calculates the negative log-likelihood of a chromosome or ploidy model
#'
#' @param log.par A vector of parameter values in ln scale. Length of the vector depends on model defined by Q.FUN
#' @param phy Binary rooted phylogenetic tree
#' @param tip.values A data frame with at least two columns first column has tip labels in the same order than phylogenetic tree, second column ploidy or chromosome numbers depending on Q.FUN
#' @param pi.0 Probability distribution to be used at the root of phylogeny
#' @param Q.FUN A function that calculates the Q-matrix of the parameters of interest
#' @param Q.ARGS A list of arguments to be passed to Q.FUN
#' @return Result is the negative log-likelihood of parameters log.par given the model in Q.FUN
#' @author Rosana Zenil-Ferguson
#' @export
chromploid_nllike <- function(log.par, phy, tip.values, pi.0, Q.FUN, Q.ARGS) {

    if (!is.list(Q.ARGS)) stop("Q.ARGS must be a list")
    Q.FUN <- match.fun(Q.FUN)
    Q.mat <-do.call(Q.FUN, args=c(list(log.par),Q.ARGS)) #do.call(Q.FUN, args = list(log.par, unlist(Q.ARGS)))

    if (any(abs(rowSums(Q.mat)) > 1e-8)==TRUE){
		print("At least one row in the Q-matrix does not add to zero")
      break
	}
	charnum=1
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nl <- nrow(Q.mat)
	#Now we need to build the matrix of likelihoods to pass to dev.raydisc:
	liks <- matrix(0, nb.tip + nb.node, nl)
	#Now loop through the tips.
	for(i in 1:nb.tip){
		#The codon at a site for a species is not NA, then just put a 1 in the appropriate column.
		#Note: We add charnum+1, because the first column in the data is the species labels:
		if(!is.na(tip.values[i,charnum+1])){
			liks[i,tip.values[i,charnum+1]] <- 1
		}else{
			#If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can modified later:
			liks[i,] <- 1
		}
	}
	#The result here is just the likelihood:
	result <- pruning_like(p=NULL, phy=phy, liks=liks, Q=Q.mat, rate=NULL, root.p=pi.0)
	return(result)
}
