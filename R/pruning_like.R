# Most recent version of Jeremy Beaulieu's prunning algorithm with calculation of distribution for the root 11/11/15 
# Modified by RZF to correctly calculate exponential of matrices, internally used by negloglike
pruning_like<-function(p,phy,liks,Q,rate,root.p){



	nb.tip <- length(phy$tip.label)

	nb.node <- phy$Nnode

	TIPS <- 1:nb.tip

	comp <- numeric(nb.tip + nb.node)

	phy <- reorder(phy, "pruningwise")

	#Obtain an object of all the unique ancestors

	anc <- unique(phy$edge[,1])


	#This bit is to allow packages like "selac" the ability to deal with this function directly:

	if(is.null(rate)){

		Q=Q

	}else{

		if (any(is.nan(p)) || any(is.infinite(p))) return(1000000)

		Q[] <- c(p, 0)[rate]

		diag(Q) <- -rowSums(Q)

	}



	for (i  in seq(from = 1, length.out = nb.node)) {

		#the ancestral node at row i is called focal

		focal <- anc[i]

		#Get descendant information of focal

		desRows<-which(phy$edge[,1]==focal)

		desNodes<-phy$edge[desRows,2]

		v <- 1



		for (desIndex in sequence(length(desRows))){

			v<-v*expm:::expm(Q * phy$edge.length[desRows[desIndex]], method=c("Higham08")) %*% liks[desNodes[desIndex],]

		}

		comp[focal] <- sum(v)

		liks[focal, ] <- v/comp[focal]

	}

	#Specifies the root:

	root <- nb.tip + 1L

	#If any of the logs have NAs restart search:

	if (is.na(sum(log(comp[-TIPS])))){return(1000000)}

	else{

        equil.root <- NULL

        for(i in 1:ncol(Q)){

            posrows <- which(Q[,i] >= 0)

            rowsum <- sum(Q[posrows,i])

            poscols <- which(Q[i,] >= 0)

            colsum <- sum(Q[i,poscols])

            equil.root <- c(equil.root,rowsum/(rowsum+colsum))

        }

        if (is.null(root.p)){

            flat.root = equil.root

            k.rates <- 1/length(which(!is.na(equil.root)))

            flat.root[!is.na(flat.root)] = k.rates

            flat.root[is.na(flat.root)] = 0

            loglik<- -(sum(log(comp[-TIPS])) + log(sum(flat.root * liks[root,])))

        }

        else{

            if(is.character(root.p)){

                # root.p==yang will fix root probabilities based on the inferred rates: q10/(q01+q10)

                if(root.p == "yang"){

                    diag(Q) = 0

                    equil.root = colSums(Q) / sum(Q)

                    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(equil.root)+log(liks[root,])))))

                    if(is.infinite(loglik)){

                        return(1000000)

                    }

                }else{

                    # root.p==maddfitz will fix root probabilities according to FitzJohn et al 2009 Eq. 10:

                    root.p = liks[root,] / sum(liks[root,])

                    loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))

                }

            }

            # root.p!==NULL will fix root probabilities based on user supplied vector:

            else{

                loglik <- -(sum(log(comp[-TIPS])) + log(sum(exp(log(root.p)+log(liks[root,])))))

                if(is.infinite(loglik)){

                    return(1000000)

                }

            }

        }

	}

	loglik

}
