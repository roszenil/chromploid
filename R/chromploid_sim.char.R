#' Simulates a a sample of chromosome number or ploidy change for models with large sparse matrices
#'
#' @param phy Bifurcating phylogenetic tree
#' @param par Q-matrix model to be used
#' @param nsim number of simulations by default nsim=1
#' @param model by default model="discrete" for other usage see geiger:::sim.char
#' @param root discrete value at the root, by default 1
#' @usage For BiChroM simulation sim.char(phy=mytree, par=Q_bichrom, nsim=1, model="discrete", root=56)
#' @export

chromploid_sim.char<-function (phy, par, nsim = 1, model = "discrete", root = 1)
{
  model = match.arg(model, c("BM", "speciational", "discrete"))
  model.matrix = geiger:::.make.modelmatrix(par, model)
  nbranches <- nrow(phy$edge)
  nspecies <- Ntip(phy)
  if (length(root) > 1)
    stop("'root' should be a single value")
  if (model %in% c("BM", "speciational")) {
    m <- .get.simulation.matrix(phy)
    if (model == "speciational") {
      m[m > 0] <- 1
    }
    nchar <- nrow(model.matrix)
    rnd <- t(mvrnorm(nsim * nbranches, mu = rep(0, nchar),
                     Sigma = model.matrix))
    rnd <- array(rnd, dim = c(nchar, nbranches, nsim))
    simulate <- function(v, root) (m %*% as.matrix(v)) +
      root
    result <- apply(rnd, 1, simulate, root)
    result <- aperm(array(result, dim = c(nspecies, nsim,
                                          nchar)), c(1, 3, 2))
    rownames(result) <- phy$tip.label
  }
  else {
    rt = nspecies + 1
    zphy = reorder.phylo(phy, "postorder")
    el = zphy$edge.length
    nchar <- length(model.matrix)
    result <- array(0, dim = c(nspecies, nchar, nsim))
    .get.state = function(s, p) {
      pp = cumsum(p[s, ])
      min(which(runif(1) < pp))
    }
    for (j in 1:nchar) {
      m = model.matrix[[j]]
      if (!root %in% c(1:nrow(m)))
        stop(paste("'root' must be a character state from 1 to ",
                   nrow(m), sep = ""))
      p = lapply(el, function(l) expm:::expm(m * l, method="Higham08"))
      for (k in 1:nsim) {
        node.value <- numeric(nspecies + Nnode(zphy))
        node.value[rt] <- root
        for (i in nbranches:1) {
          cur = zphy$edge[i, 2]
          anc = zphy$edge[i, 1]
          curp = p[[i]]
          s = node.value[anc]
          node.value[cur] = .get.state(s, curp)
        }
        result[, j, k] <- node.value[1:nspecies]
      }
    }
    rownames(result) <- zphy$tip.label
  }
  result<-data.frame(taxa=zphy$tip.label,sample=as.numeric(result))
  return(result)
}
