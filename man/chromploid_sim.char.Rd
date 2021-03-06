% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chromploid_sim.char.R
\name{chromploid_sim.char}
\alias{chromploid_sim.char}
\title{Simulates a a sample of chromosome number or ploidy change for models with large sparse matrices}
\usage{
For BiChroM simulation sim.char(phy=mytree, par=Q_bichrom, nsim=1, model="discrete", root=56)
}
\arguments{
\item{phy}{Bifurcating phylogenetic tree}

\item{par}{Q-matrix model to be used}

\item{nsim}{number of simulations by default nsim=1}

\item{model}{by default model="discrete" for other usage see geiger:::sim.char}

\item{root}{discrete value at the root, by default 1}
}
\description{
Simulates a a sample of chromosome number or ploidy change for models with large sparse matrices
}
