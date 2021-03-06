chromploid Package
================

Overview of Chromploid
----------------------

Chromploid is a package that calculates the negative log-likelihood function forrates of chromosome number or ploidy change models in phylogenetic trees. Currently, chromploid includes four different models of chromosome number or ploidy evolution: BiChroM, BiChroM linear (beta), PloidEvol, and ChromEvol M3. Users can choose among these three models but if interested Chromploid package can be extended to accomodate custom chromosome number or ploidy change models.

About the models included

1.  BiChroM: A continuous time Markov chain that links the evolution of chromosome number change to a binary trait. This model contains 10 parameters that represent rates (change per time).[1]
2.  BiChroM linear: BETA VERSION. A continuous time Markov chain that links the evolution of chromosome number change to a binary trait where the rates of binary trait change are linear functions of chromosome number. This model is define via 12 parameters that represent rates of evolution (change per time).
3.  PloidEvol: A continuous time Markov chain that models the evolution of ploidy change. This model contains 6 parameters that represent rates of ploidy evolution (change per time)[2].
4.  ChromEvol M3: A continuous time Markov chain that models the evolution of chromosome number change. The model has 4 parameters that represent rates of chromosome number change. This model was proposed by Mayrose et al. (2010)[3]

Installation
------------

You'll need package `devtools`. Chromploid is easily installed from github

``` r
library("devtools")
devtools::install_github("roszenil/chromploid")
library("chromploid")
library("geiger") 
library("nloptr")
#chromploid depends on geiger and nloptr so don't forget to load the packages
```

Negative log likelihood
-----------------------

The function `chromploid_nllike()` is possibly the most important and critical of the package. It calculates the negative log-likelihood of any discrete trait model for chromosome number or ploidy change. Minimizing negative log-likelihoods allows us to infer model parameters, that is, obtaining maximum likelihood estimates (MLEs), profile likelihoods, likelihood-confidence intervals, and likelihood ratio tests necessary to answer biological questions. I will describe how to use this important function after introducing the models coded in chromploid.

BiChroM
-------

BiChroM -Binary trait associated with chromosome number change model- is a continuous time Markov chain (CTMC) that allows for estimation and testing of chromosomal change rates correlated with a binary trait. The model is defined via 10 parameters that represent rates of chromosomal change linked to a value of a binary trait. BiChroM allows for testing the association of chromosomal change (polyploidy) to phenotypic change.

Description of parameters in BiChroM:

-   *λ*<sub>0</sub>, *λ*<sub>1</sub>: Rate of gain in one chromosome for taxa with binary trait 0, or binary trait 1 respectively.
-   *μ*<sub>0</sub>, *μ*<sub>1</sub>: Rate of loss in one chromosome for taxa with binary trait 0, or binary trait 1 respectively.
-   *ρ*<sub>0</sub>, *ρ*<sub>1</sub>: Rate of chromosome doubling for taxa with binary trait 0, or binary trait 1 respectively.
-   *q*<sub>01</sub>, *q*<sub>10</sub>: Rate of transition between binary state values, from 0 to 1 or vice versa.
-   *ϵ*<sub>0</sub>, *ϵ*<sub>1</sub>: Ancillary parameters of chromosomal changes. Rate of chromosome number changes after `size` has been determined.

### Full inference for BiChroM model

Inferences using BiChoM model require different functions that allow the calculation of full and partial likelihoods for the parameters. For example, when the goal is to estimate rates, the full likelihood of 10 parameters should be calculated and optimized, but if the goal of the inference is to assess the uncertainty about a single parameter then a univariate likelihood should be the calculated.

Each of the functions listed below creates a Q-matrix, a matrix that defines the rates of chromosomal change linked to the binary trait. The difference among these functions is that each one builds the necessary mathematics to create a different likelihood function, thus, allowing for different types of inference.

-   `Q_bichrom()`: Creates the Q-matrix (instantaneous probability rate matrix) for the full model, that allows for calculation of maximum likelihood estimates.
-   `Q_unichrom()`: Creates the marginal Q-matrix necessary for calculation of univariate profile likelihoods. Used to assess the uncertainty surrounding an individual parameter.
-   `Q_bibichrom()`: Creates the marginal Q-matrix necessary for calculation of bivariate profile likelihoods.Used to assess the uncertainty of two parameters at a time, allows to detect problems with parameter identifiability.[4]
-   `Q_reducedbichrom()`: Creates the marginal Q-matrix necessary for calculation of likelihood ratio tests.
-   `bichrom_dataset()`: Transforms a data set of chromosome and binary traits to input in the functions above

1.  BiChroM simulation and estimation

-   Choose a set of fixed parameter values. For example (*λ*<sub>0</sub>, *λ*<sub>1</sub>, *μ*<sub>0</sub>, *μ*<sub>1</sub>, *ρ*<sub>0</sub>, *ρ*<sub>1</sub>, *q*<sub>01</sub>, *q*<sub>10</sub>, *ϵ*<sub>0</sub>, *ϵ*<sub>1</sub>)= (0.12, 0.001, 0.25, 0.002, 0.036, 0.006, 0.04, 0.02, 1.792317852, 1.57*e* − 14)

-   Calculate BiChroM Q-matrix All parameters in BiChroM are used in log scale because numerical calculations are easier than in original scale.

``` r
log.params <- log(c(0.12, 0.001, 0.25, 0.002, 0.036, 0.006, 0.04, 
    0.02, 1.792317852, 1.57e-14))
N <- 50
mymatrix <- Q_bichrom(log.theta = log.params, size = N)
```

-   Simulate a tree with 500 taxa and BiChroM values Assuming that the root of the tree has a value 56 means that the root is has 5 haploid chromosomes and its binary type is 1.

``` r
mytree<- sim.bdtree(b=0.055, stop="taxa", n=500, seed=1) 
# Seed was fixed to obtain always the same tree for this tutorial
plot(mytree,type="fan",cex=0.2, no.margin=TRUE)
```

![](chromploidmd_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-3-1.png)

``` r
set.seed(1) #Seed was fixed to obtain always the same sample for this tutorial
mysample<- chromploid_sim.char(phy=mytree, par=mymatrix,model="discrete", root=56)
```

-   Calculate the likelihood of the tree, assuming that the root has a uniform distribution for chromosome number and binary state `pi.0=NULL` (other options for the root are available, like Maddison and Fitzjohn (2009)). The value of the negative log-likelihood `nllike` in this example is 943.9293, since we fixed the random seed, and it is calculated as follows.

``` r
nllike <- chromploid_nllike(log.par = log.params, phy = mytree, 
    tip.values = mysample, pi.0 = NULL, Q.FUN = Q_bichrom, Q.ARGS = list(size = 50))
nllike
```

    ## [1] 943.9293

If you have the following error appearing

``` r
# Error in tip.values[i, charnum + 1] : incorrect number of
# dimensions
```

this indicates that you don't have in the first column the species names and in the second column the sample for BiChrom. An easy way to fix it is converting your sample into a data frame with first column having taxa names and second column having the sample needed in BiChroM just as `chromploid_sim.char` does.

``` r
mysample <- data.frame(taxa = rownames(mysample), sample = as.numeric(mysample))
head(mysample)
```

1.  BiChrom optimization Usually you don't know the values of the parameters and what you need is toto find the maximum likelihood estimates of the parameters determined by BiChroM. In this case, you can optimize (minimize) the negative log-likelihood with any R optimizer given that you have a phylogenetic tree and a sample of chromosome numbers and binary trait using the function `chromploid_nllike`. I have written a shortcut in the optimization using function `chromploid_optim` but if you think that the likelihood surface is complicated I suggest using your own customiza optimizer with the negative-log likelihood function. In addition, loading the package `nloptr` will help time consuming calculations because it makes very accurate searches in the parametric space. The cost of using `nloptr` is computational time.

``` r
# Optimizations take hours to run, I suggest you run this
# routine in a cluster. For a 4700 tip tree starting with a
# x.0 close to the maximum it can take up to 96 hrs.

library(nloptr)
x.0 <- log(c(0.12, 0.001, 0.25, 0.002, 0.036, 0.006, 0.04, 0.02, 
    1.792317852, 1.57e-14))
# value where the optimization algorithm is going to start. A
# vector of 10 values for the parameters in log scale

model.args <- list(size = 50)
optimization.bichrom <- chromploid_optim(phy.tree = mytree, tip.values = mysample, 
    model = "bichrom", model.args = model.args, starting.val = x.0, 
    root.dist = NULL)

print(optimization.bichrom)
```

The vector `optimization.bicrhom` contains in the first 10 entries maximum likelihood estimates of parameters (*λ*<sub>0</sub>, *λ*<sub>1</sub>, *μ*<sub>0</sub>, *μ*<sub>1</sub>, *ρ*<sub>0</sub>, *ρ*<sub>1</sub>, *q*<sub>01</sub>, *q*<sub>10</sub>, *ϵ*<sub>0</sub>, *ϵ*<sub>1</sub>). Then the negative log likelihood value at the MLE, and finally the status of the convergence to the optimum inherited from package nloptr. Do not expect good estimates for *ϵ*<sub>0</sub> and *ϵ*<sub>1</sub> since by definition these are ancillary parameters that prevent biases for the other 8 that are key to understand chromosomal change (see Zenil-Ferguson et al. 2017 for details). The eleventh entry in `optimization.bicrhom` is the value of the negative log-likelihood at the maximum likelihood estimates and it is a useful when you want to plot relative profile likelihoods or calculate likelihood ratio tests.

1.  Likelihood Ratio Test To calculate the likelihood of the reduced model needed for likelihood ratio tests you can use the function `chromploid_nllike` with Q-matrix `Q_reducedbichrom` or use the wrapper optimization function `chromploid_optim`. The example here tests if the rates associated with chromosome doubling under the different values of the binary trait are equal that is *H*<sub>0</sub> : *ρ* = *ρ*<sub>0</sub> = *ρ*<sub>1</sub>

``` r
x.0 <- log(c(0.12, 0.001, 0.25, 0.002, 0.01, 0.04, 0.02, 1.792317852, 
    1.57e-14))
# Value where the optimization algorithm is going to start. A
# vector of 9 values for the parameters in log scale. The
# value of the hypothesis here is rho=0.010

model.args = list(size = 50, equal.param = c("rho0", "rho1"), 
    location.par = 5)
# Q.ARGS is a list that has all the arguments included in
# Q_reducedbichrom that are not the parameters size= maximum
# number of haploid chromosome numbers to consider in your
# sample, equal.params=which parameters are equal based on
# the hypothesis H0, location.par is the position in the
# vector where rho value appears

optimization.reducedbichrom <- chromploid_optim(phy.tree = mytree, 
    tip.values = mysample, model = "reducedbichrom", model.args = model.args, 
    starting.val = x.0, root.dist = NULL)
```

The vector `optimization.reducedbicrhom` contains in the first 9 entries maximum likelihood estimates of parameters (*λ*<sub>0</sub>, *λ*<sub>1</sub>, *μ*<sub>0</sub>, *μ*<sub>1</sub>, *ρ*, *q*<sub>01</sub>, *q*<sub>10</sub>, *ϵ*<sub>0</sub>, *ϵ*<sub>1</sub>). Then the negative log likelihood value at the MLE, and finally the status of the convergence to the optimum inherited from package nloptr.

print(optimization.reducedbichrom) \`\`\`

Finally the value of the likelihood ratio test for *H*<sub>0</sub> : *ρ*<sub>0</sub> = *ρ*<sub>1</sub> hypothesis is simply *D* = 2 × (reduced bichrom nloglike − full model nloglike) that has a *χ*<sub>(1)</sub><sup>2</sup> distribution

``` r
neglog.red <- optimization.reducedbichrom[10]
neglog.full <- optimization.bichrom[11]
alpha <- 0.05
D = 2 * (neglog.red - neglog.full)
p.value <- pchisq(D, lower.tail = FALSE, df = 1)
if (p.value > 0.05) {
    reject <- 0
} else {
    reject <- 1
}
```

1.  Profile likelihoods (this will change soon when `chromploid_optim` is fully defined, meanwhile the process here works)

Point estimations are useful but so are interval estimations. Profile likelihoods allow for better understanding of the precision, estimability, and identifiability of estimates, given BiChroM model, the sample size, and the phylogenetic tree.

-   Univariate profile likelihoods To understand the behavior of one parameter of a time under BiChroM model the function you need is `Q_unibichrom`. Assume that you want to see the profile likelihood of parameter *ρ*<sub>0</sub>. You have to choose a grid (range of values) for which the likelihood could be interesting. I suggest again that you do this procedure in a cluster. In the example here shown, I do a sequential way of calculating the profile in this example but more efficiently you could program this in parallel, sending one optimization per processor, retrieving all the results in 24hrs, so you don't have to wait for months. I also suggest using as your starting point x.0 the maximum likelihood value (in log-scale) that you obtained from optimizing `Q_bichrom`, optimizations will take shorter if your starting point is good.

``` r
logrho0.values <- log(seq(0.02, 0.04, 0.005))
# Calculate the likelihood values in a grid 0.02, 0.025, ...,
# 0.035, 0.04 (in log scale). True value is around 0.36 so
# values surrounding the true value are interesting
long <- length(logrho0.values)
x.0 <- log(c(0.12, 0.001, 0.25, 0.002, 0.006, 0.04, 0.02, 1.792317852, 
    1.57e-14))
# value where the optimization algorithm is going to start. A
# vector of 9 values for the nuisance parameters in log scale
# (everything except rho0)
profile.values <- rep(0, long)  # Vector useful to save optimization results
my.options <- list(algorithm = "NLOPT_LN_SBPLX", ftol_rel = 1e-08, 
    print_level = 1, maxtime = 1.7e+08, maxeval = 1000)
# Options used in nloptr, for more information ?nloptr. This
# is a subplex algorithm, with high tolerance to fine a very
# precise optimum.
for (i in 1:long) {
    mle.profilerho0 <- nloptr(x0 = x.0, eval_f = chromploid_nllike, 
        opts = my.options, phy = mytree, tip.values = mysample, 
        pi.0 = NULL, Q.FUN = Q_unibichrom, Q.ARGS = list(log.theta0 = logrho0.values[i], 
            size = 50, param.name = "rho0"))
    # Q.ARGS is a list that has all the arguments included in
    # Q_unibichrom Type ?Q_unibichrom to see a longer explanation
    # of each argument
    print(mle.profilerho0)
    profile.values[i] <- mle.profilerho0$objective
}

plot(exp(logrho0.values), profile.values, xlab = expression(rho_0), 
    ylab = "Negative log-likelihood", type = "l", lwd = 2)
```

-   Bivariate profile likelihood

Analogous to the example above but using a grid for two parameters and calculating using the function `Q_bibichrom`.

BiChroM Linear
--------------

THIS IS A BETA VERSION BiChroM linear assumes the same chromosome number changes than BiChroM but the binary trait change is a linear function of chromosome number. This model can be useful when trying to associate if polyploidy is associated with the speed of binary trait change. The model is defined via 12 parameters, chromosomal change rates *λ*<sub>0</sub>, *λ*<sub>1</sub>, *μ*<sub>0</sub>, *μ*<sub>1</sub>, *ρ*<sub>0</sub>, *ρ*<sub>1</sub> are exactly the same as in the simple BiChroM but binary trait change rates are defined as follows.

-   *κ*<sub>0</sub> + *i**q*<sub>01</sub>: Chromosome number *i* dependent transition rate from state 0 to state 1.
-   *κ*<sub>1</sub> + *i**q*<sub>10</sub>: Chromosome number *i* dependent transition rate from state 1 to state 0.

### Usage

Substitute `Q_bichrom()` for `Q_bichromlinear()` but the arguments for this new function are a vector log.par of size 12, with parameters *λ*<sub>0</sub>, *λ*<sub>1</sub>, *μ*<sub>0</sub>, *μ*<sub>1</sub>, *ρ*<sub>0</sub>, *ρ*<sub>1</sub>, *κ*<sub>0</sub>, *κ*<sub>1</sub>, *q*<sub>01</sub>, *q*<sub>10</sub>, *ϵ*<sub>0</sub>, *ϵ*<sub>1</sub> in this order and in log-scale, and `size` the maximum number of haploid chromosome numbers allowed.

BiChroM example of power analysis
=================================

Libraries needed for the simulation example

``` r
library("geiger")
library("chromploid")
```

Equal chromosome number doubling rate in BiChroM simulation
-----------------------------------------------------------

First choose some parameter values in log scale for the 10 parameters of BiChroM and the maximum number of haploid chromosomes to build the Q-matrix.

``` r
# l.0, l.1, m.0, m.1, r.0, r.1, prob.01, prob.10, e.0, e.1
# input in log-scale
params <- log(c(0.01, 0.005, 0.01, 0.005, 0.01, 0.01, 0.01, 0.005, 
    1e-06, 1e-06))
# Maximum number of haploid chromosomes
N <- 25
```

Then, simulate a phylogeny with 250 tips and BiChroM process on top of it

``` r
# simulating a phylogeny with 250 tips (birth process)
tree <- sim.bdtree(b = 0.03, stop = "taxa", n = 250)  # Approximatedly 110 million years

# Bichrom.model defined via Q-matrix
bichrom.model <- Q_bichrom(log.theta = params, size = N)
# Simulating chromosome number and binary trait change
# states, root 27 reflects that the simulation starts at the
# root having binary trait 1 with 1 haploid chromosome number
res <- chromploid_sim.char(phy = tree, par = bichrom.model, nsim = 1, 
    root = 27)
```

Now you have some data.

Optimization of full and reduced models of BiChoM
-------------------------------------------------

The optimization process requires you to assign a probability distribution at the root of the tree. In this example is uniform (we really don't know the chromosome number or binary trait at the root). Then, you can use a wrapper function like `bichrom_optim` to optimize full and reduced BiChroM models with the data you just generated. If you are interested in optimizing your own model you can just follow what is inside this function but with your own custom Q-matrix (see solanum example).

``` r
## Uniform distribution of values at the root
p.0 <- rep(1, 2 * (N + 1))/(2 * (N + 1))
model.args <- list(size = 25)
optimization.bichrom <- bichrom_optim(phy.tree = tree, tip.values = res, 
    model = "bichrom", model.args = model.args, starting.val = params, 
    root.dist = p.0)

model.args = list(size = 50, equal.param = c("rho0", "rho1"), 
    location.par = 5)

optimization.reducedbichrom <- bichrom_optim(phy.tree = tree, 
    tip.values = res, model = "reducedbichrom", model.args = model.args, 
    startin.val = params[-6], root.dist = p.0)
```

Likelihood ratio test
---------------------

Since you already calculated the negative log-likelihoods using `bichrom_optim` function you only need to calculate statistic `D` and the p-value as follows.

``` r
neglog.red <- optimization.reducedbichrom$nloglike
neglog.full <- optimization.bichrom$nloglike
alpha <- 0.05
D = 2 * (neglog.red - neglog.full)
p.value <- pchisq(D, lower.tail = FALSE, df = 1)
if (p.value < 0.05) {
    reject <- 0
} else {
    reject <- 1
}
```

Example of a custom chromosome number model in Solanum
======================================================

Libraries needed for the Solanum example

``` r
library("geiger")
library("nloptr")
library("chromploid")
```

Definition of full and reduced models via Q-matrix
--------------------------------------------------

The full model of Solanum is defined via a Q-matrix with six parameters. Two polyploidy rates (*ρ*<sub>*H*</sub>, *ρ*<sub>*W*</sub>), two demiploidy rates (*ϵ*<sub>*H*</sub>, *ϵ*<sub>*W*</sub>), and two transition rates between binary traits (*q*<sub>*H**W*</sub>, *q*<sub>*W**H*</sub>). The model has 5 states for herbaceous and 5 states for woody because chromosome numbers in the sample aree 12, 18, 24, 36, and 48. The Q-matrix is simply defined as follows:

``` r
Q_solanum <- function(size = 5, log.theta) {
    theta <- exp(log.theta)
    # parameters are in alphabetical order, it holds true for all
    # chromploid functions
    epsilon.0 <- theta[1]
    epsilon.1 <- theta[2]
    rho.0 <- theta[3]
    rho.1 <- theta[4]
    prob.01 <- theta[5]
    prob.10 <- theta[6]
    C = size
    # pure duplications
    Q <- matrix(rep(0, 4 * C * C), nrow = 2 * size)
    for (i in 1:3) {
        Q[i, (i + 2)] <- rho.0
        Q[(5 + i), (i + 7)] <- rho.1
    }
    # demiploidies
    for (i in 1:4) {
        Q[i, (i + 1)] <- epsilon.0
        Q[(5 + i), (i + 6)] <- epsilon.1
    }
    # transitions between binary state
    for (i in 1:5) {
        Q[i, (i + 5)] <- prob.01
        Q[(5 + i), i] <- prob.10
    }
    # diagonal of Q-matrix sums to zero
    diagQ <- apply(Q, 1, sum)
    for (i in 1:10) {
        Q[i, i] <- -diagQ[i]
    }
    return(Q)
}
```

To defined the Q-matrix for the reduced model where *ρ* = *ρ*<sub>*H*</sub> = *ρ*<sub>*W*</sub> we do exactly the same as before but we equal the vector entry number 3 to the two parameters, that is

``` r
Q_solanumred <- function(size = 5, log.theta) {
    theta <- exp(log.theta)
    epsilon.0 <- theta[1]
    epsilon.1 <- theta[2]
    rho.0 <- theta[3]
    # these two rhos have the same entry of theta which is size 5
    # instead of 6
    rho.1 <- theta[3]
    prob.01 <- theta[4]
    prob.10 <- theta[5]
    C = size
    Q <- matrix(rep(0, 4 * C * C), nrow = 2 * size)
    for (i in 1:3) {
        Q[i, (i + 2)] <- rho.0
        Q[(5 + i), (i + 7)] <- rho.1
    }
    for (i in 1:4) {
        Q[i, (i + 1)] <- epsilon.0
        Q[(5 + i), (i + 6)] <- epsilon.1
    }
    for (i in 1:5) {
        Q[i, (i + 5)] <- prob.01
        Q[(5 + i), i] <- prob.10
    }
    diagQ <- apply(Q, 1, sum)
    for (i in 1:10) {
        Q[i, i] <- -diagQ[i]
    }
    return(Q)
}
```

Those are the only two definitions needed for chromploid to perform a likelihood ratio test.

Optimization of the full and reduced models
-------------------------------------------

The optimization process requires you to assign a probability distribution at the root of the tree. In this example is uniform (we really don't know the chromosome number or binary trait at the root). We also use some starting values in the optimization called `x.0`.

Here are the steps to optimize the full model

``` r
mytree<-read.tree("solanumwoody2.tre") 
mysample<-read.csv("solanumtransformed.csv", header=FALSE)
#this contains the dataset transformed, the orginal dataset is in solanumsimpledataset.csv 
#5 states for woody and 5 for herbaceous coded 0=Herbaceous 1=woody 
last.state=5 

#Uniform distribution for root
p.0<-rep(1,2*(last.state))/(2*(last.state)) 

#Initial parameter values for optimization (random)
x.0<-log(c(0.182742343,0.214736895,0.787883643,0.550296163,0.553774817,0.368289079)) 

#store results
resultsfull<-rep(0,8) 
# Options used in nloptr, for more information ?nloptr. This
# is a subplex algorithm, with high tolerance to fine a very
# precise optimum.
my.options<-
list("algorithm"="NLOPT_LN_SBPLX","ftol_rel"=1e-08,"print_level"=1,
     "maxtime"=170000000, "maxeval"=5000) 
#Full model optimization
full.model<- nloptr(x0 = x.0, eval_f = chromploid_nllike, opts = my.options,
    phy = mytree, tip.values = mysample, pi.0 =p.0, Q.FUN = Q_solanum,
Q.ARGS = list(size = 5)) 
print(full.model)
resultsfull[1:6]<-mle$solution
        resultsfull[7]<-mle$objective
        resultsfull[8]<-mle$status}
        resultsfull<-as.dataframe(resultsfull)
        names(results)<-
          c("epsilon0","epsilon1","rho0","rho1","q01","q10","nloglike","convergencestatus")
```

Now the reduced model optimization works exactly the same but replacing `x.0` for a vector of 5 parameter (because the polyploidy rates are equal), and `Q_solanum` for `Q_solanum_red`.

``` r
#Initial parameter values for optimization (random)
x.0<-log(c(0.626796332,0.889894765,0.78642404,0.786166935,0.634590639)) 

#store results
resultsred<-rep(0,7) 
# Options used in nloptr, for more information ?nloptr. This
# is a subplex algorithm, with high tolerance to fine a very
# precise optimum.
my.options<-
  list("algorithm"="NLOPT_LN_SBPLX","ftol_rel"=1e-08,"print_level"=1,
       "maxtime"=170000000, "maxeval"=5000) 
#Full model optimization
reduced.model<- nloptr(x0 = x.0, eval_f = chromploid_nllike, opts = my.options,
    phy = mytree, tip.values = mysample, pi.0 =p.0, Q.FUN = Q_solanum,
Q.ARGS = list(size = 5)) 
print(reduced.model)
resultsred[1:5]<-mle$solution
        resultsred[6]<-mle$objective
        resultsred[7]<-mle$status}
        resultsred<-as.dataframe(resultsred)
        names(resultsred)<-
          c("epsilon0","epsilon1","rho","q01","q10","nloglike","convergencestatus")
```

Likelihood ratio test
---------------------

Since you already calculated the negative log-likelihoods using `nloptr` in `chromploid_nllike` function you only need to calculate statistic `D` and the p-value as follows.

``` r
neglog.red <- resultsred$nloglike
neglog.full <- resultsfull$nloglike
alpha <- 0.05
D = 2 * (neglog.red - neglog.full)
p.value <- pchisq(D, lower.tail = FALSE, df = 1)
if (p.value < 0.05) {
    reject <- 0
} else {
    reject <- 1
}
```

Some more details
-----------------

The dataset has chromosome numbers ranging from 12 to 48. To make this models work and not having to specify really large matrices with lots of zeros it is important to recode the original dataset into 10 states as follows:

| Herbaceous Csomes | Woody Csomes |
|-------------------|--------------|
| 12=1              | 12=6         |
| 18=2              | 18=7         |
| 24=3              | 24=8         |
| 36=4              | 18=9         |
| 36=5              | 18=10        |

And those are the 10 states used by `Q_solanum` and what you would find in the transformed dataset.

[1] Zenil‐Ferguson, R., Ponciano, J.M. and Burleigh, J.G., 2017. Testing the association of phenotypes with polyploidy: An example using herbaceous and woody eudicots. Evolution. 71(5):1138-1148.

[2] Zenil-Ferguson et al. in prep

[3] Mayrose, I., Barker, M.S. and Otto, S.P., 2010. Probabilistic models of chromosome number evolution and the inference of polyploidy. Systematic Biology, 59(2), pp.132-144.

[4] Ponciano, J.M., Burleigh, J.G., Braun, E.L. and Taper, M.L., 2012. Assessing parameter identifiability in phylogenetic models using data cloning. Systematic biology, p.sys055.
