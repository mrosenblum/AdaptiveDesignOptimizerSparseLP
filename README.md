
<!-- README.md is generated from README.Rmd. Please edit that file -->
AdaptiveDesignOptimizerSparseLP
===============================

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/mrosenblum/AdaptiveDesignOptimizerSparseLP.svg?branch=master)](https://travis-ci.com/mrosenblum/AdaptiveDesignOptimizerSparseLP) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/mrosenblum/AdaptiveDesignOptimizerSparseLP?branch=master&svg=true)](https://ci.appveyor.com/project/mrosenblum/AdaptiveDesignOptimizerSparseLP) <!-- badges: end -->

The goal of AdaptiveDesignOptimizerSparseLP is to optimize the decision rule and multiple testing procedure for two-stage, two subpopulation, adaptive enrichment trial designs. The methodology used is available here: https://biostats.bepress.com/jhubiostat/paper273/

Installation
------------

The software requires that you have installed a linear program solver, preferably Cplex or Gurobi since they are faster at solving large linear programs (both of which can be obtained free by academics from IBM and Gurobi, respectively), but also possible is to install Matlab or open-source GLPK. If using Cplex, you also need to have Matlab installed. For whichever linear program solving software you use, you need to ensure that its path is available to R. 
You select your solver when calling the main function (optimize_design) in this R package.

You can install the released version of AdaptiveDesignOptimizerSparseLP from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("AdaptiveDesignOptimizerSparseLP")
```

Example
-------

This is an example, which involves solving the problem from Example 3.2 as described in Section 5.2 of the manuscript. 

``` r
#Install R package if not already done so:
#remotes::install_github("mrosenblum/AdaptiveDesignOptimizerSparseLP")

# Load R package
library(AdaptiveDesignOptimizerSparseLP)

#This problem was solved in 5 iterations, with iteration based on refining the decision/rejection regions partition based on the solution from the previous iteration. Each iteration involves specifying and then solving the corresponding linear program. 

# Set all problem input parameters based on Example 3.2, but using explicit choices of sample sizes defined by setting n=200:
subpopulation.1.proportion=0.5;
total.alpha=0.05-(1e-4);
data.generating.distributions=matrix(data=c(0,0,1,1,1,1,
                                           0,sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5,1,1,1,1,
                                           sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5,0,1,1,1,1,
                                           sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5,sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5,1,1,1,1),nrow=4,ncol=6,byrow=TRUE,dimnames=list(c(),c("Delta1","Delta2","Variance10","Variance11","Variance20","Variance21")));
stage.1.sample.sizes=c(50,50);
stage.2.sample.sizes.per.enrollment.choice=matrix(c(50,50,
                                                   0,0,
                                                   150,0,
                                                   0,150),nrow=4,ncol=2,byrow=TRUE,dimnames=list(c(),c("Subpopulation1Stage2SampleSize","Subpopulation2Stage2SampleSize")));
objective.function.weights=0.25*c(1,1,1,1);
prior.covariance.matrix=diag(2);
power.constraints=matrix(c(0,0,0,
			   0,0.83,0,
			   0.83,0,0,
			   0,0,0.83),nrow=4,ncol=3,byrow=TRUE,dimnames=list(c(),c("PowerH01","PowerH02","PowerH0C")));
type.of.LP.solver="matlab";
discretization.parameter=c(1,1,10);
number.cores=30;

# Run first iteration solving sparse linear program at slightly stronger power constraint (0.83) than the intended final output (0.82)
###optimize_design(subpopulation.1.proportion,total.alpha,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter,number.cores,ncp.list=c(),list.of.rectangles.dec=c(),LP.iteration=1,prior.covariance.matrix)
# Construct new list of familywise Type I error constraints focusing on neighborhoods of active constraints from previous iteration 
# The Active Type I error constraints are contained in the output file optimized.design1.rdata and stored in the variable ncp.active.FWER.constraints, and in this case are the following 3 pairs: 0 0; 0.000000 1.960894; 1.960894 0.000000
ncp.list <- list(c(0,0))
for(z in seq(1.5,3,length=30)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)))}
for(z in seq(-9,9,length=30)) {ncp.list <- c(ncp.list,list(c(z,-z)))}

# Run second iteration solving sparse linear program as in first iteration but using Type I error constraints defined above
###optimize_design(subpopulation.1.proportion,total.alpha,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter,number.cores,ncp.list,list.of.rectangles.dec=c(),LP.iteration=2,prior.covariance.matrix)
# Refine Decision Region Rectangle Set
load("optimized.design2.rdata")
list.of.rectangles.dec <- refine_decision_rectangles(subpopulation.1.proportion=0.5,stage.1.sample.sizes=stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice=stage.2.sample.sizes.per.enrollment.choice,list.of.rectangles.dec=list.of.rectangles.dec,sln=sln,set.rectangles.with.identically.valued.neighbors.and.split.others=TRUE,LP.iteration=2)
# Construct new list of familywise Type I error constraints focusing on neighborhoods of active constraints from previous iteration 
# The Active Type I error constraints are contained in the output file optimized.design2.rdata and stored in the variable ncp.active.FWER.constraints
ncp.list <- list(c(0,0))
for(z in seq(1.5,3,length=30)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)))}
for(z in seq(-9,9,length=10)) {ncp.list <- c(ncp.list,list(c(z,-z)))}

# Run third iteration solving sparse linear program as in first iteration but using Type I error constraints defined above
###optimize_design(subpopulation.1.proportion,total.alpha,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter,number.cores,ncp.list,list.of.rectangles.dec,LP.iteration=3,prior.covariance.matrix)
# Refine Decision Region Rectangle Set
load("optimized.design3.rdata")
list.of.rectangles.dec <- refine_decision_rectangles(subpopulation.1.proportion=0.5,stage.1.sample.sizes=stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice=stage.2.sample.sizes.per.enrollment.choice,list.of.rectangles.dec=list.of.rectangles.dec,sln=sln,set.rectangles.with.identically.valued.neighbors.and.split.others=TRUE,LP.iteration=3)
# Construct new list of familywise Type I error constraints focusing on neighborhoods of active constraints from previous iteration 
# The Active Type I error constraints are contained in the output file optimized.design3.rdata and stored in the variable ncp.active.FWER.constraints
ncp.list <- list(c(0,0))
for(z in seq(1.8,2.5,length=20)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)))}

# Run fourth iteration solving sparse linear program as in first iteration but using Type I error constraints defined above
###optimize_design(subpopulation.1.proportion,total.alpha,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter,number.cores,ncp.list,list.of.rectangles.dec,LP.iteration=4,prior.covariance.matrix)
# Set Decision Region Rectangles 
load("optimized.design4.rdata")
list.of.rectangles.dec <- refine_decision_rectangles(subpopulation.1.proportion=0.5,stage.1.sample.sizes=stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice=stage.2.sample.sizes.per.enrollment.choice,list.of.rectangles.dec=list.of.rectangles.dec,sln=sln,set.rectangles.with.identically.valued.neighbors.and.split.others=FALSE,round.each.decision.rectangle.to.integer=TRUE,LP.iteration=4)
# Construct new list of familywise Type I error constraints focusing on neighborhoods of active constraints from previous iteration 
# The Active Type I error constraints are contained in the output file optimized.design4.rdata and stored in the variable ncp.active.FWER.constraints
ncp.list <- list(c(0,0))
for(z in seq(2,2.5,length=500)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)))}
for(z in seq(-9,9,length=950)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)),list(c(z,-z)))}
for(z in c(seq(-3,-1,length=500),seq(1,3,length=500))) {ncp.list <- c(ncp.list,list(c(z,-z)))}

# Run fifth iteration solving sparse linear program as in first iteration but using Type I error constraints defined above
optimize_multiple_testing_procedure(subpopulation.1.proportion,total.alpha=0.049,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter=c(1,0.25,10),number.cores,ncp.list,list.of.rectangles.dec,LP.iteration=5,prior.covariance.matrix,round.each.multiple.testing.procedure.rectangle.to.integer=TRUE)
```
