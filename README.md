---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# AdaptiveDesignOptimizerSparseLP

The goal of AdaptiveDesignOptimizerSparseLP is to construct an optimal two-Stage, two subpopulation, adaptive enrichment designs for a given problem. The problem inputs are the desired familywise Type I error rate and power, and set of scenarios (data generating distributions) of interest. The software aims to minimize the expected sample size under the power and Type I error constraints.

## Installation

You can install AdaptiveDesignOptimizerSparseLP using the remotes R package, which can be obtained by typing the following in your R session:

``` r
source("https://install-github.me/r-lib/remotes")
```

and then by entering the following in your R session: 

``` r
remotes::install_github("mrosenblum/AdaptiveDesignOptimizerSparseLP")
```

## Examples and Replication of Key Results from Manuscript

The computations for the key results from the paper (Examples 3.1 and 3.2 as described in Section 5.2) can be reproduced by the code in the following 2 files in the R project's inst/examples directory: replicate.results.example.3.1.R and 
 replicate.results.example.3.2.R. We used Cplex to solve these problems, as noted in the paper.
 
Below is a simplified example that can be run in 10 minutes using the GLPK solver, which involves solving a modified version of the problem from Example 3.2 as described in Section 5.2 of the manuscript; the main modifications are that we  use a coarsened partition of the decision region and rejection regions in order to speed up the computation. 


```r
#Install R package if not already done so using the following command:
#remotes::install_github("mrosenblum/AdaptiveDesignOptimizerSparseLP")

# Load R package if not already done so by the following command: library(AdaptiveDesignOptimizerSparseLP)
# For reproducibility, set the random number generator seed:
set.seed(32515);
# Set all problem parameters based on Example 3.2, and using explicit choices of the following input parameters:
# The proportion of the population in subpopulation 1:
subpopulation.1.proportion=0.5;
# Sample size in stage 1 for each subpopulation: 50; when both subpopulations enrolled in stage 2, then 50 new participants are enrolled from each subpopulation; when only 1 subpopulation is enrolled in stage 2 then 150 new participants are enrolled from that subpopulation. This is equivalent to setting n=200 in our adaptive design template n^(1b) from Section 3.2 of the paper. These sample sizes are encoded below:
stage.1.sample.sizes=c(50,50);
stage.2.sample.sizes.per.enrollment.choice=matrix(c(50,50,
                                                   0,0,
                                                   150,0,
                                                   0,150),nrow=4,ncol=2,byrow=TRUE,dimnames=list(c(),c("Subpopulation1Stage2SampleSize","Subpopulation2Stage2SampleSize")));
# Minimum, clinically meaningful treatment effect size: Delta^min=sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5=0.465.
Delta_min = 1.2*sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5;
# The corresponding data generating distributions are encoded as follows (where we set the outcome variance to 1 for each subpopulation bby study arm combination):
data.generating.distributions=matrix(data=c(0,0,1,1,1,1,
                                           0,Delta_min,1,1,1,1,
                                           Delta_min,0,1,1,1,1,
                                           Delta_min,Delta_min,1,1,1,1),nrow=4,ncol=6,byrow=TRUE,dimnames=list(c(),c("Delta1","Delta2","Variance10","Variance11","Variance20","Variance21")));
# The resulting non-centrality parameter (see Section 5.1 of the paper) matches that used in the paper computations.
# Required Familywise Type I error:
total.alpha=0.05;
# We set the desired power for each of the power constraint to be 0.83 in the first four iterations of our algorithm, rather than the intended final output of 0.82, since we anticipate rounding in the final iteration that leads to a slight decrease in power.
desired.power=0.8;
power.constraints=matrix(c(0,0,0,
			   0,desired.power,0,
			   desired.power,0,0,
			   0,0,desired.power),nrow=4,ncol=3,byrow=TRUE,dimnames=list(c(),c("PowerH01","PowerH02","PowerH0C")));
objective.function.weights=0.25*c(1,1,1,1);
prior.covariance.matrix=diag(2);
type.of.LP.solver="glpk";
discretization.parameter=c(3,3,1);
number.cores=1;

# Run first iteration solving sparse linear program
optimized.policy <- optimize_design(subpopulation.1.proportion,total.alpha,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter,number.cores,prior.covariance.matrix=prior.covariance.matrix)
#> [1] "Adaptive Design Optimization Completed. Optimal design is stored in the file: optimized_design1.rdata"
#> [1] "Feasible Solution was Found and Optimal Expected Sample Size is 202.02"
#> [1] "Fraction of solution components with integral value solutions"
#> [1] 0.952
#> [1] "User defined power constraints (desired power); each row corresponds to a data generating distribution; each column corresponds to H01, H02, H0C desired power, respectively."
#>            Delta1 Delta2 Variance10 Variance11 Variance20 Variance21
#> Scenario 1  0.000  0.000          1          1          1          1
#> Scenario 2  0.000  0.558          1          1          1          1
#> Scenario 3  0.558  0.000          1          1          1          1
#> Scenario 4  0.558  0.558          1          1          1          1
#>            PowerH01 PowerH02 PowerH0C
#> Scenario 1      0.0      0.0      0.0
#> Scenario 2      0.0      0.8      0.0
#> Scenario 3      0.8      0.0      0.0
#> Scenario 4      0.0      0.0      0.8
#> [1] "Probability of rejecting each null hypothesis (last 3 columns) under each data generating distribution (row)"
#>            Delta1 Delta2 Variance10 Variance11 Variance20 Variance21   H01
#> Scenario 1  0.000  0.000          1          1          1          1 0.030
#> Scenario 2  0.000  0.558          1          1          1          1 0.046
#> Scenario 3  0.558  0.000          1          1          1          1 0.800
#> Scenario 4  0.558  0.558          1          1          1          1 0.624
#>              H02   H0C
#> Scenario 1 0.027 0.022
#> Scenario 2 0.800 0.236
#> Scenario 3 0.049 0.323
#> Scenario 4 0.647 0.800

print("The optimized.policy returned by our software consists of the pair of functions pi_1 and pi_2 as defined in Section 2.1 of the manuscript, along with state spaces and actions spaces for the sequential decision problem denoted by S1, A1, S2, A2. Each state in S1 is a rectangle with lower-left coordinates (x1,y1) encoded as optimized.policy$lower_boundaries and upper right coordinates (x2,y2) encoded as optimized.policy$upper_boundaries. Each action in A1 is an integer among 1,2,...,d corresponding to each of the d stage 2 enrollment choices that were input by the user. Since we allow the stage 2 state space to depend on end of stage 1 action in A1, we have that S2 is a function from an action in A1. E.g., the stage 2 state space after stage 1 action a1=3 is encoded as optimized.policy$S2[[3]]. Each action in A2 is an integer among 1,2,...,7 encoding the following outcome of the multiple testing procedure: Reject none,Reject H01,Reject H02,Reject H0C,Reject H01 and H0C,Reject H02 and H0C,Reject all, respectively.")
#> [1] "The optimized.policy returned by our software consists of the pair of functions pi_1 and pi_2 as defined in Section 2.1 of the manuscript, along with state spaces and actions spaces for the sequential decision problem denoted by S1, A1, S2, A2. Each state in S1 is a rectangle with lower-left coordinates (x1,y1) encoded as optimized.policy$lower_boundaries and upper right coordinates (x2,y2) encoded as optimized.policy$upper_boundaries. Each action in A1 is an integer among 1,2,...,d corresponding to each of the d stage 2 enrollment choices that were input by the user. Since we allow the stage 2 state space to depend on end of stage 1 action in A1, we have that S2 is a function from an action in A1. E.g., the stage 2 state space after stage 1 action a1=3 is encoded as optimized.policy$S2[[3]]. Each action in A2 is an integer among 1,2,...,7 encoding the following outcome of the multiple testing procedure: Reject none,Reject H01,Reject H02,Reject H0C,Reject H01 and H0C,Reject H02 and H0C,Reject all, respectively."

print("Lower-left and upper right coordinates of 28th rectangle in list of stage 1 states S1:")
#> [1] "Lower-left and upper right coordinates of 28th rectangle in list of stage 1 states S1:"
optimized.policy$S1[[28]]$lower_boundaries
#> [1] 1.5 1.5
optimized.policy$S1[[28]]$upper_boundaries
#> [1] 3 3

print("Stochastic policy pi_1 evaluated at 28th rectangle in S1, i.e., the probabilities of taking each action in A1 given that the first stage z-statistics are in the rectangle corresponding to this state:")
#> [1] "Stochastic policy pi_1 evaluated at 28th rectangle in S1, i.e., the probabilities of taking each action in A1 given that the first stage z-statistics are in the rectangle corresponding to this state:"
optimized.policy$pi_1(28)
#> [1] "Probabilities of enrollment decisions 1 through  4  respectively:"
#> [1] 0 1 0 0

print("Lower-left and upper right coordinates of 13th rectangle in list of stage 2 states S2 following action a1=2")
#> [1] "Lower-left and upper right coordinates of 13th rectangle in list of stage 2 states S2 following action a1=2"
optimized.policy$S2[[2]][[13]]$lower_boundaries
#> [1] 0 0
optimized.policy$S2[[2]][[13]]$upper_boundaries
#> [1] 3 3

print("Stochastic policy pi_2 evaluated at 28th rectangle in S1, action a1=2, and 13th  rectangle in S2:")
#> [1] "Stochastic policy pi_2 evaluated at 28th rectangle in S1, action a1=2, and 13th  rectangle in S2:"
optimized.policy$pi_2(28,2,13)
#> [1] "Inputs correspond to the following:"
#> [1] "s1 is rectangle defined as the Cartesian product [ 1.5 , 3 ] x [ 1.5 , 3 ]"
#> [1] "a1 is enrollment decision  2"
#> [1] "s2 is rectangle defined as the Cartesian product [ 0 , 3 ] x [ 0 , 3 ]"
#>        Reject none         Reject H01         Reject H02 
#>                  0                  0                  0 
#>         Reject H0C Reject H01 and H0C Reject H02 and H0C 
#>                  0                  0                  0 
#>         Reject all 
#>                  1
```
