---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# AdaptiveDesignOptimizerSparseLP

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/mrosenblum/AdaptiveDesignOptimizerSparseLP.svg?branch=master)](https://travis-ci.com/mrosenblum/AdaptiveDesignOptimizerSparseLP)
<!-- badges: end -->

The goal of the R package AdaptiveDesignOptimizerSparseLP is to construct an optimal two-stage, two subpopulation, adaptive enrichment design for a given problem. The problem inputs are the desired familywise Type I error rate and power, and the set of scenarios (data generating distributions) of interest. The software minimizes the expected sample size under the power and Type I error constraints.
We assume that the user has read our manuscript describing the general method (currently under review).

The user chooses the following (which are inputs to the software): the subpopulation 1 proportion, the power constraints, the familywise Type I error rate, the prior distribution used to define the objective function (which can be any finite mixture of point masses or bivariate normal distributions), and the adaptive design template (consisting of the stage 1 sample sizes and the stage 2 sample sizes under each possible stage 2 enrollment choice). The optimization goal is to minimize  expected sample size under the power and Type I error constraints.
     The software optimizes over the discretized policies defined in Section 4.2 to produce an optimized, two-stage adaptive trial design tailored to the user's inputs. 
     
  The R package calls a linear program solver and is compatible with the solvers in Matlab, Cplex, and  Gurobi (three commercial solvers, with the latter two free for academic use) and also with the open-source GLPK solver. We wrote scripts that  make it seamless to use our package with any of these solvers, though we recommend Cplex or Gurobi due to their high performance.
  Our software reproduces our examples as well. 

  If you have any questions or want help using the software, please email Michael Rosenblum at mrosen at jhu dot edu.
    
## R Package Installation

You can install AdaptiveDesignOptimizerSparseLP using the remotes R package, which can be obtained by typing the following in your R session:

``` r
source("https://install-github.me/r-lib/remotes")
```

and then by entering the following in your R session: 

``` r
remotes::install_github("mrosenblum/AdaptiveDesignOptimizerSparseLP")
```

## Linear Program Solver Installation:

You need to have installed a linear program solver to use our software.
Our software is compatible with Cplex, Gurobi, Matlab, and GLPK solvers.

We recommend using Cplex or Gurobi whenever possible (due to their improved speed
compared to the GLPK and Matlab solvers). 
Cplex and Gurobi are commercial solvers, but are free for academic researchers. 
See the following webpages to obtain these, respectively:

Cplex: https://www.ibm.com/support/pages/ibm-ilog-optimization-academic-initiative

Gurobi: https://www.gurobi.com/academia/academic-program-and-licenses/

The Matlab solver comes with Matlab: https://www.mathworks.com/products/matlab.html

GLPK can be downloaded free here: https://www.gnu.org/software/glpk/ 

Important Notes: 

GLPK needs to be installed after it is downloaded, using instructions at the link above.

Our software calls Cplex through a Matlab interface, so having Matlab is also required if one plans to use the Cplex linear program solver; the other linear program solvers (in GLPK, Gurobi, Matlab) do not have this constraint, i.e., each solver can be used without requiring any of the other software. 

If using Gurobi, the gurobi R package is required; this can be obtained and installed by following the instructions from Gurobi here: https://www.gurobi.com/documentation/8.1/quickstart_mac/r_installing_the_r_package.html

## Examples and Replication of Key Results from Manuscript

The computations for the key results from the paper (Examples 3.1 and 3.2 as described in Section 5.2) can be reproduced by the code in the following 2 files in the R project's inst/examples directory: replicate.results.example.3.1.R and 
 replicate.results.example.3.2.R. We used Cplex to solve these problems, as noted in the paper.
 
Below is a simplified example that can be run in 4 minutes using a 4 core, 2.8 GHz processor on a Macbook laptop using the GLPK solver, which involves solving a modified version of the problem from Example 3.2 as described in Section 5.2 of the manuscript; the main modifications are that we use a coarsened partition of the decision region and rejection regions in order to speed up the computation for illustration purposes. 

To obtain definitions of each input and output argument in our main function optimize_design,
type help(optimize_design) after installing the R package. 


```r
# Install R package if not already done so using the following command:
# remotes::install_github("mrosenblum/AdaptiveDesignOptimizerSparseLP")
# Load R package:
library(AdaptiveDesignOptimizerSparseLP)
# For reproducibility, set the random number generator seed:
set.seed(32515)

# Problem Inputs:
#  For illustration below, we set all problem inputs based on Example 3.2.
#  However, our software is general in that users can set the 
#   inputs based on their own problems, and then run our trial design 
#   optimizer.

#  The proportion of the population in subpopulation 1 is 1/2:
subpopulation.1.proportion = 0.5

# Sample Sizes:
#  We set n=200 in our adaptive design template n^(1b). 
#  This corresponds to total sample size 100 in stage 1, 
#  and total sample size 100 in stage 2 if both subpopulations get enrolled
#  in stage 2.

#  Stage 1 Sample Size for each subpopulation:
#  Since p1=1/2 and total sample size in stage 1 is 100,
#  we set the stage 1 sample size in each subpopulation to be 50 participants.
#  This is encoded as follows (where first entry is subpopulation 1 sample size,
# second entry is subpopulation 2 sample size):
stage.1.sample.sizes = c(50, 50)

# In general, one can specify any finite number of sample size options for 
#  stage 2 enrollment. However, at some point the 
#  resulting problem size will be too large to solve, depending on the 
#  computational resources one has available.
#  Here we consider 4 options as in adaptive design template n^(1b)
#  (see Fig. 1b in the manuscript).
#  Having set n=200 in adaptive design template n^(1b) as described above, 
#  this corresponds to the following four choices for stage 2 enrollment, 
#  where each row corresponds to a different stage 2 enrollment choice,
#  encoded as (subpopulation 1 sample size, subpopulation 2 sample size)
stage.2.sample.sizes.per.enrollment.choice = matrix(
  c(50, 50,  # Stage 2: enroll 50 from each subpopulation
    0, 0,    # Stop trial after stage 1
    150, 0,  # Stage 2: enroll 150 from subpopulation 1 and 0 from subpopulation 2
    0, 150), # Stage 2: enroll 0 from subpopulation 1 and 150 from subpopulation 2
  nrow = 4,
  ncol = 2,
  byrow = TRUE,
  dimnames = list(
    c(),
    c(
      "Subpopulation1Stage2SampleSize", # Note: these labels are optional
      "Subpopulation2Stage2SampleSize"  # and just for illustration
    )
  )
)

# We next set the minimum, clinically meaningful treatment effect size.
#  to beDelta_min = 4*qnorm(0.95)/sqrt(n) = 0.465. We explain this choice
#  below, which was made to match the distributions in Example 3.2 and 
#  the non-centrality parameter choice in Section 5.1.
Delta_min = 0.465

# The data generating distributions where the trial design will be 
#  evaluated are input next, encoded as the matrix data.generating.distributions
#  Each row of the matrix corresponds to a different data generating distribution. 
#  The user decides on how many rows to put; we use 4 rows below.
#  Each column gives a feature of the corresponding distribution, where 
#  the first two columns represent (Delta_1,Delta_2)
#  and the last 4 columns are the outcome variances 
#  sigma_10^2, sigma_11^2, sigma_20^2, sigma_21^2,
#  where sigma_sa^2 is the outcome variance for subpopulation s and arm a.
#  In our example, we set each variance equal to 1.
data.generating.distributions = matrix(
  data = c(
    # Data generating distribution with (Delta_1,Delta_2)=(0,0):
    0,0,1,1,1,1,
    # Data generating distribution with (Delta_1,Delta_2)=(Delta_min,0):
    Delta_min,0,1,1,1,1,
    # Data generating distribution with (Delta_1,Delta_2)=(0,Delta_min):
    0,Delta_min,1,1,1,1,
    # Data generating distribution with (Delta_1,Delta_2)=(Delta_min,Delta_min):
    Delta_min,Delta_min,1,1,1,1
  ),
  nrow = 4,
  ncol = 6,
  byrow = TRUE,
  dimnames = list(
    c(),
    c(
      "Delta1",  # Note: these labels are optional
      "Delta2",
      "Variance10",
      "Variance11",
      "Variance20",
      "Variance21"
    )
  )
)

# The above choice of Delta_min=0.465, sigma_sa=1 for each subpopulation s and arm a,
# and n=200 were selected in order for the non-centrality parameter (defined in 
# Section 5.1) to equal 2.33, so that our problem formulation here matches 
# Example 3.2 under the setup in Section 5 of the manuscript.
# The difference here is that we use a 
# coarse discretization so that the computation runs relatively quickly on a single
# processor using the open-source GLPK solver for illustration. We used multiple 
# processors, used Cplex (which is substantially faster than GLPK for our problems), 
# and ran our computations for longer durations to solve the problems in the paper; 
# this allowed finer discretizations. 

# Required Familywise Type I error:
total.alpha = 0.05

# Power Constraints: 
#  The power constraints are represented by a matrix power.constraints 
#  with the same number ofrows as data.generating.distributions.
#  The 3 columns correspond to the following:
#  Column 1: Power for H01; Column 2: Power for H02; Column 3: Power for H0C.
#  Each row of power.constraints represents the minimum required power to reject 
#  H01, H02, H0C, respectively, when the data generating distribution is 
#  the corresponding row of data.generating.distributions. 
#  There are four power constraints below (the first of which is vacuous)
#  For illustration these are set lower than typical (only 60%), since 
#  we use a coarse discretization that cannot achieve much higher power.
power.constraints = matrix(
  c(
    # No power requirements under first data generating distribution:
    0,0,0,
    # 60% power required for rejecting H01 under 2nd data generating distribution:
    0.6,0,0,
    # 60% power required for rejecting H02 under 3rd data generating distribution:
    0,0.6,0,
    # 60% power required for rejecting H0C under 4th data generating distribution:
    0,0,0.6
  ),
  nrow = 4,
  ncol = 3,
  byrow = TRUE,
  dimnames = list(c(), c("PowerH01", "PowerH02", "PowerH0C"))
)

# Objective Function:
#  The prior distribution Lambda used to define the objective function 
#  is a mixture of k distributions, where k is the number of rows of 
#  data.generating.distributions, and with weights
#  defined by objective.function.weights (which must sum to 1).
#  We next describe how each component distribution of the mixture is encoded.
#  Each component distribution is a bivariate normal distribution 
#  (where we allow degenerate distributions, e.g., a point mass).
#  Each component distribution is centered at (Delta_1,Delta_2) from 
#  the corresponding row of data.generating.distributions.
#  The covariance matrix of the each component distribution is defined by
#  prior.covariance.matrix, e.g., we set 
prior.covariance.matrix = diag(2) # the 2x2 identity matrix
#  When prior.covariance.matrix is input as a single 2x2 matrix, it is assumed that
#  this is the same covariance matrix for each component distribution. In this way,
#  if we also set 
objective.function.weights = 0.25 * c(1, 1, 1, 1)
#  then we obtain the Lambda corresponding to the objective function in Example 3.2.
#  The user can instead decide to set each component distribution 
#  to have a different covariance matrix,
#  which can be input by setting prior.covariance.matrix to a list of k 
#  2x2 matrices, where k is the number of rows of data.generating.distributions.
#  The special case of point masses can be encoded by setting the covariance matrix
#  (or matrices) to contain all 0's, e.g., if we set 
#  prior.covariance.matrix to the 2x2 matrix of all 0's, then Lambda 
#  in our example here would correspond to an equally weighted mixture of 
#  4 point masses as in Example 3.1. One does not have to use all 
#  of the rows of data.generating.distructions, e.g., one could set some 
#  components of objective.function.weights to 0, which effectively 
#  ignores those rows in defining the objective function prior Lambda.

#  Our software is compatible with 4 linear program solvers: Matlab, Cplex, Gurobi,
#  GLPK. Here we use GLPK for illustration.
type.of.LP.solver = "glpk"

# discretization.parameter sets how fine/coarse the rectangle partitions are
#  for the decision and rejection regions (which are exactly as in paragraph 2 in Section 5.2 
#  except that here we allow to multiply the side-lengths of squares by constants defined next), 
#  and also how fine the grid G of Type I error constraints should be. 
#  The first component is twice the side-length of squares in the box [-3,3]x[-3,3] and equal to
#  the side-lengths of squares in the remaining region of [-6,6] x [-6,6];
#  in the decision region partition;
#  the second component is the side-length of squares in the rejection region partitions;
#  the third component (after we multiply by theconstant 54) 
#  gives the number of grid points in G. 
#  There is also the option to set more complex partitions and grids; please 
#  run help(optimize_design) and look at the arguments: list.of.rectangles.dec
#   and ncp.list for explanation of available options for this.
#  For illustration below,
#  we set a coarse distribution since the corresponding problem can be generated
#  and solved in about 4 minutes on 4 core, 2.8 GHz processor on a Macbook laptop.
discretization.parameter = c(3, 3, 1)

#  Number of available cores (the more available, the faster the algorithm will work)
number.cores = min(parallel::detectCores(), 4)

# The following call to the main function optimize_design in our package
#  solves the adaptive design optimization problem using the above inputs:
optimized.policy <- optimize_design(
  subpopulation.1.proportion,
  total.alpha,
  data.generating.distributions,
  stage.1.sample.sizes,
  stage.2.sample.sizes.per.enrollment.choice,
  objective.function.weights,
  power.constraints,
  type.of.LP.solver,
  discretization.parameter,
  number.cores,
  prior.covariance.matrix = prior.covariance.matrix
)
#> [1] "Adaptive Design Optimization Completed. Optimal design is stored in the file: optimized_design1.rdata"
#> [1] "Feasible Solution was Found and Optimal Expected Sample Size is 181.23"
#> [1] "Fraction of solution components with integral value solutions"
#> [1] 0.957
#> [1] "User defined power constraints (desired power); each row corresponds to a data generating distribution; each column corresponds to H01, H02, H0C desired power, respectively."
#>            Delta1 Delta2 Variance10 Variance11 Variance20 Variance21
#> Scenario 1  0.000  0.000          1          1          1          1
#> Scenario 2  0.465  0.000          1          1          1          1
#> Scenario 3  0.000  0.465          1          1          1          1
#> Scenario 4  0.465  0.465          1          1          1          1
#>            PowerH01 PowerH02 PowerH0C
#> Scenario 1      0.0      0.0      0.0
#> Scenario 2      0.6      0.0      0.0
#> Scenario 3      0.0      0.6      0.0
#> Scenario 4      0.0      0.0      0.6
#> [1] "Probability of rejecting each null hypothesis (last 3 columns) under each data generating distribution (row)"
#>            Delta1 Delta2 Variance10 Variance11 Variance20 Variance21   H01
#> Scenario 1  0.000  0.000          1          1          1          1 0.028
#> Scenario 2  0.465  0.000          1          1          1          1 0.600
#> Scenario 3  0.000  0.465          1          1          1          1 0.050
#> Scenario 4  0.465  0.465          1          1          1          1 0.525
#>              H02   H0C
#> Scenario 1 0.028 0.023
#> Scenario 2 0.050 0.287
#> Scenario 3 0.600 0.173
#> Scenario 4 0.525 0.617
```

We next explain the output of the above function. It prints (as shown above)
a message indicating whether the problem was feasible, and if so, what the 
expected sample size is. When the problem is feasible, the power achieved under
each data generating distribution for each null hypothesis is given.
Also, the optimzed policy is returned by the software. 

The optimized.policy returned by our software consists of the pair of functions pi_1 and pi_2 as defined in Section 2.1 of the manuscript, along with state spaces and actions spaces for the sequential decision problem denoted by S1, A1, S2, A2. Each state in S1 is a rectangle with lower-left coordinates (x1,y1) encoded as optimized.policy\$lower_boundaries and upper right coordinates (x2,y2) encoded as optimized.policy\$upper_boundaries. Each action in A1 is an integer among 1,2,...,d corresponding to each of the d stage 2 enrollment choices that were input by the user. Since we allow the stage 2 state space to depend on end of stage 1 action in A1, we have that S2 is a function from an action in A1. E.g., the stage 2 state space after stage 1 action a1=3 is encoded as optimized.policy\$S2[[3]]. Each action in A2 is an integer among 1,2,...,7 encoding the following outcome of the multiple testing procedure: Reject none,Reject H01,Reject H02,Reject H0C,Reject H01 and H0C,Reject H02 and H0C,Reject all, respectively.

Below are examples illustrating how to read the output of our solver, i.e., the optimized adaptive enrichment design encoded as optimized.policy:


```r
#Lower-left and upper right coordinates of 20th rectangle in list of stage 1 states S1:
optimized.policy$S1[[20]]$lower_boundaries
#> [1] -1.5  1.5
optimized.policy$S1[[20]]$upper_boundaries
#> [1] 0 3

# Stochastic policy pi_1 evaluated at 20th rectangle in S1, i.e., the probabilities 
# of taking each action in A1 given that the first stage z-statistics are in the 
# rectangle corresponding to this state:
round(optimized.policy$pi_1(20),4)
#> [1] "Probabilities of enrollment decisions 1 through  4  respectively:"
#> [1] 0 0 0 1

# Lower-left and upper right coordinates of 4th rectangle in list of stage 2 states
# S2 following action a1=4:
a1 <- 4
optimized.policy$S2[[4]][[a1]]$lower_boundaries
#> [1] -3  0
optimized.policy$S2[[4]][[a1]]$upper_boundaries
#> [1] 0 3

# Stochastic policy pi_2 evaluated at 20th rectangle in S1, action a1=4, and 4th  rectangle in S2:
optimized.policy$pi_2(20,a1,4)
#> [1] "Inputs correspond to the following:"
#> [1] "s1 is rectangle defined as the Cartesian product [ -1.5 , 0 ] x [ 1.5 , 3 ]"
#> [1] "a1 is enrollment decision  4"
#> [1] "s2 is rectangle defined as the Cartesian product [ -3 , 0 ] x [ 0 , 3 ]"
#>        Reject none         Reject H01         Reject H02 
#>                  1                  0                  0 
#>         Reject H0C Reject H01 and H0C Reject H02 and H0C 
#>                  0                  0                  0 
#>         Reject all 
#>                  0
```
