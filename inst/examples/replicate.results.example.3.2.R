#The example below replicates the results from Example 3.2. 
#You will need to input the path to your linear program solver below
# or set LP.solver.path=NULL if this is unnecessary

LP.solver.path='path_to_your_linear_program_solver'

#Install R package if not already done so:
#remotes::install_github("mrosenblum/AdaptiveDesignOptimizerSparseLP")

# Load R package
library(AdaptiveDesignOptimizerSparseLP)
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
Delta_min = sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5; 
# The corresponding data generating distributions are encoded as follows (where we set the outcome variance to 1 for each subpopulation bby study arm combination):
data.generating.distributions=matrix(data=c(0,0,1,1,1,1,
                                           0,Delta_min,1,1,1,1,
                                           Delta_min,0,1,1,1,1,
                                           Delta_min,Delta_min,1,1,1,1),nrow=4,ncol=6,byrow=TRUE,dimnames=list(c(),c("Delta1","Delta2","Variance10","Variance11","Variance20","Variance21")));
# The resulting non-centrality parameter (see Section 5.1 of the paper) matches that used in the paper computations.
# Required Familywise Type I error:
total.alpha=0.05;
# We set the desired power for each of the power constraint to be 0.83 in the first four iterations of our algorithm, rather than the intended final output of 0.82, since we anticipate rounding in the final iteration that leads to a slight decrease in power. 
desired.power=0.83;
power.constraints=matrix(c(0,0,0,
			   0,desired.power,0,
			   desired.power,0,0,
			   0,0,desired.power),nrow=4,ncol=3,byrow=TRUE,dimnames=list(c(),c("PowerH01","PowerH02","PowerH0C")));
objective.function.weights=0.25*c(1,1,1,1);
prior.covariance.matrix=diag(2);
type.of.LP.solver="cplex";
discretization.parameter=c(1,1,10);
number.cores=30;

# Run first iteration solving sparse linear program
optimal_policy1 <- optimize_design(subpopulation.1.proportion,total.alpha,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter,number.cores,ncp.list=c(),list.of.rectangles.dec=c(),LP.iteration=1,prior.covariance.matrix,LP.solver.path=LP.solver.path)
# Construct new list of familywise Type I error constraints focusing on neighborhoods of active constraints from previous iteration 
# The Active Type I error constraints are contained in the output file optimized.design1.rdata and stored in the variable ncp.active.FWER.constraints, and in this case are the following 3 pairs: 0 0; 0.000000 1.960894; 1.960894 0.000000
# We construct a list of new Type I error constraints in neighborhoods around these active constraints, for use in the next iteration (below)
ncp.list <- list(c(0,0))
for(z in seq(1.5,3,length=30)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)))}
for(z in seq(-9,9,length=30)) {ncp.list <- c(ncp.list,list(c(z,-z)))}

# Run second iteration solving sparse linear program as in first iteration but using Type I error constraints defined above
optimal_policy2 <- optimize_design(subpopulation.1.proportion,total.alpha,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter,number.cores,ncp.list,list.of.rectangles.dec=c(),LP.iteration=2,prior.covariance.matrix,LP.solver.path=LP.solver.path)
# Refine Decision Region Rectangle Set
load("optimized.design2.rdata")
list.of.rectangles.dec <- refine_decision_rectangles(subpopulation.1.proportion=0.5,stage.1.sample.sizes=stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice=stage.2.sample.sizes.per.enrollment.choice,list.of.rectangles.dec=list.of.rectangles.dec,sln=sln,set.rectangles.with.identically.valued.neighbors.and.split.others=TRUE,LP.iteration=2)
save(list.of.rectangles.dec,file="list.of.rectangles.dec2.rdata")
# Construct new list of familywise Type I error constraints focusing on neighborhoods of active constraints from previous iteration 
# The Active Type I error constraints are contained in the output file optimized.design2.rdata and stored in the variable ncp.active.FWER.constraints
ncp.list <- list(c(0,0))
for(z in seq(1.5,3,length=30)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)))}
for(z in seq(-9,9,length=10)) {ncp.list <- c(ncp.list,list(c(z,-z)))}

# Run third iteration solving sparse linear program as in first iteration but using Type I error constraints defined above
optimal_policy3 <- optimize_design(subpopulation.1.proportion,total.alpha,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter,number.cores,ncp.list,list.of.rectangles.dec,LP.iteration=3,prior.covariance.matrix,LP.solver.path=LP.solver.path)
# Refine Decision Region Rectangle Set
load("optimized.design3.rdata")
list.of.rectangles.dec <- refine_decision_rectangles(subpopulation.1.proportion=0.5,stage.1.sample.sizes=stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice=stage.2.sample.sizes.per.enrollment.choice,list.of.rectangles.dec=list.of.rectangles.dec,sln=sln,set.rectangles.with.identically.valued.neighbors.and.split.others=TRUE,LP.iteration=3)
save(list.of.rectangles.dec,file="list.of.rectangles.dec3.rdata")
# Construct new list of familywise Type I error constraints focusing on neighborhoods of active constraints from previous iteration 
# The Active Type I error constraints are contained in the output file optimized.design3.rdata and stored in the variable ncp.active.FWER.constraints
ncp.list <- list(c(0,0))
for(z in seq(1.8,2.5,length=20)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)))}

# Run fourth iteration solving sparse linear program as in first iteration but using Type I error constraints defined above
optimal_policy4 <-optimize_design(subpopulation.1.proportion,total.alpha,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter,number.cores,ncp.list,list.of.rectangles.dec,LP.iteration=4,prior.covariance.matrix,LP.solver.path=LP.solver.path)
# Set Decision Region Rectangles 
load("optimized.design4.rdata")
list.of.rectangles.dec <- refine_decision_rectangles(subpopulation.1.proportion=0.5,stage.1.sample.sizes=stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice=stage.2.sample.sizes.per.enrollment.choice,list.of.rectangles.dec=list.of.rectangles.dec,sln=sln,set.rectangles.with.identically.valued.neighbors.and.split.others=FALSE,round.each.decision.rectangle.to.integer=TRUE,LP.iteration=4)
save(list.of.rectangles.dec,file="list.of.rectangles.dec4.rdata")
# Construct new list of familywise Type I error constraints focusing on neighborhoods of active constraints from previous iteration 
# The Active Type I error constraints are contained in the output file optimized.design4.rdata and stored in the variable ncp.active.FWER.constraints
ncp.list <- list(c(0,0))
for(z in seq(2,2.5,length=500)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)))}
for(z in seq(-9,9,length=950)) {ncp.list <- c(ncp.list,list(c(0,z)),list(c(z,0)),list(c(z,-z)))}
for(z in c(seq(-3,-1,length=500),seq(1,3,length=500))) {ncp.list <- c(ncp.list,list(c(z,-z)))}

# Run fifth iteration solving sparse linear program as in first iteration but using Type I error constraints defined above
optimal_policy5 <- optimize_multiple_testing_procedure(subpopulation.1.proportion,total.alpha=0.049,data.generating.distributions,stage.1.sample.sizes,stage.2.sample.sizes.per.enrollment.choice,objective.function.weights,power.constraints,type.of.LP.solver,discretization.parameter=c(1,0.25,10),number.cores,ncp.list,list.of.rectangles.dec,LP.iteration=5,prior.covariance.matrix,round.each.multiple.testing.procedure.rectangle.to.integer=TRUE,plots.to.round.simply = c(1,2),rounding.threshold.H01 = 1-1e-10,rounding.threshold.H02 = 1-1e-10,rounding.threshold.H0C = 1-1e-10,power.constraint.tolerance = 0.01,LP.solver.path=LP.solver.path)