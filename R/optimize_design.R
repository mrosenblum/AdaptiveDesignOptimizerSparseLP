#' Adaptive Enrichment Design Optimization Using Sparse Linear Programming
#' Authors: Michael Rosenblum, Ethan Fang, Han Liu
#'
#' @param subpopulation.1.proportion Proportion of overall population in subpopulation 1. Must be between 0 and 1.
#' @param total.alpha Familywise Type I error rate (1-sided)
#' @param data.generating.distributions Matrix encoding data generating distributions (defined in terms of treatment effect pairs and outcome variances) used to define power constraints and  objective function; each row defines the pair (Delta_1,Delta_2) of subpopulation 1 and 2 average treatment effects, followed by outcome variances for the four combinations of subpouplation (1 and 2) by study arm (0 and 1).
#' @param stage.1.sample.sizes Vector with 2 entries representing stage 1 sample sizes for subpopulations 1 and 2, respectively
#' @param stage.2.sample.sizes.per.enrollment.choice Matrix with number.choices.end.of.stage.1 rows and 2 columns, where the (i,j) entry represents the stage 2 sample size under enrollment choice i for subpopulation j.
#' @param objective.function.weights Vector with length equal to number of rows of population.parameters, representing weights used to define the objective function
#' @param power.constraints Matrix with same number of rows as population.parameters (each representing a data generating distribution) and three columns corresponding to the required power to reject (at least) H_01, H_02, H_0C, respectively.
#' @param type.of.LP.solver "matlab", "cplex", "GLPK", or "Gurobi" The linear program solve that you want to use; assumes that you have installed this already and that path is set
#' @param discretization.parameter vector with 3 elements representing initial discretization of decision region, rejection regions, and grid representing Type I error constraints
#' @param number.cores the number of cores available for parallelization using the parallel R package
#' @param ncp.list list of pairs of real numbers representing the non-centrality parameters to be used in the Type I error constraints; if list is empty, then default list is used.
#' @param list.of.rectangles.dec list of rectangles representing decision region partition, encoded as a list with each element of the list having fields $lower_boundaries (pair of real numbers representing coordinates of lower left corner of rectangle), $upper_boundaries (pair of real numbers representing upper right corner of rectangle), $allowed_decisions (subset of stage.2.sample.sizes.per.enrollment.choice representing which decisions allowed if first stage z-statistics are in corresponding rectangle; default is entire list stage.2.sample.sizes.per.enrollment.choice), $preset_decision (indicator of whether the decision probabilities are hard-coded by the user; default is 0), $d_probs (empty unless $preset_decision==1, in which case it is a vector representing the probabilities of each decision); if list.or.rectangles.dec is empty, then a default partition is used based on discretization.parameter.
#' LP.iteration positive integer used in file name to store output; can be used to avoid overwriting previous computations
#' round.each.decision.rectangle.to.integer TRUE/FALSE indicator of whether decision probabilities encoded in list.of.rectangles.dec should be rounded to integer values
#' set.rectangles.with.identically.valued.neighbors.and.split.others  TRUE/FALSE indicator of whether decision probabilities encoded in list.of.rectangles.dec should be modified for use in next iteration
#' @return 4 element list containing optimized designs from four classes (with increasing complexity):
#' @section Output
#' The software computes and optimized design saved as "optimized_design.rdata" and the corresponding expected sample size is
#' saved as "optimized_design_expected_sample_size.rdata".
#' @examples
#' #For demonstration purposes, the examples below use a coarse discretization.
#' optimize_design(discretization.parameter=c(3,3,1),number.cores=1)
#' @export
optimize_design <- function(subpopulation.1.proportion=0.5,
		total.alpha=0.05-(1e-4),
		data.generating.distributions=matrix(data=c(0,0,1,1,1,1,
					       0,sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5,1,1,1,1,
					       sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5,0,1,1,1,1,
					       sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5,sqrt(1/2)*(qnorm(0.95+1e-4)+qnorm(0.95))/5,1,1,1,1),nrow=4,ncol=6,byrow=TRUE,dimnames=list(c(),c("Delta1","Delta2","Variance10","Variance11","Variance20","Variance21"))),
		stage.1.sample.sizes=c(50,50),
		stage.2.sample.sizes.per.enrollment.choice=matrix(c(50,50,
								    0,0,
								    150,0,
		    					            0,150),nrow=4,ncol=2,byrow=TRUE,dimnames=list(c(),c("Subpopulation1Stage2SampleSize","Subpopulation2Stage2SampleSize"))),
	  objective.function.weights=0.25*c(1,1,1,1),
		power.constraints=matrix(c(0,0,0,
					   0,0.82,0,
					   0.82,0,0,
					   0,0,0.82),nrow=4,ncol=3,byrow=TRUE,dimnames=list(c(),c("PowerH01","PowerH02","PowerH0C"))),
	        type.of.LP.solver="matlab",
		discretization.parameter=c(1,1,10),
		number.cores=30,
		ncp.list=c(),
		list.of.rectangles.dec=c(),
		LP.iteration=1,
		round.each.decision.rectangle.to.integer=FALSE,
		set.rectangles.with.identically.valued.neighbors.and.split.others=FALSE
		){

max_error_prob <- 0 # track approximation errors in problem construction; initialize to 0 here
covariance_Z_1_Z_2 <-  0 # covariance_due_to overlap of subpopulations (generally we assume no overlap)
p1 <- subpopulation.1.proportion;
p2 <- 1-p1 # proportion of population in subpopulation 2 (assumes non-overlapping subpopulations that partition entire population)
rho1 <- sqrt(p1)
rho2 <- sqrt(p2)

actions <- c(1,2,3,4,5,6,7)
number_actions <- length(actions) # correspond to rejecting nothing, H01, H02, H0C, H01 and H0C, H02 and H0C, all

## ncp (shorthand for non-centrality parameter) is a 2 component vector equal to c(sqrt(2*p_1*sum(stage.1.sample.sizes))Delta_1/(2sigma_1) , sqrt(2*p_2*sum(stage.1.sample.sizes)Delta_2/(2sigma_2))) for sigma^2_s=(sigma^2_(s0)+sigma^2_(s1))/2.

## It represents c(EZ_1,EZ_2) for the fixed design that enrolls 2*p_s*sum(stage.1.sample.sizes) from each subpopulation s.

map_from_P_to_type_I_error_indicator_over_set_of_actions <- function(ncp)
{
	return(c(0,ncp[1]<=0,ncp[2]<=0,rho1*ncp[1]+rho2*ncp[2]<=0, ncp[1]<=0 || rho1*ncp[1]+rho2*ncp[2]<=0, ncp[2] <= 0 || rho1*ncp[1]+rho2*ncp[2]<=0, ncp[1]<=0 || ncp[2]<=0 || rho1*ncp[1]+rho2*ncp[2]<=0))
}

indicator_contribute_to_H01_power <- c(0,1,0,0,1,0,1)
indicator_contribute_to_H02_power <- c(0,0,1,0,0,1,1)
indicator_contribute_to_H0C_power <- c(0,0,0,1,1,1,1)

# Sample Sizes per stage and decision:
# where sample sizes are in units of sum(stage.1.sample.sizes)/2
n_stage1_subpopulation1 <- stage.1.sample.sizes[1];
n_stage1_subpopulation2 <- stage.1.sample.sizes[2];
n_stage2_subpopulation1_decision <- as.vector(stage.2.sample.sizes.per.enrollment.choice[,1]) #Sample size in stage 2 under each decision rule for subpopulation 1
n_stage2_subpopulation2_decision <- as.vector(stage.2.sample.sizes.per.enrollment.choice[,2]) #Sample size in stage 2 under each decision rule for subpopulation 1
number_decisions <- ifelse(length(n_stage2_subpopulation1_decision)==length(n_stage2_subpopulation2_decision),length(n_stage2_subpopulation1_decision),exit())
decisions <- (1:number_decisions)

## Set loss function to be sample size; can modify if desired--general format is matrix with number_decision rows and number_actions columns, and entry is loss function value at corresponding (decision,action pair). Can also generalize to make it depend on the ncp value as well but if so need to modify generalized_generate... objective function construction
loss_function_value <- array(n_stage1_subpopulation1+n_stage1_subpopulation2+n_stage2_subpopulation1_decision+n_stage2_subpopulation2_decision,c(number_decisions,number_actions))

# Convert data.generating.distributions to list of non-centrality parameter vectors
prior_mean_support <- c()
for(count in 1:dim(data.generating.distributions)[1]){
  ncp <- c(data.generating.distributions[count,1]*sqrt(p1*sum(stage.1.sample.sizes)/(data.generating.distributions[count,3]+data.generating.distributions[count,4])),
           data.generating.distributions[count,2]*sqrt(p2*sum(stage.1.sample.sizes)/(data.generating.distributions[count,5]+data.generating.distributions[count,6])))
  names(ncp) <- c("NonCentralityParameter1","NonCentralityParameter2")
  prior_mean_support <- c(prior_mean_support,list(ncp))
}

# Prior information used in constructing objective function
prior_weights <- objective.function.weights
prior_covariance_matrix <- diag(2)*0
# List of non-centrality parameters to use for power constraints:
# WARNING: If following line is modified, need to modify the line that follows it as well, and also the line:
# switch(power_constraint_counter,power_constraint_vector_subpopulation1<-power_constraint_vector,power_constraint_vector_subpopulation2<-power_constraint_vector,power_constraint_vector_combined_population<-power_constraint_vector)
# in the generalized_generate... file in order to modify the format of power constraint output
power_constraint_list <- prior_mean_support
## The type of power to put constraint on for each scenario in the power_constraint list:
#power_constraint_null_hyp_contribution_list <- c(list(indicator_contribute_to_subpopulation1_power),list(indicator_contribute_to_subpopulation2_power),list(indicator_contribute_to_combined_population_power))

# Decision Type: indicates one of 4 types of stage 2 enrollment choices:
# 1: both subpopulations enrolled
# 2: neither subpopulation enrolled
# 3: just subpopulation 1 enrolled
# 4: just subpopulation 2 enrolled
# Each is handled differently, e.g., decisions 3 and 4 use a 3x3 covariance matrix while decision 1 uses 4x4.

d_type <- ifelse(n_stage2_subpopulation1_decision>0 & n_stage2_subpopulation2_decision>0,1,
          ifelse(n_stage2_subpopulation1_decision==0 & n_stage2_subpopulation2_decision==0,2,
          ifelse(n_stage2_subpopulation1_decision>0 & n_stage2_subpopulation2_decision==0,3,
          ifelse(n_stage2_subpopulation1_decision==0 & n_stage2_subpopulation2_decision>0,4,5))))

# Create covariance matrices under each decision rule
covariance_matrix <- list()
for(d in decisions){
  if(d_type[d]==1){
  # Corresponds to decision to enroll from both subpopulations in stage 2
  gamma1 <- sqrt(n_stage1_subpopulation1/(n_stage1_subpopulation1+n_stage2_subpopulation1_decision[d]))
  gamma2 <- sqrt(n_stage1_subpopulation2/(n_stage1_subpopulation2+n_stage2_subpopulation2_decision[d]))
  ## Statistics: Z^(1)_1,Z^(1)_2,Z^(C)_1,Z^(C)_2,
  covariance_matrix <- c(covariance_matrix,list(
                         array(c(1,covariance_Z_1_Z_2,gamma1,covariance_Z_1_Z_2*gamma2,
                      					    covariance_Z_1_Z_2,1,covariance_Z_1_Z_2*gamma1,gamma2,
                     					    gamma1,covariance_Z_1_Z_2*gamma1,1,covariance_Z_1_Z_2,
                     					    covariance_Z_1_Z_2*gamma2,gamma2,covariance_Z_1_Z_2,1),c(4,4))))
  }else if(d_type[d]==2){
  # Corresponds to decision to stop early
    covariance_matrix <- c(covariance_matrix,list(
                           array(c(1,covariance_Z_1_Z_2,
                      					    covariance_Z_1_Z_2,1),c(2,2))))
  }else if(d_type[d]==3){
  # Corresponds to decision to enrich subpopulation 1 only
  gamma1 <- sqrt(n_stage1_subpopulation1/(n_stage1_subpopulation1+n_stage2_subpopulation1_decision[d]))
  ## Statistics: Z^(1)_1,Z^(1)_2,Z^(C)_1,Z^(C)_2,
  covariance_matrix <- c(covariance_matrix,list(
                         array(c(1,covariance_Z_1_Z_2,gamma1,
                      					 covariance_Z_1_Z_2,1,covariance_Z_1_Z_2*gamma1,
                     					  gamma1,covariance_Z_1_Z_2*gamma1,1),c(3,3))))
  }else if(d_type[d]==4){
  # Corresponds to decision to enrich subpopulation 2 only
  gamma2 <- sqrt(n_stage1_subpopulation2/(n_stage1_subpopulation2+n_stage2_subpopulation2_decision[d]))
  ## Statistics: Z^(1)_1,Z^(1)_2,Z^(C)_1,Z^(C)_2,
  covariance_matrix <- c(covariance_matrix,list(
                         array(c(1,covariance_Z_1_Z_2,covariance_Z_1_Z_2*gamma2,
                      					    covariance_Z_1_Z_2,1,gamma2,
                     					    covariance_Z_1_Z_2*gamma2,gamma2,1),c(3,3))))
  }
}

mean_vector <- function(ncp,d){
    if(d_type[d]==1){
    return(c(
       ncp[1]*sqrt(n_stage1_subpopulation1/(2*p1*sum(stage.1.sample.sizes))),ncp[2]*sqrt(n_stage1_subpopulation2/(2*p2*sum(stage.1.sample.sizes))),
       ncp[1]*sqrt((n_stage1_subpopulation1+n_stage2_subpopulation1_decision[d])/(2*p1*sum(stage.1.sample.sizes))),ncp[2]*sqrt((n_stage1_subpopulation2+n_stage2_subpopulation2_decision[d])/(2*p2*sum(stage.1.sample.sizes)))))} else if(d_type[d]==2){
    return(c(
       ncp[1]*sqrt(n_stage1_subpopulation1/(2*p1*sum(stage.1.sample.sizes))),ncp[2]*sqrt(n_stage1_subpopulation2/(2*p1*sum(stage.1.sample.sizes)))))} else if(d_type[d]==3){
    return(c(
       ncp[1]*sqrt(n_stage1_subpopulation1/(2*p1*sum(stage.1.sample.sizes))),ncp[2]*sqrt(n_stage1_subpopulation2/(2*p2*sum(stage.1.sample.sizes))),
       ncp[1]*sqrt((n_stage1_subpopulation1+n_stage2_subpopulation1_decision[d])/(2*p1*sum(stage.1.sample.sizes)))))} else if(d_type[d]==4){
    return(c(
       ncp[1]*sqrt(n_stage1_subpopulation1/(2*p1*sum(stage.1.sample.sizes))),ncp[2]*sqrt(n_stage1_subpopulation2/(2*p2*sum(stage.1.sample.sizes))),
       ncp[2]*sqrt((n_stage1_subpopulation2+n_stage2_subpopulation2_decision[d])/(2*p2*sum(stage.1.sample.sizes)))))
       }
}

# L <- function(ncp){return(
#    -(as.numeric(ncp[1]>delta_1_min)*indicator_contribute_to_subpopulation1_power
#    +as.numeric(ncp[2]>delta_2_min)*indicator_contribute_to_subpopulation2_power
#    +as.numeric(rho1*ncp[1]+rho2*ncp[2]>rho1*delta_1_min+rho2*delta_2_min)*indicator_contribute_to_combined_population_power))}

## Handles case of prior distribution used in objective function
modified_joint_distribution <- function(prior_component_index,decision){
	modified_mean_vector <- mean_vector(ncp=prior_mean_support[[prior_component_index]],d=decision)
	modified_covariance_matrix <-  (array(c(mean_vector(ncp=c(1,0),d=decision),mean_vector(ncp=c(0,1),d=decision)),c(length(mean_vector(ncp=c(1,0),d=decision)),2))  %*% prior_covariance_matrix %*% t(array(c(mean_vector(ncp=c(1,0),d=decision),mean_vector(ncp=c(0,1),d=decision)),c(length(mean_vector(ncp=c(1,0),d=decision)),2)))) +  covariance_matrix[[decision]]
	return(list(modified_mean_vector,modified_covariance_matrix))
}

number_reference_rectangles <- number_decisions

## Below set the FWER constraints and the discretization of the decision and rejection regions

#set widths of large square outside of which multiple testing procedure rejects no null hypotheses (corresponds to w in text)
w1 <- 6
w2 <- 6
# dimension of small squares used in discretization
tau <- discretization.parameter[1]
tau_mtp <- discretization.parameter[2]
# Settings for integration region and precision in computing lower bound to original problem using dual solution to discretized problem
max_eval_iters <- 100000
w1_unconst <- 5
w2_unconst <- 5
constraints_per_A1_file <- 1 # Set number of familywise Type I error constraints to encode per file written

if(is.null(ncp.list)){
  # list of pairs of non-centrality parameters in G_{tau,w}
  ncp.list <- list()
  for(z in seq(-9,9,length=18*discretization.parameter[3])) {ncp.list <- c(ncp.list,list(c(0,z)))}
  for(z in seq(-9,9,length=18*discretization.parameter[3])) {ncp.list <- c(ncp.list,list(c(z,0)))}
  for(z in seq(-9,9,length=18*discretization.parameter[3])) {ncp.list <- c(ncp.list,list(c(z,-(rho1/rho2)*z)))}
  ncp.list <- c(ncp.list,list(c(0,0)))
  ncp.list <- unique(ncp.list)
}

  # construct list of rectangles in set R
  ## List of rectangles defining decision boundaries
  ## Each rectangle is encoded by c(lower left (x,y) coordinates, upper right (x,y) coordinates)

if(is.null(list.of.rectangles.dec)){
  list.of.rectangles.dec <- list()

  number_preset_decision_rectangles <- 0

  for(x in c(seq(-w1,-3-tau,by=tau),seq(3,w1-tau,by=tau)))
  {
	for(y in seq(-w2,w2-tau,by=tau))
	{
		list.of.rectangles.dec<- c(list.of.rectangles.dec,list(list(lower_boundaries=c(x,y),upper_boundaries=c(x+tau,y+tau),allowed_decisions=decisions,preset_decision=0)))
	}
  }

  for(x in seq(-3,3-tau,by=tau))
    {
    for(y in c(seq(-w1,-3-tau,by=tau),seq(3,w1-tau,by=tau)))
	{
		list.of.rectangles.dec<- c(list.of.rectangles.dec,list(list(lower_boundaries=c(x,y),upper_boundaries=c(x+tau,y+tau),allowed_decisions=decisions,preset_decision=0)))
	}
  }

  for(x in seq(-3,3-tau/2,by=tau/2))
  {
	for(y in seq(-3,3-tau/2,by=tau/2))
	{
		list.of.rectangles.dec<- c(list.of.rectangles.dec,list(list(lower_boundaries=c(x,y),upper_boundaries=c(x+tau/2,y+tau/2),allowed_decisions=decisions,preset_decision=0)))
	}
  }

  for(d in decisions){list.of.rectangles.dec <- c(list.of.rectangles.dec,list(list(lower_boundaries=c(0,0),upper_boundaries=c(0,0),allowed_decisions=d,preset_decision=0)))}# these are reference rectangles

   ## Set upper neighbors
   counter_for_r <- 1
   for(counter_for_r in 1:(length(list.of.rectangles.dec)-number_reference_rectangles)){
	r <- list.of.rectangles.dec[[counter_for_r]]
	count_value <- 1
	while(count_value <= length(list.of.rectangles.dec) && (!(r$upper_boundaries[2]== list.of.rectangles.dec[[count_value]]$lower_boundaries[2] && r$lower_boundaries[1]== list.of.rectangles.dec[[count_value]]$lower_boundaries[1] && r$upper_boundaries[1]== list.of.rectangles.dec[[count_value]]$upper_boundaries[1]))){count_value <- count_value +1}
	if(count_value <= length(list.of.rectangles.dec)){list.of.rectangles.dec[[counter_for_r]]$upper_neighbor <- count_value}
   }
   save(list.of.rectangles.dec,file=paste("list.of.rectangles.dec",LP.iteration,".rdata",sep=""))
   save(ncp.list,file=paste("ncp.list",LP.iteration,".rdata",sep=""))
}

## Set multiple testing procedure rectangle partition

stage_2_rectangle_offset_value <- 0
list.of.rectangles.mtp1 <- list()
for(x in seq(-w1,w1,by=tau_mtp))
{
	#Skip lower left quadrant
	for(y in seq(ifelse(x<0,0,-w2),w2,by=tau_mtp))
	{
		list.of.rectangles.mtp1<- c(list.of.rectangles.mtp1,list(list(lower_boundaries=c(x,y),upper_boundaries=c(x+tau_mtp,y+tau_mtp),allowed_actions=actions,stage_2_rectangle_offset=stage_2_rectangle_offset_value)))
		stage_2_rectangle_offset_value <- stage_2_rectangle_offset_value + length(actions)
	}
}
## Only one large rectangle for lower left quadrant
list.of.rectangles.mtp1<- c(list.of.rectangles.mtp1,list(list(lower_boundaries=c(-Inf,-Inf),upper_boundaries=c(0,0),allowed_actions=actions,stage_2_rectangle_offset=stage_2_rectangle_offset_value)))
stage_2_rectangle_offset_value <- stage_2_rectangle_offset_value + length(actions)

number_rectangles_by_actions_for_decision_type1  <- stage_2_rectangle_offset_value

## Set right neighbors
counter_for_rprime <- 1
for(rprime in list.of.rectangles.mtp1)
{
	count_value <- 1
	while(count_value <= length(list.of.rectangles.mtp1) && (!(rprime$upper_boundaries[1]== list.of.rectangles.mtp1[[count_value]]$lower_boundaries[1] && rprime$lower_boundaries[2]== list.of.rectangles.mtp1[[count_value]]$lower_boundaries[2] && rprime$upper_boundaries[2]== list.of.rectangles.mtp1[[count_value]]$upper_boundaries[2]))){count_value <- count_value +1}
	if(count_value <= length(list.of.rectangles.mtp1)){list.of.rectangles.mtp1[[counter_for_rprime]]$right_neighbor <- count_value}
	counter_for_rprime <- counter_for_rprime +1
}
## Set upper neighbors
counter_for_rprime <- 1
for(rprime in list.of.rectangles.mtp1)
{
	count_value <- 1
	while(count_value <= length(list.of.rectangles.mtp1) && (!(rprime$upper_boundaries[2]== list.of.rectangles.mtp1[[count_value]]$lower_boundaries[2] && rprime$lower_boundaries[1]== list.of.rectangles.mtp1[[count_value]]$lower_boundaries[1] && rprime$upper_boundaries[1]== list.of.rectangles.mtp1[[count_value]]$upper_boundaries[1]))){count_value <- count_value +1}
	if(count_value <= length(list.of.rectangles.mtp1)){list.of.rectangles.mtp1[[counter_for_rprime]]$upper_neighbor <- count_value}
	counter_for_rprime <- counter_for_rprime +1
}

##Set list.of.rectangles.mtp2
stage_2_rectangle_offset_value <- 0
list.of.rectangles.mtp2 <- list()
for(x in seq(-w1,w1,by=tau)){
    for(y in seq(-w2,w2,by=tau)){
        list.of.rectangles.mtp2<- c(list.of.rectangles.mtp2,list(list(lower_boundaries=c(x,y),upper_boundaries=c(x+tau,y+tau),allowed_actions=actions,stage_2_rectangle_offset=stage_2_rectangle_offset_value)))
        stage_2_rectangle_offset_value <- stage_2_rectangle_offset_value + length(actions)
    }
}
## IF WANT TO TEMPORARILY SET MTP2 TO SINGLE LARGE RECTANGLE, replace above loops by:
#list.of.rectangles.mtp2<- c(list.of.rectangles.mtp2,list(list(lower_boundaries=c(-Inf,-Inf),upper_boundaries=c(Inf,Inf),allowed_actions=actions,stage_2_rectangle_offset=stage_2_rectangle_offset_value)))
#stage_2_rectangle_offset_value <- stage_2_rectangle_offset_value + length(actions)

number_rectangles_by_actions_for_decision_type2  <- stage_2_rectangle_offset_value

##Set stage 2 rectangle sets for each d_type:
#Set stage 2 rectangles sets to be the same for d_type 1, 3, and 4. (Allow different for d_type 2)
list.of.rectangles.mtp3 <- list.of.rectangles.mtp4 <- list.of.rectangles.mtp1;
number_rectangles_by_actions_for_decision_type3 <- number_rectangles_by_actions_for_decision_type4 <- number_rectangles_by_actions_for_decision_type1
list.of.rectangles.mtp <- list()
number_rectangles_by_actions_for_decision <- c()
for(d in decisions){
  list.of.rectangles.mtp <- c(list.of.rectangles.mtp,list(switch(d_type[d],list.of.rectangles.mtp1,list.of.rectangles.mtp2,list.of.rectangles.mtp3,list.of.rectangles.mtp4)));
  number_rectangles_by_actions_for_decision <- c(number_rectangles_by_actions_for_decision,switch(d_type[d],number_rectangles_by_actions_for_decision_type1,number_rectangles_by_actions_for_decision_type2,number_rectangles_by_actions_for_decision_type3,number_rectangles_by_actions_for_decision_type4))
}
#Can replace above loop if desired to have tailored rectangle sets for each possible enrollment choice.
#save(list.of.rectangles.dec,file="list_of_rectangles.rdata")

# compute number of variables in linear program
number_of_variables <- 0
current_rectangle <-1
for(r in list.of.rectangles.dec){
    if(r$preset_decision == 0){# we do not encode variables corresponding to decision region rectangles that are preset, i.e, where preset_decision==1
        list.of.rectangles.dec[[current_rectangle]]$stage_1_rectangle_offset <- number_of_variables
        for(d in decisions){
            if(sum(d==r$allowed_decisions)>0){
                list.of.rectangles.dec[[current_rectangle]]$stage_1_rectangle_and_decision_offset[[d]] <- number_of_variables
                number_of_variables <- number_of_variables+number_rectangles_by_actions_for_decision[d]
            } else {list.of.rectangles.dec[[current_rectangle]]$stage_1_rectangle_and_decision_offset[[d]] <- NA}
        }
    }
    current_rectangle <- current_rectangle+1
}

variable_location <- function(r,d,rprime,action){return(r$stage_1_rectangle_and_decision_offset[d]+rprime$stage_2_rectangle_offset+which(rprime$allowed_actions==action))}

#compute number_equality_constraints_part2
number_equality_constraints_part2 <- 0
for(r in list.of.rectangles.dec){
    if(r$preset_decision==0){
        for(d in r$allowed_decisions){
            number_equality_constraints_part2<-number_equality_constraints_part2+(length(list.of.rectangles.mtp[[d]])-1)
	}
    }
}

print("number of variables")
print(number_of_variables)
print("number of familywise Type I error constraints")
print(length(ncp.list))
number_equality_constraints_part1 <- length(list.of.rectangles.dec)-number_preset_decision_rectangles
write(number_equality_constraints_part1,f=paste("number_equality_constraints_of_first_type.txt"))
print("number of equality constraints of first type")
print(number_equality_constraints_part1)
print("number of equality constraints of second type")
print(number_equality_constraints_part2)
write(number_equality_constraints_part2,f=paste("number_equality_constraints_of_second_type.txt"))
write(length(ncp.list),f=paste("number_A1_constraints.txt"))
write(ceiling(length(ncp.list)/constraints_per_A1_file),f=paste("number_A1_files.txt"))
power.constraints <- as.vector(power.constraints)
save(power.constraints,file="power_constraints.rdata")
save(list.of.rectangles.mtp,file=paste("list.of.rectangles.mtp",LP.iteration,".rdata",sep=""))

## Includes constraints the restrict multiple testing procedure not depend on stage 1 statistics (given cumulative statistics Z^C and decision d)
## Implements unequal rectangle sizes, with setting of rectangles to specific values
## Generate Discretized LP based on settings in file

generate_LP <- function(task_id,...){
    max_task_id_for_computing_FWER_constraints <- ceiling(length(ncp.list)/constraints_per_A1_file)
    print("max_task_id_for_computing_FWER_constraints")
    print(max_task_id_for_computing_FWER_constraints)

    if(task_id <= max_task_id_for_computing_FWER_constraints){
        ## Construct Familywise Type I error constraints
        constraint_list <- c()
        counter <- 1
        max_error_prob <- 0

for(ncp in ncp.list[((task_id-1)*constraints_per_A1_file+1):min(length(ncp.list),((task_id)*constraints_per_A1_file))])
{
                                        # construct vectors for storing reallocated probabilities for preset decision rectangles
    for(d in decisions){
        for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
            list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- rep(0,length(list.of.rectangles.mtp[[d]][[rprime_counter]]$allowed_actions))
        }
    }

    rvec <- c()
    ## Construct constraint that P_{ncp}(Type I error) \leq \alpha
    true_nulls_violated_by_action_set <- map_from_P_to_type_I_error_indicator_over_set_of_actions(ncp)
    for(r in list.of.rectangles.dec[1:(length(list.of.rectangles.dec)-number_reference_rectangles)]){
        for(d in r$allowed_decisions){
            if(d_type[d]==1){
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## Get P_{ncp_vec}(Z \in rectangle r)
                    prob_r <- mvtnorm::pmvnorm(mean=mean_vector(ncp,d),sigma=covariance_matrix[[d]],lower=c(r$lower_boundaries,rprime$lower_boundaries),upper=c(r$upper_boundaries,rprime$upper_boundaries),algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                    if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    if(r$preset_decision==0){rvec <- c(rvec,prob_r*true_nulls_violated_by_action_set[rprime$allowed_actions])}else{
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles + r$preset_decision_value[d]*prob_r*true_nulls_violated_by_action_set[rprime$allowed_actions]
                    }
                }
            } else if(d_type[d]==2) {
                ## Get P_{ncp_vec}(Z \in rectangle r)
                 for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## compute overlap in rectangles in terms of z_2:
                    intersection_rectangle_lower_boundaries <- pmax(r$lower_boundaries,rprime$lower_boundaries)
                    intersection_rectangle_upper_boundaries <- pmin(r$upper_boundaries,rprime$upper_boundaries)
                    if(sum(intersection_rectangle_lower_boundaries < intersection_rectangle_upper_boundaries)==2){ #case where intersection rectangle non-empty
                        prob_r <- mvtnorm::pmvnorm(mean=mean_vector(ncp,d),sigma=covariance_matrix[[d]],lower=intersection_rectangle_lower_boundaries,upper=intersection_rectangle_upper_boundaries,algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                        if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    } else {# case where rectangle empty
                        prob_r <- 0
                    }
                    if(r$preset_decision==0){rvec <- c(rvec,prob_r*true_nulls_violated_by_action_set[rprime$allowed_actions])}else{
                            list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles + r$preset_decision_value[d]*prob_r*true_nulls_violated_by_action_set[rprime$allowed_actions]
                        }
                }
            } else if(d_type[d]==3) {
                r_lower_boundary_z_2 <- r$lower_boundaries[2]
                r_upper_boundary_z_2 <- r$upper_boundaries[2]
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## compute overlap in rectangles in terms of z_2:
                    effective_z_2_lower_boundary <- max(c(r_lower_boundary_z_2,rprime$lower_boundaries[2]))
                    effective_z_2_upper_boundary <- min(c(r_upper_boundary_z_2,rprime$upper_boundaries[2]))
                    if(effective_z_2_lower_boundary < effective_z_2_upper_boundary){ # case where rectangle non-empty
                        ## Get P_{ncp_vec}(Z \in rectangle r)
                        prob_r <- mvtnorm::pmvnorm(mean=mean_vector(ncp,d),sigma=covariance_matrix[[d]],lower=c(r$lower_boundaries[1],effective_z_2_lower_boundary,rprime$lower_boundaries[1]),upper=c(r$upper_boundaries[1],effective_z_2_upper_boundary,rprime$upper_boundaries[1]),algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                        if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    } else {# case where rectangle empty
                        prob_r <- 0
                    }
                    if(r$preset_decision==0){rvec <- c(rvec,prob_r*true_nulls_violated_by_action_set[rprime$allowed_actions])}else{
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles + r$preset_decision_value[d]*prob_r*true_nulls_violated_by_action_set[rprime$allowed_actions]
                    }
                }
            } else if(d_type[d]==4) { #only subpopulation 2 enrolled further
                r_lower_boundary_z_1 <- r$lower_boundaries[1]
                r_upper_boundary_z_1 <- r$upper_boundaries[1]
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## compute overlap in rectangles in terms of z_2:
                    effective_z_1_lower_boundary <- max(c(r_lower_boundary_z_1,rprime$lower_boundaries[1]))
                    effective_z_1_upper_boundary <- min(c(r_upper_boundary_z_1,rprime$upper_boundaries[1]))
                    if(effective_z_1_lower_boundary < effective_z_1_upper_boundary){ # case where rectangle non-empty
                        ## Get P_{ncp_vec}(Z \in rectangle r)
                        prob_r <- mvtnorm::pmvnorm(mean=mean_vector(ncp,d),sigma=covariance_matrix[[d]],lower=c(effective_z_1_lower_boundary,r$lower_boundaries[2],rprime$lower_boundaries[2]),upper=c(effective_z_1_upper_boundary,r$upper_boundaries[2],rprime$upper_boundaries[2]),algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                        if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    } else {# case where rectangle empty
                        prob_r <- 0
                    }
                    if(r$preset_decision==0){rvec <- c(rvec,prob_r*true_nulls_violated_by_action_set[rprime$allowed_actions])}else{
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles + r$preset_decision_value[d]*prob_r*true_nulls_violated_by_action_set[rprime$allowed_actions]
                    }
                }
            }
        }
    }
    #Handle reference rectangle contributions
    for(r in list.of.rectangles.dec[(length(list.of.rectangles.dec)+1-number_reference_rectangles):length(list.of.rectangles.dec)]){
        for(d in r$allowed_decisions){
            for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                rvec <- c(rvec,list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles)
            }
        }
    }

    constraint_list <- rbind(constraint_list,rvec)
    print(counter)
                                        #	if(floor(counter/10) == counter/10) {print(paste(counter,"out of ",length(ncp.list),"familywise Type I error constraints generated"))}
                                        #print(sum(constraint_list[counter,]))
    counter <- counter + 1
}

# print estimated maximum error in computing bivariate normal probabilities using mvtnorm::pmvnorm:
print("Max Error in Multivariate Normal Computation")
print(max_error_prob); save(max_error_prob,file=paste("max_error_prob",task_id,".rdata",sep=""));

# record components of discretized linear program
save(constraint_list,file=paste("A1",task_id,".rdata",sep=""))

if(task_id==1 && dim(constraint_list)[2]!=number_of_variables){print("error in construction of constraints");write(1,file="error_flag")} else{write(dim(constraint_list)[2],f=paste("number_variables.txt"))}

} else if(task_id==max_task_id_for_computing_FWER_constraints+1){# construct objective function vector


objective_function_vector <- rep(0,number_of_variables)
counter <- 1
for(ncp in prior_mean_support)
{
                                        # construct vectors for storing reallocated probabilities for preset decision rectangles
    for(d in decisions){
        for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
            list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- rep(0,length(list.of.rectangles.mtp[[d]][[rprime_counter]]$allowed_actions))
        }
    }

    contribution_to_risk <- c()
    for(r in list.of.rectangles.dec[1:(length(list.of.rectangles.dec)-number_reference_rectangles)]){
        for(d in r$allowed_decisions){
            joint_distribution_parameters <- modified_joint_distribution(prior_component_index=counter,decision=d)
            modified_mean_vector <- joint_distribution_parameters[[1]]
            modified_covariance_matrix <- joint_distribution_parameters[[2]]
            if(d_type[d]==1){
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## Get P_{ncp_vec}(Z \in rectangle r)
                    prob_r <- mvtnorm::pmvnorm(mean=modified_mean_vector,sigma=modified_covariance_matrix,lower=c(r$lower_boundaries,rprime$lower_boundaries),upper=c(r$upper_boundaries,rprime$upper_boundaries),algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                    if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    if(r$preset_decision==0){contribution_to_risk <- c(contribution_to_risk,prob_r*prior_weights[counter]*loss_function_value[d,rprime$allowed_actions])}else{
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles + r$preset_decision_value[d]*prob_r*prior_weights[counter]*loss_function_value[d,rprime$allowed_actions]}
                }
            } else if(d_type[d]==2) {
                ## Get P_{ncp_vec}(Z \in rectangle r)
                 for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## compute overlap in rectangles in terms of z_2:
                    intersection_rectangle_lower_boundaries <- pmax(r$lower_boundaries,rprime$lower_boundaries)
                    intersection_rectangle_upper_boundaries <- pmin(r$upper_boundaries,rprime$upper_boundaries)
                    if(sum(intersection_rectangle_lower_boundaries < intersection_rectangle_upper_boundaries)==2){ #case where intersection rectangle non-empty
                        prob_r <- mvtnorm::pmvnorm(mean=modified_mean_vector,sigma=modified_covariance_matrix,lower=intersection_rectangle_lower_boundaries,upper=intersection_rectangle_upper_boundaries,algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                        if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    } else {# case where rectangle empty
                        prob_r <- 0
                    }
                    if(r$preset_decision==0){contribution_to_risk <- c(contribution_to_risk,prob_r*prior_weights[counter]*loss_function_value[d,rprime$allowed_actions])}else{
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles + r$preset_decision_value[d]*prob_r*prior_weights[counter]*loss_function_value[d,rprime$allowed_actions]}
                }
            } else if(d_type[d]==3) {
                r_lower_boundary_z_2 <- r$lower_boundaries[2]
                r_upper_boundary_z_2 <- r$upper_boundaries[2]
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## compute overlap in rectangles in terms of z_2:
                    effective_z_2_lower_boundary <- max(c(r_lower_boundary_z_2,rprime$lower_boundaries[2]))
                    effective_z_2_upper_boundary <- min(c(r_upper_boundary_z_2,rprime$upper_boundaries[2]))
                    if(effective_z_2_lower_boundary < effective_z_2_upper_boundary){ # case where rectangle non-empty
                        ## Get P_{ncp_vec}(Z \in rectangle r)
                        prob_r <- mvtnorm::pmvnorm(mean=modified_mean_vector,sigma=modified_covariance_matrix,lower=c(r$lower_boundaries[1],effective_z_2_lower_boundary,rprime$lower_boundaries[1]),upper=c(r$upper_boundaries[1],effective_z_2_upper_boundary,rprime$upper_boundaries[1]),algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                        if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    } else {# case where rectangle empty
                        prob_r <- 0
                    }
                    if(r$preset_decision==0){contribution_to_risk <- c(contribution_to_risk,prob_r*prior_weights[counter]*loss_function_value[d,rprime$allowed_actions])}else{
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles + r$preset_decision_value[d]*prob_r*prior_weights[counter]*loss_function_value[d,rprime$allowed_actions]}
                }
            } else if(d_type[d]==4) { #only subpopulation 2 enrolled further
                r_lower_boundary_z_1 <- r$lower_boundaries[1]
                r_upper_boundary_z_1 <- r$upper_boundaries[1]
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## compute overlap in rectangles in terms of z_2:
                    effective_z_1_lower_boundary <- max(c(r_lower_boundary_z_1,rprime$lower_boundaries[1]))
                    effective_z_1_upper_boundary <- min(c(r_upper_boundary_z_1,rprime$upper_boundaries[1]))
                    if(effective_z_1_lower_boundary < effective_z_1_upper_boundary){ # case where rectangle non-empty
                        ## Get P_{ncp_vec}(Z \in rectangle r)
                        prob_r <- mvtnorm::pmvnorm(mean=modified_mean_vector,sigma=modified_covariance_matrix,lower=c(effective_z_1_lower_boundary,r$lower_boundaries[2],rprime$lower_boundaries[2]),upper=c(effective_z_1_upper_boundary,r$upper_boundaries[2],rprime$upper_boundaries[2]),algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                        if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    } else {# case where rectangle empty
                        prob_r <- 0
                    }
                    if(r$preset_decision==0){contribution_to_risk <- c(contribution_to_risk,prob_r*prior_weights[counter]*loss_function_value[d,rprime$allowed_actions])}else{list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles + r$preset_decision_value[d]*prob_r*prior_weights[counter]*loss_function_value[d,rprime$allowed_actions]}
                }
            }
        }
    }

    #Handle reference rectangle contributions
    for(r in list.of.rectangles.dec[(length(list.of.rectangles.dec)+1-number_reference_rectangles):length(list.of.rectangles.dec)]){
        for(d in r$allowed_decisions){
            for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                contribution_to_risk <- c(contribution_to_risk,list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles)
            }
        }
    }

    objective_function_vector <- objective_function_vector + contribution_to_risk
    counter <- counter + 1
}

if(length(objective_function_vector) !=number_of_variables){print("error in construction of constraints");write(2,file="error_flag")} else{save(objective_function_vector,file=paste("c.rdata"))}

# print estimated maximum error in computing bivariate normal probabilities using mvtnorm::pmvnorm:
print("Max Error in Multivariate Normal Computation")
print(max_error_prob); save(max_error_prob,file=paste("max_error_prob",task_id,".rdata",sep=""));


} else if(task_id==max_task_id_for_computing_FWER_constraints+2){
##
## Construct Sparse Constraints corresponding to (26) and (27) in paper
##

# Generate matrix A^2:
# Format: (row,column,value)
# First generate equality constraints representing (24)
equality_constraints_part1 <- array(0,c(length(list.of.rectangles.dec)*number_decisions*number_actions,3))
counter <- 1
constraint_number <- 1
for(r in list.of.rectangles.dec){
    if(r$preset_decision==0){
        for(d in r$allowed_decisions){
                                        # use first rectangle in list.of.rectangles.mtp[[d]] as rprime_d in constraint
            rprime <- list.of.rectangles.mtp[[d]][[1]]
            for(action in rprime$allowed_actions){
                equality_constraints_part1[counter,] <- c(constraint_number,variable_location(r,d,rprime,action),1)
                counter <- counter+1
            }
        }
        constraint_number <- constraint_number+1
    }
}
equality_constraints_part1 <- equality_constraints_part1[1:(counter-1),]

# Next generate equality constraints representing (25)
total_number_stage2_rectangles <- 0
for(d in decisions){total_number_stage2_rectangles <- total_number_stage2_rectangles + length(list.of.rectangles.mtp[[d]])}

equality_constraints_part2 <- array(0,c(length(list.of.rectangles.dec)*total_number_stage2_rectangles*2*number_actions,3))
counter <- 1
for(r in list.of.rectangles.dec){
    if(r$preset_decision==0){
        for(d in r$allowed_decisions){
                                        # use first rectangle in list.of.rectangles.mtp[[d]] as rprime_d in constraint
            rprime_1 <- list.of.rectangles.mtp[[d]][[1]]
            for(rprime_2 in list.of.rectangles.mtp[[d]][2:length(list.of.rectangles.mtp[[d]])]){
                for(action in rprime_1$allowed_actions){
                    equality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime_1,action),1)
                    counter <- counter+1
                }
                for(action in rprime_2$allowed_actions){
                    equality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime_2,action),-1)
                    counter <- counter+1
                }
                constraint_number <- constraint_number+1
            }
        }
    }
}

equality_constraints_part2 <- equality_constraints_part2[1:(counter-1),]

equality_constraints <- rbind(equality_constraints_part1,equality_constraints_part2)

# Check FWER constraints have correct number of variables
if((constraint_number-1)!=(number_equality_constraints_part1+number_equality_constraints_part2)){print("error in construction of constraints");write(3,file="error_flag")} else if(max(equality_constraints[,2]) > number_of_variables){print("ERROR in Generating Equality Constraints: max variable number exceeded");write(4,file="error_flag")} else{
	save(equality_constraints,file=paste("A2.rdata")) #note: these are in sparse matrix format
	}
} else if(task_id==max_task_id_for_computing_FWER_constraints+3){ #Power constraint (uses vector indicator_contribute_to_combined_pop_power)

power_constraint_matrix_H01 <- c()
power_constraint_matrix_H02 <- c()
power_constraint_matrix_H0C <- c()

power_constraint_counter <- 1
for(ncp in power_constraint_list){
    power_constraint_vector_H01 <- power_constraint_vector_H02 <- power_constraint_vector_H0C <- c();
#    indicator_contribute_to_power <- power_constraint_null_hyp_contribution_list[[power_constraint_counter]]
    # construct vectors for storing reallocated probabilities for preset decision rectangles
    for(d in decisions){
        for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
            list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01 <- rep(0,length(list.of.rectangles.mtp[[d]][[rprime_counter]]$allowed_actions));
            list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02 <- rep(0,length(list.of.rectangles.mtp[[d]][[rprime_counter]]$allowed_actions));
            list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C <- rep(0,length(list.of.rectangles.mtp[[d]][[rprime_counter]]$allowed_actions));
        }
    }
    for(r in list.of.rectangles.dec[1:(length(list.of.rectangles.dec)-number_reference_rectangles)]){
        for(d in r$allowed_decisions){
            if(d_type[d]==1){
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## Get P_{ncp_vec}(Z \in rectangle r)
                    prob_r <- mvtnorm::pmvnorm(mean=mean_vector(ncp,d),sigma=covariance_matrix[[d]],lower=c(r$lower_boundaries,rprime$lower_boundaries),upper=c(r$upper_boundaries,rprime$upper_boundaries),algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000)); if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    if(r$preset_decision==0){
                      power_constraint_vector_H01 <- c(power_constraint_vector_H01,prob_r*indicator_contribute_to_H01_power[rprime$allowed_actions]);
                      power_constraint_vector_H02 <- c(power_constraint_vector_H02,prob_r*indicator_contribute_to_H02_power[rprime$allowed_actions]);
                      power_constraint_vector_H0C <- c(power_constraint_vector_H0C,prob_r*indicator_contribute_to_H0C_power[rprime$allowed_actions]);
                    }else{
                      list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01 <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01 + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H01_power[rprime$allowed_actions];
                      list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02 <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02 + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H02_power[rprime$allowed_actions];
                      list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H0C_power[rprime$allowed_actions];
                      }
                }
            } else if(d_type[d]==2) {
                ## Get P_{ncp_vec}(Z \in rectangle r)
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## compute overlap in rectangles in terms of z_2:
                    intersection_rectangle_lower_boundaries <- pmax(r$lower_boundaries,rprime$lower_boundaries)
                    intersection_rectangle_upper_boundaries <- pmin(r$upper_boundaries,rprime$upper_boundaries)
                    if(sum(intersection_rectangle_lower_boundaries < intersection_rectangle_upper_boundaries)==2){ #case where intersection rectangle non-empty
                        prob_r <- mvtnorm::pmvnorm(mean=mean_vector(ncp,d),sigma=covariance_matrix[[d]],lower=intersection_rectangle_lower_boundaries,upper=intersection_rectangle_upper_boundaries,algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                        if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    } else {# case where rectangle empty
                        prob_r <- 0
                    }
                    if(r$preset_decision==0){
                      power_constraint_vector_H01 <- c(power_constraint_vector_H01,prob_r*indicator_contribute_to_H01_power[rprime$allowed_actions]);
                      power_constraint_vector_H02 <- c(power_constraint_vector_H02,prob_r*indicator_contribute_to_H02_power[rprime$allowed_actions]);
                      power_constraint_vector_H0C <- c(power_constraint_vector_H0C,prob_r*indicator_contribute_to_H0C_power[rprime$allowed_actions]);} else{
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01 <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01 + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H01_power[rprime$allowed_actions];
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02 <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02 + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H02_power[rprime$allowed_actions];
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H0C_power[rprime$allowed_actions];
                      }
                }
            } else if(d_type[d]==3) {
                r_lower_boundary_z_2 <- r$lower_boundaries[2]
                r_upper_boundary_z_2 <- r$upper_boundaries[2]
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## compute overlap in rectangles in terms of z_2:
                    effective_z_2_lower_boundary <- max(c(r_lower_boundary_z_2,rprime$lower_boundaries[2]))
                    effective_z_2_upper_boundary <- min(c(r_upper_boundary_z_2,rprime$upper_boundaries[2]))
                    if(effective_z_2_lower_boundary < effective_z_2_upper_boundary){ # case where rectangle non-empty
                        ## Get P_{ncp_vec}(Z \in rectangle r)
                        prob_r <- mvtnorm::pmvnorm(mean=mean_vector(ncp,d),sigma=covariance_matrix[[d]],lower=c(r$lower_boundaries[1],effective_z_2_lower_boundary,rprime$lower_boundaries[1]),upper=c(r$upper_boundaries[1],effective_z_2_upper_boundary,rprime$upper_boundaries[1]),algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                        if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    } else {# case where rectangle empty
                        prob_r <- 0
                    }
                    if(r$preset_decision==0){
                      power_constraint_vector_H01 <- c(power_constraint_vector_H01,prob_r*indicator_contribute_to_H01_power[rprime$allowed_actions]);
                      power_constraint_vector_H02 <- c(power_constraint_vector_H02,prob_r*indicator_contribute_to_H02_power[rprime$allowed_actions]);
                      power_constraint_vector_H0C <- c(power_constraint_vector_H0C,prob_r*indicator_contribute_to_H0C_power[rprime$allowed_actions]);} else{
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01 <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01 + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H01_power[rprime$allowed_actions];
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02 <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02 + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H02_power[rprime$allowed_actions];
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H0C_power[rprime$allowed_actions];
                      }
                }
            } else if(d_type[d]==4) { #only subpopulation 2 enrolled further
                r_lower_boundary_z_1 <- r$lower_boundaries[1]
                r_upper_boundary_z_1 <- r$upper_boundaries[1]
                for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    rprime <- list.of.rectangles.mtp[[d]][[rprime_counter]]
                    ## compute overlap in rectangles in terms of z_2:
                    effective_z_1_lower_boundary <- max(c(r_lower_boundary_z_1,rprime$lower_boundaries[1]))
                    effective_z_1_upper_boundary <- min(c(r_upper_boundary_z_1,rprime$upper_boundaries[1]))
                    if(effective_z_1_lower_boundary < effective_z_1_upper_boundary){ # case where rectangle non-empty
                        ## Get P_{ncp_vec}(Z \in rectangle r)
                        prob_r <- mvtnorm::pmvnorm(mean=mean_vector(ncp,d),sigma=covariance_matrix[[d]],lower=c(effective_z_1_lower_boundary,r$lower_boundaries[2],rprime$lower_boundaries[2]),upper=c(effective_z_1_upper_boundary,r$upper_boundaries[2],rprime$upper_boundaries[2]),algorithm=mvtnorm::GenzBretz(abseps = 0.000000001,maxpts=100000))
                        if(attr(prob_r,"error") > max_error_prob){max_error_prob <- attr(prob_r,"error")}
                    } else {# case where rectangle empty
                        prob_r <- 0
                    }
                    if(r$preset_decision==0){
                      power_constraint_vector_H01 <- c(power_constraint_vector_H01,prob_r*indicator_contribute_to_H01_power[rprime$allowed_actions]);
                      power_constraint_vector_H02 <- c(power_constraint_vector_H02,prob_r*indicator_contribute_to_H02_power[rprime$allowed_actions]);
                      power_constraint_vector_H0C <- c(power_constraint_vector_H0C,prob_r*indicator_contribute_to_H0C_power[rprime$allowed_actions]);} else{
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01 <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01 + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H01_power[rprime$allowed_actions];
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02 <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02 + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H02_power[rprime$allowed_actions];
                        list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C <- list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C + r$preset_decision_value[d]*prob_r*indicator_contribute_to_H0C_power[rprime$allowed_actions];
                      }
                }
            }
        }
    }
    #Handle reference rectangle contributions
    for(r in list.of.rectangles.dec[(length(list.of.rectangles.dec)+1-number_reference_rectangles):length(list.of.rectangles.dec)]){
        for(d in r$allowed_decisions){
            for(rprime_counter in 1:length(list.of.rectangles.mtp[[d]])){
                    power_constraint_vector_H01 <- c(power_constraint_vector_H01,list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H01);
                    power_constraint_vector_H02 <- c(power_constraint_vector_H02,list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H02);
                    power_constraint_vector_H0C <- c(power_constraint_vector_H0C,list.of.rectangles.mtp[[d]][[rprime_counter]]$reallocated_probability_from_preset_decision_rectangles_H0C)
            }
        }
    }
    power_constraint_matrix_H01 <- rbind(power_constraint_matrix_H01,power_constraint_vector_H01);
    power_constraint_matrix_H02 <- rbind(power_constraint_matrix_H02,power_constraint_vector_H02);
    power_constraint_matrix_H0C <- rbind(power_constraint_matrix_H0C,power_constraint_vector_H0C);
    power_constraint_counter <- power_constraint_counter+1
}
save(power_constraint_matrix_H01,power_constraint_matrix_H02,power_constraint_matrix_H0C,file=paste("A3.rdata"))

# print estimated maximum error in computing bivariate normal probabilities using mvtnorm::pmvnorm:
print("Max Error in Multivariate Normal Computation")
print(max_error_prob); save(max_error_prob,file=paste("max_error_prob",task_id,".rdata",sep=""));

} else if(task_id==max_task_id_for_computing_FWER_constraints+4){
# Generate sparse inequality constraints (32) and (33) in Section 4.1 that restrict multiple testing procedure not depend on stage 1 statistics (given cumulative statistics Z^C and decision d)

additional_inequality_constraints_part1 <- array(0,c(length(list.of.rectangles.dec)*sum(number_rectangles_by_actions_for_decision)*(2+number_actions*(number_decisions-1)),3))
constraint_number <- 1
counter <- 1
for(r in list.of.rectangles.dec[1:(length(list.of.rectangles.dec)-number_reference_rectangles)]){
    if(r$preset_decision==0){
        for(d in r$allowed_decisions){
            ## first set r_reference to be reference rectangle corresponding to r
            r_reference<-list.of.rectangles.dec[[length(list.of.rectangles.dec)-number_reference_rectangles+d]]
            for(rprime in list.of.rectangles.mtp[[d]]){
                for(action in rprime$allowed_actions){
                                        # Constraint (33):
                    additional_inequality_constraints_part1[counter,] <- c(constraint_number,variable_location(r_reference,d,rprime,action),1)
                    counter <- counter +1
                    additional_inequality_constraints_part1[counter,] <- c(constraint_number,variable_location(r,d,rprime,action),-1)
                    counter <- counter +1
                    r_tilde_position <- r$position_offset
                    for(d_tilde in r$allowed_decisions){
                        if(d_tilde != d){
                            rprime_tilde <- list.of.rectangles.mtp[[d_tilde]][[1]]
                            for(action in rprime_tilde$allowed_actions){
                                additional_inequality_constraints_part1[counter,] <- c(constraint_number,variable_location(r,d_tilde,rprime_tilde,action),-1)
                                counter <- counter +1
                            }
                        }
                    }
                    constraint_number <- constraint_number+1
                }
            }
        }

    }
}

if(counter>1){
# previous output in case nothing generated: {additional_inequality_constraints_part1 <- c(); save(additional_inequality_constraints_part1,file=paste("Inequality_Constraints_to_Restrict_MTP_to_Sufficient_Statistics.rdata"))}

additional_inequality_constraints_part1 <- additional_inequality_constraints_part1[1:(counter-1),]
if(max(additional_inequality_constraints_part1[,2]) > number_of_variables){print("ERROR in Generating Inequality_Constraints_to_Restrict_MTP_to_Sufficient_Statistics");write(5,file="error_flag")} else{
save(additional_inequality_constraints_part1,file=paste("Inequality_Constraints_to_Restrict_MTP_to_Sufficient_Statistics.rdata"))}
}

} else if(task_id==max_task_id_for_computing_FWER_constraints+5){                                        # Generate sparse inequality constraints (34) and (35) in Section 4.1 that represent monotonicity constraints in the multiple testing procedure

additional_inequality_constraints_part2 <- array(0,c(sum(number_rectangles_by_actions_for_decision)*2*number_actions*4,3))
constraint_number <- 1
counter <- 1
for(r in list.of.rectangles.dec[(length(list.of.rectangles.dec)-number_reference_rectangles+1):length(list.of.rectangles.dec)]){## only need to set for reference rectangles
    for(d in r$allowed_decisions){
        for(rprime in list.of.rectangles.mtp[[d]]){
            if((sum(rprime$allowed_actions==2)>0 || sum(rprime$allowed_actions==5)>0 || sum(rprime$allowed_actions==7)>0) && !is.null(rprime$right_neighbor)){ # if possitlbe to reject H01 and exists a right neighbor
                for(action in c(2,5,7)){
                    if(sum(rprime$allowed_actions==action)>0){
                        additional_inequality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime,action),1)
                        counter <- counter +1
                    }
                }
                rprime_right_neighbor <-  list.of.rectangles.mtp[[d]][[rprime$right_neighbor]]
                for(action in c(2,5,7)){
                    if(sum(rprime_right_neighbor$allowed_actions==action)>0){
                        additional_inequality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime_right_neighbor,action),-1)
                        counter <- counter +1
                    }
                }
                constraint_number <- constraint_number+1
            }
            if((sum(rprime$allowed_actions==3)>0 || sum(rprime$allowed_actions==6)>0 || sum(rprime$allowed_actions==7)>0) && !is.null(rprime$upper_neighbor)){ # if possitlbe to reject H02 and exists a upper neighbor
                for(action in c(3,6,7)){
                    if(sum(rprime$allowed_actions==action)>0){
                        additional_inequality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime,action),1)
                        counter <- counter +1
                    }
                }
                rprime_upper_neighbor <-  list.of.rectangles.mtp[[d]][[rprime$upper_neighbor]]
                for(action in c(3,6,7)){
                    if(sum(rprime_upper_neighbor$allowed_actions==action)>0){
                        additional_inequality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime_upper_neighbor,action),-1)
                        counter <- counter +1
                    }
                }
                constraint_number <- constraint_number+1
            }

            if((sum(rprime$allowed_actions==4)>0 || sum(rprime$allowed_actions==5)>0 || sum(rprime$allowed_actions==6)>0 || sum(rprime$allowed_actions==7)>0) && (!is.null(rprime$right_neighbor) || !!is.null(rprime$upper_neighbor))){# if possible to reject H0C and exists a right or upper neighbor
                if(!is.null(rprime$right_neighbor)){
                    for(action in c(4,5,6,7)){
                        if(sum(rprime$allowed_actions==action)>0){
                            additional_inequality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime,action),1)
                            counter <- counter +1
                        }
                    }
                    rprime_right_neighbor <-  list.of.rectangles.mtp[[d]][[rprime$right_neighbor]]
                    for(action in c(4,5,6,7)){
                        if(sum(rprime_right_neighbor$allowed_actions==action)>0){
                            additional_inequality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime_right_neighbor,action),-1)
                            counter <- counter +1
                        }
                    }
                    constraint_number <- constraint_number+1
                }

                if(!is.null(rprime$upper_neighbor)){
                    for(action in c(4,5,6,7)){
                        if(sum(rprime$allowed_actions==action)>0){
                            additional_inequality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime,action),1)
                            counter <- counter +1
                        }
                    }
                    rprime_upper_neighbor <-  list.of.rectangles.mtp[[d]][[rprime$upper_neighbor]]
                    for(action in c(4,5,6,7)){
                        if(sum(rprime_upper_neighbor$allowed_actions==action)>0){
                            additional_inequality_constraints_part2[counter,] <- c(constraint_number,variable_location(r,d,rprime_upper_neighbor,action),-1)
                            counter <- counter +1
                        }
                    }
                    constraint_number <- constraint_number+1
                }
            }
        }
    }
}
additional_inequality_constraints_part2 <- additional_inequality_constraints_part2[1:(counter-1),]
if(max(additional_inequality_constraints_part2[,2]) > number_of_variables){print("ERROR in Generating Inequality_Constraints_to_set_monotonicity_in_hypotheses_rejected");write(6,file="error_flag")} else{
save(additional_inequality_constraints_part2,file=paste("Inequality_Constraints_to_set_monotonicity_in_hypotheses_rejected.rdata"))}
}
}

#
# Generate Linear program in parallel
#

number_jobs <- ceiling(length(ncp.list)/constraints_per_A1_file)+6
parallel::mclapply(c((number_jobs-5):number_jobs,1:(number_jobs-6)),generate_LP,mc.cores=number.cores) # order of jobs puts computation of power constraints and objective function first since they take longer to compute

# Convert linear program to matlab format

R.matlab::writeMat("alphaValue.mat",alphaValue=total.alpha)

number_A1_files <- scan("number_A1_files.txt")
A1 = numeric(0)
for(i in 1:number_A1_files){
        tmp = load(paste("A1",i,".rdata",sep=""))
        tmp = constraint_list
        R.matlab::writeMat(paste("A1",i,".mat",sep=""),tmp=tmp)
}

tmp = load("A2.rdata")
tmp = equality_constraints
R.matlab::writeMat("A2.mat",A2 =tmp)

tmp = load("A3.rdata")
A3  = rbind(power_constraint_matrix_H01, power_constraint_matrix_H02,power_constraint_matrix_H0C)
R.matlab::writeMat("A3.mat",A3=A3)

tmp = load("power_constraints.rdata")
a3  = power.constraints
R.matlab::writeMat("a3.mat",a3=a3)

try((tmp = load("Inequality_Constraints_to_Restrict_MTP_to_Sufficient_Statistics.rdata")),silent=TRUE)
if(any(ls()=="additional_inequality_constraints_part1")){
tmp = additional_inequality_constraints_part1
col = read.table("number_variables.txt")
col = col$V1
R.matlab::writeMat("A4.mat",A4=tmp)
rm(additional_inequality_constraints_part1)
R.matlab::writeMat("a4status.mat",a4status=1)
} else {R.matlab::writeMat("a4status.mat",a4status=0)}

tmp11 = read.table("number_equality_constraints_of_first_type.txt")
tmp11 = tmp11$V1
R.matlab::writeMat("a21.mat",a21 = tmp11)

R.matlab::writeMat("iteration.mat",iteration = LP.iteration)

tmp = load("Inequality_Constraints_to_set_monotonicity_in_hypotheses_rejected.rdata")
tmp = additional_inequality_constraints_part2

R.matlab::writeMat("A5.mat",A5=tmp)

tmp = load("c.rdata")
obj = objective_function_vector
R.matlab::writeMat("cc.mat",cc = obj)

#
# Solve linear program by call to cplex via matlab
#

system('matlab -nojvm -r "siterprl()" > output_LP_solver')

#
# Extract results from linear program solver and examine whether feasible solution was found
#

sln = R.matlab::readMat(paste("sln2M",LP.iteration,".mat",sep=""))
save(sln,file=paste("sln2M",LP.iteration,".rdata",sep=""))
input_parameters <- as.list(environment())
print(paste("Adaptive Design Optimization Completed. Optimal design is stored in the file: optimized_design.rdata"))
save(input_parameters,list.of.rectangles.dec,list.of.rectangles.mtp,ncp.list,sln,file=paste("optimized.design",LP.iteration,".rdata"))

if(sln$status==1 || sln$status==5){
  print(paste("Feasible Solution was Found and Optimal Expected Sample Size is",sln$val))
} else {print("Problem was Infeasible")}

if(round.each.decision.rectangle.to.integer){## If Final iteration, round solution and save; only does this if decision rule was rounded and set to be deterministic
print("Fraction of solution components with integral value solutions")
print(sum(sln$z>1-1e-10 | sln$z<10e-10)/length(sln$z))

postscript(paste("rejection_regions.eps"),height=8,horizontal=FALSE,onefile=FALSE,width=8)
plot(0,type="n",xaxt="n",yaxt="n",xlim=c(-2.78,2.78),ylim=c(-2.78,2.78),main="Rejection Regions",xlab=expression(paste(Z[1])),ylab=expression(paste(Z[2])),cex.lab=2,
 cex.axis=2, cex.main=2, cex.sub=2)

z_solution <- sln$z
axis(1,at=seq(-3,3,by=1),labels=-3:3,cex.axis=2)
axis(2,at=seq(-3,3,by=1),labels=-3:3,cex.axis=2)
z_rounded <- rep(0,length(z_solution))
for(d_plot in decisions){
	rounding_threshold_H01 <- rounding_threshold <- 0.9
	rounding_threshold_H02 <- rounding_threshold <- 0.9
	rounding_threshold_H0C <- rounding_threshold <- 0.9
	r_reference_index <- length(list.of.rectangles.dec) - number_reference_rectangles + d_plot
	r_reference <- list.of.rectangles.dec[[r_reference_index]]
	for(d in r_reference$allowed_decisions){
	      for(rprime in list.of.rectangles.mtp[[d]]){
		 variable_start_position <- variable_location(r_reference,d,rprime,rprime$allowed_actions[1])
		 variable_end_position <- variable_location(r_reference,d,rprime,rprime$allowed_actions[length(rprime$allowed_actions)])
		 action_indicator <-z_solution[variable_start_position:variable_end_position]
             
		 if( sum(action_indicator[c(2,5,7)])>rounding_threshold_H01){H01_reject <- 1}else{H01_reject <- 0}
		 if( sum(action_indicator[c(3,6,7)])>rounding_threshold_H02){H02_reject <- 1}else{H02_reject <- 0}
	    	 if( ((H01_reject && H02_reject) || sum(action_indicator[c(4,5,6,7)])>rounding_threshold_H0C)){H0C_reject <- 1}else{H0C_reject <- 0}
					   
          	 col_value <- ifelse((!H01_reject) && (!H02_reject) && (!H0C_reject),1,
		              ifelse((H01_reject) && (!H02_reject) && (!H0C_reject),2,  
            		      ifelse((!H01_reject) && (H02_reject) && (!H0C_reject),3,
          		      ifelse((!H01_reject) && (!H02_reject) && (H0C_reject),4,  
             		      ifelse((H01_reject) && (!H02_reject) && (H0C_reject),5,
            		      ifelse((!H01_reject) && (H02_reject) && (H0C_reject),6,  
             		      ifelse((H01_reject) && (H02_reject) && (H0C_reject),7,8)))))))
                 
                 z_rounded[variable_start_position:variable_end_position] <- rep(0,7)
                 z_rounded[variable_start_position+(col_value-1)] <- 1
	         rect(max(rprime$lower_boundaries[1]-tau,-10),max(rprime$lower_boundaries[2]-tau,-10),min(rprime$upper_boundaries[1]+tau,10),min(rprime$upper_boundaries[2]+tau,10),col=col_value-1,border=NA) 
             }
        }
}
par(las=0)
dev.off()	

save(z_rounded,file="z_rounded.rdata")

max_FWER <- 0
for(task_id in 1:length(ncp.list)){
   load(file=paste("A1",task_id,".rdata",sep=""))
   #print(ncp.list[[which((constraint_list %*% sln$z) ==max(constraint_list %*% sln$z))]])
   fwer_candidates <- constraint_list %*% z_rounded
   if(max(fwer_candidates)>max_FWER){
	max_FWER <- max(max_FWER,max(constraint_list %*% z_rounded))
	print(constraint_list %*% z_rounded)
	#print(ncp.list[(1+(counter-1)*6):(counter*6)])
   }
}
print("Maximum familywise Type I error rate among Type I error constraints")
print(max_FWER)

load("A3.rdata")
print("User defined power constraints (desired power); each row corresponds to a scenario and columns correspond to H01, H02, H0C desired power")
print("power.constraints")
print("Power achieved for each null hypothesis under each power constraint scenario (row)")
print(cbind(power_constraint_matrix_H01 %*% z_rounded,power_constraint_matrix_H02 %*% z_rounded,power_constraint_matrix_H03 %*% z_rounded))

load("c.rdata")
print("Objective function value")
print(objective_function_vector %*% z_rounded)
}
					      
# Clean up files used to specify LP
system('rm A*.rdata')
system('rm A*.mat')
system('rm Inequality_Constraints_to_Restrict_MTP_to_Sufficient_Statistics.rdata')
system('rm Inequality_Constraints_to_set_monotonicity_in_hypotheses_rejected.rdata')
system('rm number_equality_constraints_of_first_type.txt')
system('rm number_equality_constraints_of_second_type.txt')
system('rm number_A1_constraints.txt')
system('rm number_A1_files.txt')
system('rm power_constraints.rdata')

  #Clean up files
  system('rm max_error_prob*')
}