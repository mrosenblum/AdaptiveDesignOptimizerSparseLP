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
#' @param LP.iteration positive integer used in file name to store output; can be used to avoid overwriting previous computations
#' @param round.each.decision.rectangle.to.integer TRUE/FALSE indicator of whether decision probabilities encoded in list.of.rectangles.dec should be rounded to integer values
#' @param set.rectangles.with.identically.valued.neighbors.and.split.others  TRUE/FALSE indicator of whether decision probabilities encoded in list.of.rectangles.dec should be modified for use in next iteration
#' @param sln solution to linear program computed previously
#' @return 4 element list containing optimized designs from four classes (with increasing complexity):
#' @section Output
#' The software computes and optimized design saved as "optimized_design.rdata" and the corresponding expected sample size is
#' saved as "optimized_design_expected_sample_size.rdata".
#' @examples
#' #For demonstration purposes, the examples below use a coarse discretization.
#' optimize_design(discretization.parameter=c(3,3,1),number.cores=1)
#' @export
refine_decision_rectangles <- function(subpopulation.1.proportion=0.5,
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
		set.rectangles.with.identically.valued.neighbors.and.split.others=TRUE,
		sln
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

  # construct list of rectangles in set R
  ## List of rectangles defining decision boundaries
  ## Each rectangle is encoded by c(lower left (x,y) coordinates, upper right (x,y) coordinates)

if(!is.null(list.of.rectangles.dec)){list.of.rectangles.dec <- list.of.rectangles.dec} else {
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
   #save(list.of.rectangles.dec,file=paste("list.of.rectangles.dec",LP.iteration,".rdata",sep=""))
   #save(ncp.list,file=paste("ncp.list",LP.iteration,".rdata",sep=""))
}

## Set multiple testing procedure rectangle partition

stage_2_rectangle_offset_value <- 0
list_of_rectangles_mtp1 <- list()
for(x in seq(-w1,w1,by=tau_mtp))
{
	#Skip lower left quadrant
	for(y in seq(ifelse(x<0,0,-w2),w2,by=tau_mtp))
	{
		list_of_rectangles_mtp1<- c(list_of_rectangles_mtp1,list(list(lower_boundaries=c(x,y),upper_boundaries=c(x+tau_mtp,y+tau_mtp),allowed_actions=actions,stage_2_rectangle_offset=stage_2_rectangle_offset_value)))
		stage_2_rectangle_offset_value <- stage_2_rectangle_offset_value + length(actions)
	}
}
## Only one large rectangle for lower left quadrant
list_of_rectangles_mtp1<- c(list_of_rectangles_mtp1,list(list(lower_boundaries=c(-Inf,-Inf),upper_boundaries=c(0,0),allowed_actions=actions,stage_2_rectangle_offset=stage_2_rectangle_offset_value)))
stage_2_rectangle_offset_value <- stage_2_rectangle_offset_value + length(actions)

number_rectangles_by_actions_for_decision_type1  <- stage_2_rectangle_offset_value

## Set right neighbors
counter_for_rprime <- 1
for(rprime in list_of_rectangles_mtp1)
{
	count_value <- 1
	while(count_value <= length(list_of_rectangles_mtp1) && (!(rprime$upper_boundaries[1]== list_of_rectangles_mtp1[[count_value]]$lower_boundaries[1] && rprime$lower_boundaries[2]== list_of_rectangles_mtp1[[count_value]]$lower_boundaries[2] && rprime$upper_boundaries[2]== list_of_rectangles_mtp1[[count_value]]$upper_boundaries[2]))){count_value <- count_value +1}
	if(count_value <= length(list_of_rectangles_mtp1)){list_of_rectangles_mtp1[[counter_for_rprime]]$right_neighbor <- count_value}
	counter_for_rprime <- counter_for_rprime +1
}
## Set upper neighbors
counter_for_rprime <- 1
for(rprime in list_of_rectangles_mtp1)
{
	count_value <- 1
	while(count_value <= length(list_of_rectangles_mtp1) && (!(rprime$upper_boundaries[2]== list_of_rectangles_mtp1[[count_value]]$lower_boundaries[2] && rprime$lower_boundaries[1]== list_of_rectangles_mtp1[[count_value]]$lower_boundaries[1] && rprime$upper_boundaries[1]== list_of_rectangles_mtp1[[count_value]]$upper_boundaries[1]))){count_value <- count_value +1}
	if(count_value <= length(list_of_rectangles_mtp1)){list_of_rectangles_mtp1[[counter_for_rprime]]$upper_neighbor <- count_value}
	counter_for_rprime <- counter_for_rprime +1
}

##Set list_of_rectangles_mtp2
stage_2_rectangle_offset_value <- 0
list_of_rectangles_mtp2 <- list()
for(x in seq(-w1,w1,by=tau)){
    for(y in seq(-w2,w2,by=tau)){
        list_of_rectangles_mtp2<- c(list_of_rectangles_mtp2,list(list(lower_boundaries=c(x,y),upper_boundaries=c(x+tau,y+tau),allowed_actions=actions,stage_2_rectangle_offset=stage_2_rectangle_offset_value)))
        stage_2_rectangle_offset_value <- stage_2_rectangle_offset_value + length(actions)
    }
}
## IF WANT TO TEMPORARILY SET MTP2 TO SINGLE LARGE RECTANGLE, replace above loops by:
#list_of_rectangles_mtp2<- c(list_of_rectangles_mtp2,list(list(lower_boundaries=c(-Inf,-Inf),upper_boundaries=c(Inf,Inf),allowed_actions=actions,stage_2_rectangle_offset=stage_2_rectangle_offset_value)))
#stage_2_rectangle_offset_value <- stage_2_rectangle_offset_value + length(actions)

number_rectangles_by_actions_for_decision_type2  <- stage_2_rectangle_offset_value

##Set stage 2 rectangle sets for each d_type:
#Set stage 2 rectangles sets to be the same for d_type 1, 3, and 4. (Allow different for d_type 2)
list_of_rectangles_mtp3 <- list_of_rectangles_mtp4 <- list_of_rectangles_mtp1;
number_rectangles_by_actions_for_decision_type3 <- number_rectangles_by_actions_for_decision_type4 <- number_rectangles_by_actions_for_decision_type1
list_of_rectangles_mtp <- list()
number_rectangles_by_actions_for_decision <- c()
for(d in decisions){
  list_of_rectangles_mtp <- c(list_of_rectangles_mtp,list(switch(d_type[d],list_of_rectangles_mtp1,list_of_rectangles_mtp2,list_of_rectangles_mtp3,list_of_rectangles_mtp4)));
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
            number_equality_constraints_part2<-number_equality_constraints_part2+(length(list_of_rectangles_mtp[[d]])-1)
	}
    }
}

print("number of variables")
print(number_of_variables)
print("number of familywise Type I error constraints")

if(sln$status==(-9)){return("Linear Program was Infeasible; Please Try Again e.g., With Greater Sample Sizes")} else if(sln$status==1 || sln$status==5){
print(paste("Linear Program at Iteration ",LP.iteration," Solved. Now Evaluating New Linear Program with Finer Discretization of Decision Regions"))} else{return(paste("Error in linear program; see solver output: see sln2M",LP.iteration,".rdata"))}
print("Fraction of solution components with integral value solutions")
print(sum(sln$z>1-1e-10 | sln$z<10e-10)/length(sln$z))

z_solution <- sln$z
list_of_rectangles_dec <- list.of.rectangles.dec

#
# Either (i) round all decision region rectangles to be integer valued or (ii) set rectangles surrounded by identically valued rectanges and split others.
#

if(round.each.decision.rectangle.to.integer){
  generate_list_of_rectangles_dec_with_decision_probs <- function(){
    list_of_rectangles_dec_with_decision_probs <- list_of_rectangles_dec[1:(length(list_of_rectangles_dec)-number_reference_rectangles)]
    counter <- 1
    for(r_counter in 1:(length(list_of_rectangles_dec)-number_reference_rectangles)){
      r <- list_of_rectangles_dec_with_decision_probs[[r_counter]]
      if(r$preset_decision==1){
        list_of_rectangles_dec_with_decision_probs[[counter]]$d_probs <- r$preset_decision_value
      } else{
        for(d in decisions){
          if(sum(d==r$allowed_decisions)>0){
            # use first rectangle in list_of_rectangles_mtp[[d]] as rprime_d in constraint
            rprime <- list_of_rectangles_mtp[[d]][[1]]
            variable_start_position <- variable_location(r,d,rprime,rprime$allowed_actions[1])
            variable_end_position <- variable_location(r,d,rprime,rprime$allowed_actions[length(rprime$allowed_actions)])
            list_of_rectangles_dec_with_decision_probs[[counter]]$d_probs[[d]] <-sum(z_solution[variable_start_position:variable_end_position])
          } else{list_of_rectangles_dec_with_decision_probs[[counter]]$d_probs[[d]] <- 0}
        }
      }
      counter<- counter+1
    }
    return(list_of_rectangles_dec_with_decision_probs )
  }

  list_of_rectangles_dec_with_decision_probs <- generate_list_of_rectangles_dec_with_decision_probs()

  ## Code to check neighbor_check works
  #for(r in list_of_rectangles_dec_with_decision_probs){
  #plot(0,type="n",xlim=c(-10,10),ylim=c(-10,10),main=paste("Decision Rule for Stage 2 enrollment"),xlab="Z_stage_1_subpop_1",ylab="Z_stage_1_subpop_2",cex=2)
  #for(r2 in list_of_rectangles_dec_with_decision_probs){
  #	if(neighbor_check(r,r2) && r$lower_boundaries[1]>=(-3) && r$lower_boundaries[2]>=(-3) && r$upper_boundaries[1]<=3 && r$upper_boundaries[2]<=3){
  #		col_value <- 2
  #	rect(max(r$lower_boundaries[1],-10),max(r$lower_boundaries[2],-10),min(r$upper_boundaries[1],10),min(r$upper_boundaries[2],10),col=2)
  #	rect(max(r2$lower_boundaries[1],-10),max(r2$lower_boundaries[2],-10),min(r2$upper_boundaries[1],10),min(r2$upper_boundaries[2],10),col=3)
  #	browser()
  #	}
  #}}


  for(counter in 1:length(list_of_rectangles_dec_with_decision_probs)){
    list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision <- 1
    list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_value <- list_of_rectangles_dec_with_decision_probs[[counter]]$d_probs
  }

  ## Round all values
  for(counter in 1:length(list_of_rectangles_dec_with_decision_probs)){
    #if(r$preset_decision==1){
    #if(max(list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_value)>=0.8){
    d <- which(list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_value==max(list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_value))
    list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_value <- rep(0,length(decisions))
    list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_value[d] <- 1
    list_of_rectangles_dec_with_decision_probs[[counter]]$allowed_decisions <- c(d)
    list_of_rectangles_dec_with_decision_probs[[counter]]$d_probs <- rep(0,length(decisions))
    list_of_rectangles_dec_with_decision_probs[[counter]]$d_probs[d] <- 1
    #    rect(max(r$lower_boundaries[1]-tau,-10),max(r$lower_boundaries[2]-tau,-10),min(r$upper_boundaries[1]+tau,10),min(r$upper_boundaries[2]+tau,10),col=d) #, border=NA)
  }


  merge_rectangles <- function(rectangle_list_to_be_merged){
    lgth = length(rectangle_list_to_be_merged)
    xl   = rep(0,lgth)
    xu   = rep(0,lgth)
    yl   = rep(0,lgth)
    yu   = rep(0,lgth)
    dec  = rep(0,lgth)
    d_probs = list()
    allowed_decisions = list()

    for (i in 1:lgth)
    {
      xl[i]  = rectangle_list_to_be_merged[[i]]$lower_boundaries[1]
      xu[i]  = rectangle_list_to_be_merged[[i]]$upper_boundaries[1]
      yl[i]  = rectangle_list_to_be_merged[[i]]$lower_boundaries[2]
      yu[i]  = rectangle_list_to_be_merged[[i]]$upper_boundaries[2]
      dec[i] = rectangle_list_to_be_merged[[i]]$allowed_decisions
      d_probs[[i]] = rectangle_list_to_be_merged[[i]]$d_probs
      allowed_decisions[[i]] = rectangle_list_to_be_merged[[i]]$allowed_decisions
    }

    nxl  = c()
    nxu  = c()
    nyl  = c()
    nyu  = c()
    ndec = c()

    nd_probs = list()
    nallowed_decisions = list()

    done = rep(0,lgth)

    oxl  = order(xl)
    xl   = xl[oxl]
    xu   = xu[oxl]
    yl   = yl[oxl]
    yu   = yu[oxl]
    dec  = dec[oxl]
    d_probs = d_probs[oxl]
    allowed_decisions = allowed_decisions[oxl]


    w = 1
    for (i in 1:lgth)
    {
      if (done[i]==0)
      {
        nxl[w]  = xl[i]
        nxu[w]  = xu[i]
        nyl[w]  = yl[i]
        nyu[w]  = yu[i]
        ndec[w] = dec[i]
        if (length(d_probs)>=i)
        {
          if (!is.null(d_probs[[i]]))
          {
            nd_probs[[w]] = d_probs[[i]]
          }
        }
        nallowed_decisions[[w]] = allowed_decisions[[i]]
        if (i < lgth)
        {
          for (j in (i+1):lgth)
          {
            if ((yl[j]==yl[i])&(yu[j]==yu[i])&(xl[j]==nxu[w])&(abs(dec[j]-dec[i])<1e-3)&&(dec[j]>0.1))
            {
              nxu[w]  = xu[j]
              done[j] = 1
            }
          }
        }
        w = w + 1
      }
    }
    n = length(nxu)
    combined_rectangles= list()
    for (i in 1:n)
    {
      tmp = list()
      tmp$lower_boundaries       = c(nxl[i],nyl[i])
      tmp$upper_boundaries       = c(nxu[i],nyu[i])
      if (length(nd_probs)>=i)
      {
        if (!is.null(nd_probs[[i]]))
        {
          tmp$d_probs                = nd_probs[[i]]
        }
      }
      tmp$allowed_decisions      = nallowed_decisions[[i]]
      if (ndec[i]>0){
        tmp$preset_decision = ndec[i]
      }
      combined_rectangles[[i]]   = tmp
    }
    return(combined_rectangles)

  }

  list_of_rectangles_dec_with_decision_probs_merged <- merge_rectangles(list_of_rectangles_dec_with_decision_probs)

  for(counter in 1:length(list_of_rectangles_dec_with_decision_probs_merged)){
    list_of_rectangles_dec_with_decision_probs_merged[[counter]]$preset_decision_value <- list_of_rectangles_dec_with_decision_probs_merged[[counter]]$d_probs
  }

#save(list.of.rectangles.dec.with.decision.probs.merged,file=paste("list.of.rectangles.dec.rounded.to.integer.rdata",sep=""))
list.of.rectangles.dec <- list_of_rectangles_dec_with_decision_probs_merged

## To view results of rounding and merging
postscript(paste("decision_rectangles.eps"),height=8,horizontal=FALSE,onefile=FALSE,width=8)
plot(0,type="n",xlim=c(-8,8),ylim=c(-8,8),main=paste("Decision Rule for Stage 2 enrollment"),xlab="Z_stage_1_subpop_1",ylab="Z_stage_1_subpop_2",cex=2)

for(counter in 1:length(list.of.rectangles.dec))
{
  r <- list.of.rectangles.dec[[counter]]
  if(r$preset_decision>0){color_value <- r$preset_decision}else{color_value <- 5}
  rect(max(r$lower_boundaries[1],-10),max(r$lower_boundaries[2],-10),min(r$upper_boundaries[1],10),min(r$upper_boundaries[2],10),col=color_value)
}
dev.off()

#
## run set and split to construct new list.of.rectangles.dec, save it, and proceed
#

} else if(set.rectangles.with.identically.valued.neighbors.and.split.others){
   #source("set_fixed_rectangles_and_split_border_rectangles.R")
  generate_list_of_rectangles_dec_with_decision_probs <- function(){
    list_of_rectangles_dec_with_decision_probs <- list_of_rectangles_dec[1:(length(list_of_rectangles_dec)-number_reference_rectangles)]
    counter <- 1
    for(r_counter in 1:(length(list_of_rectangles_dec)-number_reference_rectangles)){
      r <- list_of_rectangles_dec_with_decision_probs[[r_counter]]
      if(r$preset_decision==1){
        list_of_rectangles_dec_with_decision_probs[[counter]]$d_probs <- r$preset_decision_value
      } else{
        for(d in decisions){
          if(sum(d==r$allowed_decisions)>0){
            # use first rectangle in list_of_rectangles_mtp[[d]] as rprime_d in constraint
            rprime <- list_of_rectangles_mtp[[d]][[1]]
            variable_start_position <- variable_location(r,d,rprime,rprime$allowed_actions[1])
            variable_end_position <- variable_location(r,d,rprime,rprime$allowed_actions[length(rprime$allowed_actions)])
            list_of_rectangles_dec_with_decision_probs[[counter]]$d_probs[[d]] <-sum(z_solution[variable_start_position:variable_end_position])
          } else{list_of_rectangles_dec_with_decision_probs[[counter]]$d_probs[[d]] <- 0}
        }
      }
      counter<- counter+1
    }
    return(list_of_rectangles_dec_with_decision_probs )
  }

  list_of_rectangles_dec_with_decision_probs <- generate_list_of_rectangles_dec_with_decision_probs()
  counter <- 1

  ## Check if surrounded by rectangles with identical setting; if so, segment the rectangle; else if in center region then split it; else leave alone

  neighbor_check <- function(r,r2){
    if(
      (r$lower_boundaries[1]==r2$upper_boundaries[1] && max(c(r$lower_boundaries[2],r2$lower_boundaries[2]))<= min(c(r$upper_boundaries[2],r2$upper_boundaries[2])))
      || (r2$lower_boundaries[1]==r$upper_boundaries[1] && max(c(r$lower_boundaries[2],r2$lower_boundaries[2]))<= min(c(r$upper_boundaries[2],r2$upper_boundaries[2])))
      || (r$lower_boundaries[2]==r2$upper_boundaries[2] && max(c(r$lower_boundaries[1],r2$lower_boundaries[1]))<= min(c(r$upper_boundaries[1],r2$upper_boundaries[1])))
      || (r2$lower_boundaries[2]==r$upper_boundaries[2] && max(c(r$lower_boundaries[1],r2$lower_boundaries[1]))<= min(c(r$upper_boundaries[1],r2$upper_boundaries[1])))
    ) {return(1)} else {return(0)}
  }

  ## Code to check neighbor_check works
  #for(r in list_of_rectangles_dec_with_decision_probs){
  #plot(0,type="n",xlim=c(-10,10),ylim=c(-10,10),main=paste("Decision Rule for Stage 2 enrollment"),xlab="Z_stage_1_subpop_1",ylab="Z_stage_1_subpop_2",cex=2)
  #for(r2 in list_of_rectangles_dec_with_decision_probs){
  #	if(neighbor_check(r,r2) && r$lower_boundaries[1]>=(-3) && r$lower_boundaries[2]>=(-3) && r$upper_boundaries[1]<=3 && r$upper_boundaries[2]<=3){
  #		col_value <- 2
  #	rect(max(r$lower_boundaries[1],-10),max(r$lower_boundaries[2],-10),min(r$upper_boundaries[1],10),min(r$upper_boundaries[2],10),col=2)
  #	rect(max(r2$lower_boundaries[1],-10),max(r2$lower_boundaries[2],-10),min(r2$upper_boundaries[1],10),min(r2$upper_boundaries[2],10),col=3)
  #	browser()
  #	}
  #}}

  ## First set preliminary decision, without looking at surrounding rectangles
  for(r in list_of_rectangles_dec_with_decision_probs){
    for(d in r$allowed_decisions){
      if(!is.null(r$d_probs) &&  r$d_probs[d]>= 0.999){
        list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_preliminary <- 1
        list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_value_preliminary <- as.numeric(d==c(1,2,3,4))
      }

    }
    if (is.null(r$d_probs) || max(r$d_probs) < 0.999)
    {
      list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_preliminary <- 0
    }
    counter <- counter+1
  }

  ## Check if surrounded by identically set rectangles

  list_of_rectangles_dec_with_decision_probs_augmented <- list()
  for(counter in 1:(length(list_of_rectangles_dec_with_decision_probs)))
  {
    r <- list_of_rectangles_dec_with_decision_probs[[counter]]
    if(r$preset_decision_preliminary == 0 && r$lower_boundaries[1]>=(-4) && r$lower_boundaries[2]>=(-4) && r$upper_boundaries[1]<=4 && r$upper_boundaries[2]<=4){
      r_augmented <- r
      r_augmented$preset_decision_preliminary <- NULL
      r_augmented$preset_decision <- 0
      r_augmented$allowed_decisions <- decisions
      width <- r$upper_boundaries-r$lower_boundaries
      splitx <- splity <- 2
      for(counter2 in 1:splitx){
        for(counter3 in 1:splity){
          r_augmented$lower_boundaries <- (r$lower_boundaries + (c(width[1]*(counter2-1)/splitx,width[2]*(counter3-1)/splity)))
          r_augmented$upper_boundaries <- (r$lower_boundaries + (c(width[1]*(counter2)/splitx,width[2]*(counter3)/splity)))
          list_of_rectangles_dec_with_decision_probs_augmented <- c(list_of_rectangles_dec_with_decision_probs_augmented,list(r_augmented))
        }
      }
    } else if(r$preset_decision_preliminary == 0){
      r$preset_decision_preliminary <- NULL
      r$preset_decision <- 0
      list_of_rectangles_dec_with_decision_probs_augmented <- c(list_of_rectangles_dec_with_decision_probs_augmented,list(r))
    } else{ ## check if surrounded by same decision rectangles
      surrounded_by_same <- 1
      for(counter2 in 1:(length(list_of_rectangles_dec_with_decision_probs))){
        r2 <- list_of_rectangles_dec_with_decision_probs[[counter2]]
        if(neighbor_check(r,r2)==1 && !(r2$preset_decision_preliminary==1 && sum(r$preset_decision_value_preliminary == r2$preset_decision_value_preliminary)==4)) {		surrounded_by_same <- 0}
      }
      if(surrounded_by_same){
        r$preset_decision_preliminary <- NULL
        r$preset_decision <- 1
        r$preset_decision_value <- r$preset_decision_value_preliminary
        r$preset_decision_value_preliminary <- NULL
        list_of_rectangles_dec_with_decision_probs_augmented <- c(list_of_rectangles_dec_with_decision_probs_augmented,list(r))
      } else if(r$lower_boundaries[1]>=(-4) && r$lower_boundaries[2]>=(-4) && r$upper_boundaries[1]<=4 && r$upper_boundaries[2]<=4){
        r_augmented <- r
        r_augmented$preset_decision_preliminary <- NULL
        r_augmented$preset_decision <- 0
        r_augmented$allowed_decisions <- decisions
        width <- r$upper_boundaries-r$lower_boundaries
        splitx <- splity <- 2
        for(counter2 in 1:splitx){
          for(counter3 in 1:splity){
            r_augmented$lower_boundaries <- (r$lower_boundaries + (c(width[1]*(counter2-1)/splitx,width[2]*(counter3-1)/splity)))
            r_augmented$upper_boundaries <- (r$lower_boundaries + (c(width[1]*(counter2)/splitx,width[2]*(counter3)/splity)))
            list_of_rectangles_dec_with_decision_probs_augmented <- c(list_of_rectangles_dec_with_decision_probs_augmented,list(r_augmented))
          }
        }
      } else {
        r$preset_decision_preliminary <- NULL
        r$preset_decision <- 0
        list_of_rectangles_dec_with_decision_probs_augmented <- c(list_of_rectangles_dec_with_decision_probs_augmented,list(r))
      }
    }
  }

  list.of.rectangles.dec <- list_of_rectangles_dec_with_decision_probs_augmented

  postscript(paste("new_rectangles.eps"),height=8,horizontal=FALSE,onefile=FALSE,width=8)
  plot(0,type="n",xlim=c(-8,8),ylim=c(-8,8),main=paste("Decision Rule for Stage 2 enrollment"),xlab="Z_stage_1_subpop_1",ylab="Z_stage_1_subpop_2",cex=2)

  for(counter in 1:length(list_of_rectangles_dec_with_decision_probs_augmented))
  {
    r <- list_of_rectangles_dec_with_decision_probs_augmented[[counter]]
    if(r$preset_decision>0){color_value <- 2}else{color_value <- 3}
    rect(max(r$lower_boundaries[1],-10),max(r$lower_boundaries[2],-10),min(r$upper_boundaries[1],10),min(r$upper_boundaries[2],10),col=color_value)
  }
  dev.off()

  for(r1_counter in 1:(length(list_of_rectangles_dec_with_decision_probs_augmented)-1)){
    r1 <- list_of_rectangles_dec_with_decision_probs_augmented[[r1_counter]]
    for(r2_counter in (r1_counter+1):length(list_of_rectangles_dec_with_decision_probs_augmented)){
      r2 <- list_of_rectangles_dec_with_decision_probs_augmented[[r2_counter]]
      if((
        max(r1$lower_boundaries[1],r2$lower_boundaries[1]) < min(r1$upper_boundaries[1],r2$upper_boundaries[1])) && (max(r1$lower_boundaries[2],r2$lower_boundaries[2]) < min(r1$upper_boundaries[2],r2$upper_boundaries[2]))){print("ERROR: OVERLAPPING DECISION RECTANGLES"); print(c(r1,r2)); break;}
    }
  }


  set_counter <- 0
  for(r1_counter in 1:(length(list_of_rectangles_dec_with_decision_probs_augmented))){if(list_of_rectangles_dec_with_decision_probs_augmented[[r1_counter]]$preset_decision>0){set_counter<-set_counter+1}}
  print(set_counter)

  for( r in list.of.rectangles.dec){if(r$lower_boundaries[1]==r$upper_boundaries[1]) print(r)}
}
## In either case (rounding to integer or set and split):

for(d in decisions){list.of.rectangles.dec <- c(list.of.rectangles.dec,list(list(lower_boundaries=c(0,0),upper_boundaries=c(0,0),allowed_decisions=d,preset_decision=0)))}# these are reference rectangles

	   ## Set upper neighbors
	        counter_for_r <- 1
  		for(counter_for_r in 1:(length(list.of.rectangles.dec)-number_reference_rectangles)){
  		  list.of.rectangles.dec[[counter_for_r]]$upper_neighbor <- NULL
  		  r <- list.of.rectangles.dec[[counter_for_r]]
			count_value <- 1
		while(count_value <= length(list.of.rectangles.dec) && (!(r$upper_boundaries[2]== list.of.rectangles.dec[[count_value]]$lower_boundaries[2] && r$lower_boundaries[1]== list.of.rectangles.dec[[count_value]]$lower_boundaries[1] && r$upper_boundaries[1]== list.of.rectangles.dec[[count_value]]$upper_boundaries[1]))){count_value <- count_value +1}
		if(count_value <= length(list.of.rectangles.dec)){list.of.rectangles.dec[[counter_for_r]]$upper_neighbor <- count_value}
		}

    #save(list.of.rectangles.dec,file=paste("list.of.rectangles.dec",LP.iteration,".rdata",sep="")) # modified list of rectangles constructed

	set_counter <- 0
	for(r1_counter in 1:(length(list.of.rectangles.dec))){if(list.of.rectangles.dec[[r1_counter]]$preset_decision>0){set_counter<-set_counter+1}}
	number_preset_decision_rectangles <- set_counter
	print(number_preset_decision_rectangles)

return(list.of.rectangles.dec)
}
