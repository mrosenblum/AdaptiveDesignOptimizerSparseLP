#' This function takes the output of the linear program solver and constructs the corresponding policy.
extract_solution <- function(list.of.rectangles.dec,decisions,list.of.rectangles.mtp,actions,sln){
S1 <- list.of.rectangles.dec[1:(length(list.of.rectangles.dec)-length(decisions))];
A1 <- decisions;
S2 <- list.of.rectangles.mtp;
A2 <- actions;
names(A2) <- c("Reject none","Reject H01","Reject H02","Reject H0C","Reject H01 and H0C","Reject H02 and H0C","Reject all");
variable_location <- function(r,d,rprime,action){return(r$stage_1_rectangle_and_decision_offset[d]+rprime$stage_2_rectangle_offset+which(rprime$allowed_actions==action))}
pi_1 <- function(s1_index,print.explanation=TRUE){# Returns vector of probabilities representing the policy pi^*_1 applied at state s1; the vector represents the probability that each possible enrollment decision is chosen given that the first stage z-statistics are in rectangle s1.
  if(s1_index < 1 || s1_index > length(S1)){print("s1_index is out of the range of the possible S1 states."); return(NULL)};
  s1 <- S1[[s1_index]];
  probability_vector <- rep(0,length(decisions));
  if(s1$preset_decision==0){
    for(d in decisions){
      s2 <- list.of.rectangles.mtp[[d]][[1]]
      variable_start_position <- variable_location(s1,d,s2,1);
      variable_end_position <- variable_location(s1,d,s2,length(actions));
      probability_vector[d] <- sum(sln$z[variable_start_position:variable_end_position]);}
  } else {
    probability_vector <- s1$preset_decision_value;
  }
  if(print.explanation==TRUE){
  print(paste("Probabilities of enrollment decisions 1 through ",length(decisions), " respectively:"))}
  return(probability_vector)
}
pi_2 <- function(s1_index,a1_index,s2_index){# Returns vector of probabilities representing the policy pi^*_2 applied after state s1, action a1, and state s2; the vector represents the probability that each possible subset of null hypotheses is rejected at the end of the trial
  if(s1_index < 1 || s1_index > length(S1)){print("s1_index is out of the range of the possible S1 states."); return(NULL)} else if(a1_index < 0 || a1_index > length(A1)){print("a1_index is out of the range of the possible A1 actions."); return(NULL);}else if(s2_index < 0 || s2_index > length(S2[[a1_index]])){print("a1_index is out of the range of the possible A1 actions."); return(NULL);}
    s1 <- S1[[s1_index]];
    a1 <- A1[[a1_index]];
    s2 <- S2[[a1]][[s2_index]];
    print("Inputs correspond to the following:")
    print(paste("s1 is rectangle defined as the Cartesian product [",s1$lower_boundaries[1],",",s1$upper_boundaries[1],"] x [",s1$lower_boundaries[2],",",s1$upper_boundaries[2],"]"));
    print(paste("a1 is enrollment decision ",a1))
    print(paste("s2 is rectangle defined as the Cartesian product [",s2$lower_boundaries[1],",",s2$upper_boundaries[1],"] x [",s2$lower_boundaries[2],",",s2$upper_boundaries[2],"]"));
    probability_of_a1_given_s1 <- pi_1(s1_index,print.explanation=FALSE)[a1];
    if(probability_of_a1_given_s1<=0){print("That sequence of states and actions is not possible under this policy since a1 has probability 0 given s1.");return(NULL)}
    if(s1$preset_decision==0){
      variable_start_position <- variable_location(s1,a1,s2,1);
      variable_end_position <- variable_location(s1,a1,s2,length(actions));
      probability_vector <- sln$z[variable_start_position:variable_end_position]/probability_of_a1_given_s1;
      probability_vector <- pmax(probability_vector,0)
      names(probability_vector) <- c("Reject none","Reject H01","Reject H02","Reject H0C","Reject H01 and H0C","Reject H02 and H0C","Reject all");
      return(probability_vector);} else if(s1$preset_decision_value[a1]==1){
    #Deterministic decision to take action a1 if first stage statistics are in s1
    reference_rectangle_s1 <- list.of.rectangles.dec[[length(list.of.rectangles.dec)-length(decisions)+which(s1$preset_decision_value==1)]];
    variable_start_position <- variable_location(reference_rectangle_s1,a1,s2,1);
    variable_end_position <- variable_location(reference_rectangle_s1,a1,s2,length(actions));
    print(paste("Probabilities of rejecting each of the following subsets of null hypotheses, respectively:"))
    probability_vector <- sln$z[variable_start_position:variable_end_position]/probability_of_a1_given_s1;
    probability_vector <- pmax(probability_vector,0);
    names(probability_vector) <- c("Reject none","Reject H01","Reject H02","Reject H0C","Reject H01 and H0C","Reject H02 and H0C","Reject all");
    return(probability_vector);} else {print(paste("That sequence of states and actions is not possible under this policy, since after state s1 the enrollment decision ",which(s1$preset_decision_value==1)," is always selected."));return(NULL)}
}
return(policy = list(S1=S1,A1=A1,S2=S2,A2=A2,pi_1=pi_1,pi_2=pi_2));
}
