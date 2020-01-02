# Object returned
extract_solution <- function(list.of.rectangles.dec,decisions,list.of.rectangles.mtp,actions){
S1 <- list.of.rectangles.dec
A1 <- decisions
S2 <- list.of.rectangles.mtp
A2 <- actions
variable_location <- function(r,d,rprime,action){return(r$stage_1_rectangle_and_decision_offset[d]+rprime$stage_2_rectangle_offset+which(rprime$allowed_actions==action))}
pi_1 <- function(s1){# Returns vector of probabilities representing the policy pi^*_1 applied at state s1; the vector represents the probability that each possible enrollment decision is chosen given that the first stage z-statistics are in rectangle s1.
  probability_vector <- rep(0,length(decisions));
  for(d in decisions){
    s2 <- list.of.rectangles.mtp[[d]][[1]]
    variable_start_position <- variable_location(s1,d,s2,1);
    variable_end_position <- variable_location(s1,d,s2,length(actions));
    probability_vector[d] <- sum(sln$z[variable_start_position:variable_end_position]);}
  print(paste("Probabilities of enrollment decisions 1 through ",length(decisions), " respectively:"))
  return(probability_vector)
}

pi_2 <- function(s1,a1,s2){# Returns vector of probabilities representing the policy pi^*_2 applied after state s1, action a1, and state s2; the vector represents the probability that each possible subset of null hypotheses is rejected at the end of the trial
  variable_start_position <- variable_location(s1,a1,s2,1);
  variable_end_position <- variable_location(s1,a1,s2,length(actions));
  print(paste("Probabilities of rejecting each of the following subsets of null hypotheses: none, H01, H02, H0C, H01 and H0C, H02 and H0C, all"))
  return(sln$z[variable_start_position:variable_end_position])
}

return(policy = list(S1,A1,S2,A2,pi_1,pi_2));
}
