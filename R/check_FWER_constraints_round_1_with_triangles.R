## Includes constraints the restrict multiple testing procedure not depend on stage 1 statistics (given cumulative statistics Z^C and decision d)
## Implements unequal rectangle sizes, with setting of rectangles to specific values
## Generate Discretized LP based on settings in file 
rm(list=ls())
job_division_number <- 1900
library(mvtnorm)
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))
# Create list of pairs of non-centrality parameters that the code below will evaluate the familywise Type I error rate (FWER) on.
## Encoding is: (ncp1_coordinate,ncp2_coordinate,width_value/2,indicator of triangle--special type); each element of this list represents a region in (ncp1,ncp2) space; if square (default) then first 2 components represent center of it and width_value represents side width; if triangle, then first 2 components represent upper right corner of right triangle and width/2 represents top length and right length of right triangle similar to triangle with vertices: (1,1),(1,0),(0,1).

ncp_fwer_check_list <- list()
for(x in seq(-15,15,length=301)){for(y in seq(-15,15,length=301)){
  ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(x,y,0.1,0))) ## 0.1 = (15 - -15)/(301-1) i.e., width of each square.
}}
#load("TOFILLIN_round_ncp_list_output_to_check_next_round.rdata")

if((1+(task_id-1)*ceiling(length(ncp_fwer_check_list)/job_division_number)) <= length(ncp_fwer_check_list)){
  ncp_fwer_check_list <- ncp_fwer_check_list[(1+(task_id-1)*ceiling(length(ncp_fwer_check_list)/job_division_number)):min(length(ncp_fwer_check_list),(task_id*ceiling(length(ncp_fwer_check_list)/job_division_number)))]
}else {ncp_fwer_check_list <- NULL}
print(ncp_fwer_check_list)

FWER_output <- array(0,c(length(ncp_fwer_check_list),7)) 
output_row_counter <- 1
source("discretized_optimization_problem_settings_adaptive_design_early_stopping_or_enrichment_two_subpopulations_normal_mixture_prior_can_set_dec_region.R")
load("z_rounded.rdata")
#Set mtp rectangle allowed decisions to be what is encoded in z_rounded
for(d_value in 1:4){
  r_reference_index <- length(list_of_rectangles_dec) - number_reference_rectangles + d_value
  r_reference <- list_of_rectangles_dec[[r_reference_index]]
  for(d in r_reference$allowed_decisions){
    for(rprime_counter in 1:length(list_of_rectangles_mtp[[d]])){
      rprime <- list_of_rectangles_mtp[[d]][[rprime_counter]]  
      variable_start_position <- variable_location(r_reference,d,rprime,rprime$allowed_actions[1])
      variable_end_position <- variable_location(r_reference,d,rprime,rprime$allowed_actions[length(rprime$allowed_actions)])
      action_indicator_vector <- z_rounded[variable_start_position:variable_end_position]
      list_of_rectangles_mtp[[d]][[rprime_counter]]$allowed_actions <- which(action_indicator_vector==1)
     #if(sum(rprime$lower_boundaries==c(-6.00,1.25))==2 && sum(rprime$upper_boundaries==c(-5.75,1.50))){browser()}
}}}

get_gradient_component <- function(upper_left_limit,lower_left_limit,upper_right_limit,lower_right_limit,xslope){
    integral_value <- 0
    ## Case 1: lower_left_limit[2] < upper_right[2]; triangle,rectangle,triangle (from top to bottom, integration regions)
  if(lower_left_limit[2] < upper_right_limit[2]){
    ## First integrate over upper triangle delineated by upper_left_limit, c(upper_left_limit[1],upper_right_limit[2]), upper_right_limit
    ## i.e., triangle with vertices: (x1,y1),(x1,y2),(x2,y1) defined as follows:
    x1<-upper_left_limit[1]
    x2<-upper_right_limit[1]
    y1<-upper_right_limit[2]
    y2<-upper_left_limit[2]
    avalue <- xslope #-(x2-x1)/(y2-y1) ##Slope for x=ay+b, upper line
    #if(!is.na(x2) && abs(xslope+(x2-x1)/(y2-y1))>0.00001){print("error");browser()}
    bvalue <- x2-y1*avalue     ##Intercept 
    
    upper_triangle_integral_value <- (1/sqrt(2*pi))*((pnorm(y2)-pnorm(y1))*exp(-x1^2/2) - exp(-(bvalue^2-avalue^2*bvalue^2/(avalue^2+1))/2)*(1/sqrt(avalue^2+1))*(pnorm((y2+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1)) -pnorm((y1+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1))))
    
    y1<-lower_left_limit[2]
    y2<-upper_right_limit[2]
    middle_rectangle_integral_value <- (1/sqrt(2*pi))*(pnorm(y2)-pnorm(y1))*(exp(-x1^2/2)-exp(-x2^2/2))
    
    if(lower_left_limit[2] != -Inf){
      y1<-lower_right_limit[2]
      y2<-lower_left_limit[2]
      avalue <- xslope #-(x2-x1)/(y2-y1) ##Slope for x=ay+b, lower line
      #if(abs(xslope+(x2-x1)/(y2-y1))>0.00001){print("error");browser()}
      bvalue <- x2-y1*avalue     ##Intercept 
      lower_triangle_value <- (1/sqrt(2*pi))*(-(pnorm(y2)-pnorm(y1))*exp(-x2^2/2) + exp(-(bvalue^2-avalue^2*bvalue^2/(avalue^2+1))/2)*(1/sqrt(avalue^2+1))*(pnorm((y2+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1)) -pnorm((y1+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1))))
    }else{lower_triangle_value <- 0}
    integral_value <- integral_value + upper_triangle_integral_value + middle_rectangle_integral_value + lower_triangle_value
    #if(is.nan(integral_value)){browser()}
  } else{  #  lower_left_limit[2] >= upper_right_limit[2]
    ## First integrate over upper triangle delineated by upper_left_limit,  lower_left_limit, upper line going through upper left limit
    ## i.e., triangle with vertices: (x1,y1),(x1,y2),(x2,y1) defined as follows:
    y1<-lower_left_limit[2]
    y2<-upper_left_limit[2]
    x1<-upper_left_limit[1]
    x2<-x1+(y2-y1)*(upper_right_limit[1]-x1)/(y2-upper_right_limit[2])    
    #if(!is.na(x2) && (x2 > upper_right_limit[1] || x2 < lower_left_limit[1]) #|| abs((-(sqrt(1-cor_value^2)/cor_value)*y1 + rprime$upper_boundaries[1]/cor_value) -x2)>1e-5
    #){print("error5");browser()}
    avalue <- xslope #-(x2-x1)/(y2-y1) ##Slope for x=ay+b, upper line
    #if(!is.na(x2) && abs(xslope+(x2-x1)/(y2-y1))>0.00001){print("error");browser()}
    bvalue <- upper_right_limit[1] - avalue*upper_right_limit[2] # x2-y1*avalue     ##Intercept 
    #if(!is.na(x2) && abs(bvalue-(x2-y1*avalue))>0.00001){print("error");browser()}
    
    upper_triangle_integral_value <- (1/sqrt(2*pi))*((pnorm(y2)-pnorm(y1))*exp(-x1^2/2) - exp(-(bvalue^2-avalue^2*bvalue^2/(avalue^2+1))/2)*(1/sqrt(avalue^2+1))*(pnorm((y2+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1)) -pnorm((y1+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1))))
    
    y1<-upper_right_limit[2]
    y2<-lower_left_limit[2]
    middle_parallelogram_integral_value <- (1/sqrt(2*pi))*(-exp(-(bvalue^2-avalue^2*bvalue^2/(avalue^2+1))/2)*(1/sqrt(avalue^2+1))*(pnorm((y2+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1)) -pnorm((y1+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1))))
    
    bvalue <- lower_right_limit[1]-lower_right_limit[2]*avalue     ##Intercept update lower line
    middle_parallelogram_integral_value <- middle_parallelogram_integral_value+(1/sqrt(2*pi))*(exp(-(bvalue^2-avalue^2*bvalue^2/(avalue^2+1))/2)*(1/sqrt(avalue^2+1))*(pnorm((y2+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1)) -pnorm((y1+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1))))  
    
    y1<-lower_right_limit[2]
    y2<-upper_right_limit[2]
    x2<-lower_right_limit[1]
    x1<-x2-(y2-y1)*(x2-lower_left_limit[1])/(lower_left_limit[2]-y1) 
    #if(!is.na(x1) && (x1 < lower_left_limit[1] || x1 > lower_right_limit[1]) #|| abs((-(sqrt(1-cor_value^2)/cor_value)*y2 + rprime$lower_boundaries[1]/cor_value)-x1)>1e-5
    #){print("error1");browser()}
    avalue <- xslope #-(x2-x1)/(y2-y1) ##Slope for x=ay+b, lower line
    #if(!is.na(x1) && !is.na(x2) && abs(xslope+(x2-x1)/(y2-y1))>0.00001){print("error");browser()}
    bvalue <- x2-y1*avalue     ##Intercept 
    lower_triangle_value <- (1/sqrt(2*pi))*(-(pnorm(y2)-pnorm(y1))*exp(-x2^2/2) + exp(-(bvalue^2-avalue^2*bvalue^2/(avalue^2+1))/2)*(1/sqrt(avalue^2+1))*(pnorm((y2+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1)) -pnorm((y1+avalue*bvalue/(avalue^2+1))*sqrt(avalue^2+1))))
    integral_value <- integral_value + upper_triangle_integral_value + middle_parallelogram_integral_value + lower_triangle_value
  }
  #if(is.nan(integral_value)){print("error6");browser()}
  #max_value <<- max(max_value,integral_value)
  return(integral_value)
}

print(date())
for(ncp_vector in ncp_fwer_check_list){
    ncp <- ncp_vector[1:2]
    gradient_vector <- array(0,c(length(decisions),4))
    Type_I_error_rate <- 0
    true_nulls_violated_by_action_set <- map_from_P_to_type_I_error_indicator_over_set_of_actions(ncp)
    for(r in list_of_rectangles_dec[1:(length(list_of_rectangles_dec)-number_reference_rectangles)]){
        for(d in r$allowed_decisions){
            if(d==1){
                for(rprime_counter in 1:length(list_of_rectangles_mtp[[d]])){
                    rprime <- list_of_rectangles_mtp[[d]][[rprime_counter]]
                    if(true_nulls_violated_by_action_set[rprime$allowed_actions]==1){# If causes a Type I error, then compute probability of this rectangle
                      prob_subpopulation_1 <- pmvnorm(mean=mean_vector(ncp,d)[c(1,3)],sigma=covariance_matrix[[d]][c(1,3),c(1,3)],lower=c(r$lower_boundaries[1],rprime$lower_boundaries[1]),upper=c(r$upper_boundaries[1],rprime$upper_boundaries[1]),algorithm=GenzBretz(abseps = 0.000000001,maxpts=100000))
                      prob_subpopulation_2 <- pmvnorm(mean=mean_vector(ncp,d)[c(2,4)],sigma=covariance_matrix[[d]][c(2,4),c(2,4)],lower=c(r$lower_boundaries[2],rprime$lower_boundaries[2]),upper=c(r$upper_boundaries[2],rprime$upper_boundaries[2]),algorithm=GenzBretz(abseps = 0.000000001,maxpts=100000))
                      combined_prob <- prob_subpopulation_1 * prob_subpopulation_2 #Uses independence of subpopulation-specific z-statistics under fixed decision d (but correlation occurs within subpopulation across stages)
                      Type_I_error_rate <- Type_I_error_rate + combined_prob
                      #Compute gradient of Type I error
                      gradient_component_z1 <- gradient_component_zF <- c(0,0)
                      for(subpopulation_number in 1:2){
                        if(subpopulation_number==1){
                          cor_value <- covariance_matrix[[d]][1,3] #correlation_z_1s_z_Fs for subpopulation s
                          subpopulation_mean_vector_values <- mean_vector(ncp,d)[c(1,3)]
                        } else {
                          cor_value <- covariance_matrix[[d]][2,4] #correlation_z_1s_z_Fs for subpopulation s
                          subpopulation_mean_vector_values <- mean_vector(ncp,d)[c(2,4)]
                        }
                        ## subtract mean vector from rectangle limits so that problem is centered:
                        z_1_left_limit_centered <- r$lower_boundaries[subpopulation_number] - subpopulation_mean_vector_values[1]
                        z_1_right_limit_centered <- r$upper_boundaries[subpopulation_number] - subpopulation_mean_vector_values[1]
                        z_F_left_limit_centered <- rprime$lower_boundaries[subpopulation_number] -  subpopulation_mean_vector_values[2]
                        z_F_right_limit_centered <- rprime$upper_boundaries[subpopulation_number] -  subpopulation_mean_vector_values[2]
                        ## First determine shape and limits of rectangle converted to z^1_1, z^2_1 plane
                        ## z^1_1 limits: r$lower_boundaries[1],r$upper_boundaries[1]
                        ## We use that cor_value z^1_1 + sqrt(1-cor_value^2) z^2_1 = z^F_1. 
                        ## z^2_1 constrained by lines: z^2_1 = -cor_value/sqrt(1-cor_value^2) z^1_1 + rprime$upper_boundaries[1]/sqrt(1-cor_value^2)
                        ##                        and  z^2_1 = -cor_value/sqrt(1-cor_value^2) z^1_1 + rprime$lower_boundaries[1]/sqrt(1-cor_value^2)
                        upper_left_limit <- c(z_1_left_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_left_limit_centered + z_F_right_limit_centered/sqrt(1-cor_value^2))
                        lower_left_limit <- c(z_1_left_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_left_limit_centered + z_F_left_limit_centered/sqrt(1-cor_value^2))
                        upper_right_limit <- c(z_1_right_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_right_limit_centered + z_F_right_limit_centered/sqrt(1-cor_value^2))
                        lower_right_limit <- c(z_1_right_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_right_limit_centered + z_F_left_limit_centered/sqrt(1-cor_value^2))
                        xslope <- (upper_left_limit[1]-upper_right_limit[1])/(upper_left_limit[2]-upper_right_limit[2]) #-(x2-x1)/(y2-y1)
                        if(upper_left_limit[2]< lower_left_limit[2] || upper_right_limit[2]< lower_right_limit[2] || upper_right_limit[2] > upper_left_limit[2] || lower_right_limit[2] > lower_left_limit[2] || xslope > 0){print("error1"); browser()}
                        gradient_component_z1[subpopulation_number] <- get_gradient_component(upper_left_limit,lower_left_limit,upper_right_limit,lower_right_limit,xslope)
                        ## z^F_1 component (final cumulative statistic, subpop 1):
                        ## transform to (z^F_1, z') plane where z'=sqrt(1-cor_value^2)z^1_1 - cor_value*z^2_1 is orthogonal to z^F_1 and normalized to variance 1
                        ## using tranformation z' = (1/sqrt(1-cor_value^2))*z^1_1 - (cor_value/sqrt(1-cor_value^2))*z^F_1
                        upper_left_limit <- c(z_F_left_limit_centered,(1/sqrt(1-cor_value^2))*z_1_right_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_left_limit_centered)
                        lower_left_limit <- c(z_F_left_limit_centered,(1/sqrt(1-cor_value^2))*z_1_left_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_left_limit_centered)
                        upper_right_limit <- c(z_F_right_limit_centered,(1/sqrt(1-cor_value^2))*z_1_right_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_right_limit_centered)
                        lower_right_limit <- c(z_F_right_limit_centered,(1/sqrt(1-cor_value^2))*z_1_left_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_right_limit_centered)
                        xslope <- -(sqrt(1-cor_value^2)/cor_value)#-(x2-x1)/(y2-y1)
                        if(upper_left_limit[2]< lower_left_limit[2] || upper_right_limit[2]< lower_right_limit[2] || upper_right_limit[2] > upper_left_limit[2] || lower_right_limit[2] > lower_left_limit[2]|| xslope>0){print("error1"); browser()}
                        gradient_component_zF[subpopulation_number] <- get_gradient_component(upper_left_limit,lower_left_limit,upper_right_limit,lower_right_limit,xslope)
                      }
                      gradient_vector[d,] <- gradient_vector[d,] + c(gradient_component_z1*c(prob_subpopulation_2,prob_subpopulation_1),gradient_component_zF*c(prob_subpopulation_2,prob_subpopulation_1))
                  }}} else if(d==2) {
                ## Get P_{ncp_vec}(Z \in rectangle r)	
                 for(rprime_counter in 1:length(list_of_rectangles_mtp[[d]])){
                    rprime <- list_of_rectangles_mtp[[d]][[rprime_counter]]
                    if(true_nulls_violated_by_action_set[rprime$allowed_actions]==1){# If causes a Type I error, then compute probability of this rectangle
                    ## compute overlap in rectangles in terms of z_2:
                    intersection_rectangle_lower_boundaries <- pmax(r$lower_boundaries,rprime$lower_boundaries)
                    intersection_rectangle_upper_boundaries <- pmin(r$upper_boundaries,rprime$upper_boundaries)
                    if(sum(intersection_rectangle_lower_boundaries < intersection_rectangle_upper_boundaries)==2){ #case where intersection rectangle non-empty
                      mean_vector_values <- mean_vector(ncp,d)
                      prob_subpopulation_1 <- pnorm(intersection_rectangle_upper_boundaries[1]-mean_vector_values[1])-pnorm(intersection_rectangle_lower_boundaries[1]-mean_vector_values[1])
                      prob_subpopulation_2 <- pnorm(intersection_rectangle_upper_boundaries[2]-mean_vector_values[2])-pnorm(intersection_rectangle_lower_boundaries[2]-mean_vector_values[2])
                      Type_I_error_rate <- Type_I_error_rate + prob_subpopulation_1*prob_subpopulation_2 #+ pmvnorm(mean=mean_vector(ncp,d),sigma=covariance_matrix[[d]],lower=intersection_rectangle_lower_boundaries,upper=intersection_rectangle_upper_boundaries,algorithm=GenzBretz(abseps = 0.000000001,maxpts=100000))
                    #Compute gradient of Type I error
                      gradient_component_z1 <- c(0,0)
                      for(subpopulation_number in 1:2){
                        subpopulation_mean_vector_values <- mean_vector(ncp,d)[subpopulation_number]
                        ## subtract mean vector from rectangle limits so that problem is centered:
                        z_1_left_limit_centered <- intersection_rectangle_lower_boundaries[subpopulation_number] - subpopulation_mean_vector_values
                        z_1_right_limit_centered <- intersection_rectangle_upper_boundaries[subpopulation_number] - subpopulation_mean_vector_values
                        gradient_component_z1[subpopulation_number] <- (1/sqrt(2*pi))*(exp(-z_1_left_limit_centered^2/2)-exp(-z_1_right_limit_centered^2/2))
                      }
                      gradient_vector[d,1:2] <- gradient_vector[d,1:2] + gradient_component_z1*c(prob_subpopulation_2,prob_subpopulation_1) 
              }}}} else if(d==3) {
                r_lower_boundary_z_2 <- r$lower_boundaries[2]
                r_upper_boundary_z_2 <- r$upper_boundaries[2]
                for(rprime_counter in 1:length(list_of_rectangles_mtp[[d]])){
                  rprime <- list_of_rectangles_mtp[[d]][[rprime_counter]]
                  if(true_nulls_violated_by_action_set[rprime$allowed_actions]==1){# If causes a Type I error, then compute probability of this rectangle
                    ## compute overlap in rectangles in terms of z_2:
                    effective_z_1_2_lower_boundary <- max(c(r_lower_boundary_z_2,rprime$lower_boundaries[2]))
                    effective_z_1_2_upper_boundary <- min(c(r_upper_boundary_z_2,rprime$upper_boundaries[2]))
                    if(effective_z_1_2_lower_boundary < effective_z_1_2_upper_boundary){ # case where rectangle non-empty
                        ## Get P_{ncp_vec}(Z \in rectangle r)	
                        prob_subpopulation_1 <- pmvnorm(mean=mean_vector(ncp,d)[c(1,3)],sigma=covariance_matrix[[d]][c(1,3),c(1,3)],lower=c(r$lower_boundaries[1],rprime$lower_boundaries[1]),upper=c(r$upper_boundaries[1],rprime$upper_boundaries[1]),algorithm=GenzBretz(abseps = 0.000000001,maxpts=100000))
                        prob_subpopulation_2 <- pnorm(effective_z_1_2_upper_boundary-mean_vector(ncp,d)[2])-pnorm(effective_z_1_2_lower_boundary-mean_vector(ncp,d)[2])
                        combined_prob <- prob_subpopulation_1 * prob_subpopulation_2 #Uses independence of subpopulation-specific z-statistics under fixed decision d (but correlation occurs within subpopulation across stages)
                        Type_I_error_rate <- Type_I_error_rate + combined_prob
                       
                        ## Compute Gradient of Type I Error component d==3 (enrich subpopulation 1)
                        gradient_component_z1 <- gradient_component_zF <- c(0,0)
                        ## First Handle Subpopulation 1 components
                        subpopulation_number <- 1
                        cor_value <- covariance_matrix[[d]][1,3] #correlation_z_1s_z_Fs for subpopulation s
                        subpopulation_mean_vector_values <- mean_vector(ncp,d)[c(1,3)]
                        ## subtract mean vector from rectangle limits so that problem is centered:
                        z_1_left_limit_centered <- r$lower_boundaries[subpopulation_number] - subpopulation_mean_vector_values[1]
                        z_1_right_limit_centered <- r$upper_boundaries[subpopulation_number] - subpopulation_mean_vector_values[1]
                        z_F_left_limit_centered <- rprime$lower_boundaries[subpopulation_number] -  subpopulation_mean_vector_values[2]
                        z_F_right_limit_centered <- rprime$upper_boundaries[subpopulation_number] -  subpopulation_mean_vector_values[2]
                        ## First determine shape and limits of rectangle converted to z^1_1, z^2_1 plane
                        ## z^1_1 limits: r$lower_boundaries[1],r$upper_boundaries[1]
                        ## We use that cor_value z^1_1 + sqrt(1-cor_value^2) z^2_1 = z^F_1. 
                        ## z^2_1 constrained by lines: z^2_1 = -cor_value/sqrt(1-cor_value^2) z^1_1 + rprime$upper_boundaries[1]/sqrt(1-cor_value^2)
                        ##                        and  z^2_1 = -cor_value/sqrt(1-cor_value^2) z^1_1 + rprime$lower_boundaries[1]/sqrt(1-cor_value^2)
                        upper_left_limit <- c(z_1_left_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_left_limit_centered + z_F_right_limit_centered/sqrt(1-cor_value^2))
                        lower_left_limit <- c(z_1_left_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_left_limit_centered + z_F_left_limit_centered/sqrt(1-cor_value^2))
                        upper_right_limit <- c(z_1_right_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_right_limit_centered + z_F_right_limit_centered/sqrt(1-cor_value^2))
                        lower_right_limit <- c(z_1_right_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_right_limit_centered + z_F_left_limit_centered/sqrt(1-cor_value^2))
                        xslope <- (upper_left_limit[1]-upper_right_limit[1])/(upper_left_limit[2]-upper_right_limit[2]) #-(x2-x1)/(y2-y1)
                        if(upper_left_limit[2]< lower_left_limit[2] || upper_right_limit[2]< lower_right_limit[2] || upper_right_limit[2] > upper_left_limit[2] || lower_right_limit[2] > lower_left_limit[2] || xslope > 0){print("error1"); browser()}
                        gradient_component_z1[subpopulation_number] <- get_gradient_component(upper_left_limit,lower_left_limit,upper_right_limit,lower_right_limit,xslope)
                        ## z^F_1 component (final cumulative statistic, subpop 1):
                        ## transform to (z^F_1, z') plane where z'=sqrt(1-cor_value^2)z^1_1 - cor_value*z^2_1 is orthogonal to z^F_1 and normalized to variance 1
                        ## using tranformation z' = (1/sqrt(1-cor_value^2))*z^1_1 - (cor_value/sqrt(1-cor_value^2))*z^F_1
                        upper_left_limit <- c(z_F_left_limit_centered,(1/sqrt(1-cor_value^2))*z_1_right_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_left_limit_centered)
                        lower_left_limit <- c(z_F_left_limit_centered,(1/sqrt(1-cor_value^2))*z_1_left_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_left_limit_centered)
                        upper_right_limit <- c(z_F_right_limit_centered,(1/sqrt(1-cor_value^2))*z_1_right_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_right_limit_centered)
                        lower_right_limit <- c(z_F_right_limit_centered,(1/sqrt(1-cor_value^2))*z_1_left_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_right_limit_centered)
                        xslope <- -(sqrt(1-cor_value^2)/cor_value)#-(x2-x1)/(y2-y1)
                        if(upper_left_limit[2]< lower_left_limit[2] || upper_right_limit[2]< lower_right_limit[2] || upper_right_limit[2] > upper_left_limit[2] || lower_right_limit[2] > lower_left_limit[2]|| xslope>0){print("error1"); browser()}
                        gradient_component_zF[subpopulation_number] <- get_gradient_component(upper_left_limit,lower_left_limit,upper_right_limit,lower_right_limit,xslope)
                        ##
                        ## Next Handle Subpopulation 2 components
                        subpopulation_number <- 2
                        subpopulation_mean_vector_values <- mean_vector(ncp,d)[subpopulation_number]
                        ## subtract mean vector from rectangle limits so that problem is centered:
                        z_1_left_limit_centered <- effective_z_1_2_lower_boundary - subpopulation_mean_vector_values
                        z_1_right_limit_centered <- effective_z_1_2_upper_boundary - subpopulation_mean_vector_values
                        gradient_component_z1[subpopulation_number] <- (1/sqrt(2*pi))*(exp(-z_1_left_limit_centered^2/2)-exp(-z_1_right_limit_centered^2/2))
                        ## Lastly, combine with subpopulation probabilities
                        gradient_vector[d,] <- gradient_vector[d,] + c(gradient_component_z1,gradient_component_zF)*c(prob_subpopulation_2,prob_subpopulation_1,prob_subpopulation_2,0)
                    }}}} else if(d==4) { #only subpopulation 2 enrolled further
                r_lower_boundary_z_1 <- r$lower_boundaries[1]
                r_upper_boundary_z_1 <- r$upper_boundaries[1]
                for(rprime_counter in 1:length(list_of_rectangles_mtp[[d]])){
                  rprime <- list_of_rectangles_mtp[[d]][[rprime_counter]]
                  if(true_nulls_violated_by_action_set[rprime$allowed_actions]==1){# If causes a Type I error, then compute probability of this rectangle
                    ## compute overlap in rectangles in terms of z_2:
                    effective_z_1_1_lower_boundary <- max(c(r_lower_boundary_z_1,rprime$lower_boundaries[1]))
                    effective_z_1_1_upper_boundary <- min(c(r_upper_boundary_z_1,rprime$upper_boundaries[1]))
                    if(effective_z_1_1_lower_boundary < effective_z_1_1_upper_boundary){ # case where rectangle non-empty
                        ## Get P_{ncp_vec}(Z \in rectangle r)	
                      ## Get P_{ncp_vec}(Z \in rectangle r)	
                       prob_subpopulation_1 <- pnorm(effective_z_1_1_upper_boundary-mean_vector(ncp,d)[1])-pnorm(effective_z_1_1_lower_boundary-mean_vector(ncp,d)[1])
                       prob_subpopulation_2 <- pmvnorm(mean=mean_vector(ncp,d)[c(2,3)],sigma=covariance_matrix[[d]][c(2,3),c(2,3)],lower=c(r$lower_boundaries[2],rprime$lower_boundaries[2]),upper=c(r$upper_boundaries[2],rprime$upper_boundaries[2]),algorithm=GenzBretz(abseps = 0.000000001,maxpts=100000))
                       combined_prob <- prob_subpopulation_1 * prob_subpopulation_2 #Uses independence of subpopulation-specific z-statistics under fixed decision d (but correlation occurs within subpopulation across stages)
                      Type_I_error_rate <- Type_I_error_rate + combined_prob
                      ## Compute Gradient of Type I Error component d==3 (enrich subpopulation 1)
                      gradient_component_z1 <- gradient_component_zF <- c(0,0)
                      ## First Handle Subpopulation 1 components
                      ##
                      subpopulation_number <- 2
                      cor_value <- covariance_matrix[[d]][2,3] #correlation_z_1s_z_Fs for subpopulation s
                      subpopulation_mean_vector_values <- mean_vector(ncp,d)[c(2,3)]
                      ## subtract mean vector from rectangle limits so that problem is centered:
                      z_1_left_limit_centered <- r$lower_boundaries[subpopulation_number] - subpopulation_mean_vector_values[1]
                      z_1_right_limit_centered <- r$upper_boundaries[subpopulation_number] - subpopulation_mean_vector_values[1]
                      z_F_left_limit_centered <- rprime$lower_boundaries[subpopulation_number] -  subpopulation_mean_vector_values[2]
                      z_F_right_limit_centered <- rprime$upper_boundaries[subpopulation_number] -  subpopulation_mean_vector_values[2]
                      ## First determine shape and limits of rectangle converted to z^1_1, z^2_1 plane
                      ## z^1_1 limits: r$lower_boundaries[1],r$upper_boundaries[1]
                      ## We use that cor_value z^1_1 + sqrt(1-cor_value^2) z^2_1 = z^F_1. 
                      ## z^2_1 constrained by lines: z^2_1 = -cor_value/sqrt(1-cor_value^2) z^1_1 + rprime$upper_boundaries[1]/sqrt(1-cor_value^2)
                      ##                        and  z^2_1 = -cor_value/sqrt(1-cor_value^2) z^1_1 + rprime$lower_boundaries[1]/sqrt(1-cor_value^2)
                      upper_left_limit <- c(z_1_left_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_left_limit_centered + z_F_right_limit_centered/sqrt(1-cor_value^2))
                      lower_left_limit <- c(z_1_left_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_left_limit_centered + z_F_left_limit_centered/sqrt(1-cor_value^2))
                      upper_right_limit <- c(z_1_right_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_right_limit_centered + z_F_right_limit_centered/sqrt(1-cor_value^2))
                      lower_right_limit <- c(z_1_right_limit_centered,-(cor_value/sqrt(1-cor_value^2))*z_1_right_limit_centered + z_F_left_limit_centered/sqrt(1-cor_value^2))
                      xslope <- (upper_left_limit[1]-upper_right_limit[1])/(upper_left_limit[2]-upper_right_limit[2]) #-(x2-x1)/(y2-y1)
                      if(upper_left_limit[2]< lower_left_limit[2] || upper_right_limit[2]< lower_right_limit[2] || upper_right_limit[2] > upper_left_limit[2] || lower_right_limit[2] > lower_left_limit[2] || xslope > 0){print("error1"); browser()}
                      gradient_component_z1[subpopulation_number] <- get_gradient_component(upper_left_limit,lower_left_limit,upper_right_limit,lower_right_limit,xslope)
                      ## z^F_1 component (final cumulative statistic, subpop 1):
                      ## transform to (z^F_1, z') plane where z'=sqrt(1-cor_value^2)z^1_1 - cor_value*z^2_1 is orthogonal to z^F_1 and normalized to variance 1
                      ## using tranformation z' = (1/sqrt(1-cor_value^2))*z^1_1 - (cor_value/sqrt(1-cor_value^2))*z^F_1
                      upper_left_limit <- c(z_F_left_limit_centered,(1/sqrt(1-cor_value^2))*z_1_right_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_left_limit_centered)
                      lower_left_limit <- c(z_F_left_limit_centered,(1/sqrt(1-cor_value^2))*z_1_left_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_left_limit_centered)
                      upper_right_limit <- c(z_F_right_limit_centered,(1/sqrt(1-cor_value^2))*z_1_right_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_right_limit_centered)
                      lower_right_limit <- c(z_F_right_limit_centered,(1/sqrt(1-cor_value^2))*z_1_left_limit_centered - (cor_value/sqrt(1-cor_value^2))*z_F_right_limit_centered)
                      xslope <- -(sqrt(1-cor_value^2)/cor_value)#-(x2-x1)/(y2-y1)
                      if(upper_left_limit[2]< lower_left_limit[2] || upper_right_limit[2]< lower_right_limit[2] || upper_right_limit[2] > upper_left_limit[2] || lower_right_limit[2] > lower_left_limit[2]|| xslope>0){print("error1"); browser()}
                      gradient_component_zF[subpopulation_number] <- get_gradient_component(upper_left_limit,lower_left_limit,upper_right_limit,lower_right_limit,xslope)
                      ##
                      ## Next Handle Subpopulation 2 components
                      subpopulation_number <- 1
                      subpopulation_mean_vector_values <- mean_vector(ncp,d)[subpopulation_number]
                      ## subtract mean vector from rectangle limits so that problem is centered:
                      z_1_left_limit_centered <- effective_z_1_1_lower_boundary - subpopulation_mean_vector_values
                      z_1_right_limit_centered <- effective_z_1_1_upper_boundary - subpopulation_mean_vector_values
                      gradient_component_z1[subpopulation_number] <- (1/sqrt(2*pi))*(exp(-z_1_left_limit_centered^2/2)-exp(-z_1_right_limit_centered^2/2))
                      ## Lastly, combine with subpopulation probabilities
                      gradient_vector[d,] <- gradient_vector[d,] + c(gradient_component_z1,gradient_component_zF)*c(prob_subpopulation_2,prob_subpopulation_1,0,prob_subpopulation_1)
                    }}}}
        }}
    # Construct gradient using formula S-6 in Web Appendix based on components computed above
    # For subpopulation 1
    # Define square root inverse of covariance matrices
    sigma_inverse <- list()
    for(d in decisions){sigma_inverse[[d]] <- solve(covariance_matrix[[d]])}
    gradient_component1 <-  ((mean_vector(c(1,0),1)%*%sigma_inverse[[1]]%*%gradient_vector[1,])
                            +(mean_vector(c(1,0),2)%*%sigma_inverse[[2]]%*%gradient_vector[2,1:2])
                            +(mean_vector(c(1,0),3)%*%sigma_inverse[[3]]%*%gradient_vector[3,1:3])
                            +(mean_vector(c(1,0),4)%*%sigma_inverse[[4]]%*%gradient_vector[4,c(1,2,4)]))
    # For subpopulation 2
    gradient_component2 <-  ((mean_vector(c(0,1),1)%*%sigma_inverse[[1]]%*%gradient_vector[1,])
                            +(mean_vector(c(0,1),2)%*%sigma_inverse[[2]]%*%gradient_vector[2,1:2])
                            +(mean_vector(c(0,1),3)%*%sigma_inverse[[3]]%*%gradient_vector[3,1:3])
                            +(mean_vector(c(0,1),4)%*%sigma_inverse[[4]]%*%gradient_vector[4,c(1,2,4)]))
    print(c(ncp,Type_I_error_rate,gradient_component1,gradient_component2))
    # Store results
    FWER_output[output_row_counter,] <- c(ncp,Type_I_error_rate,gradient_component1,gradient_component2,ncp_vector[3],ncp_vector[4])
    output_row_counter <- output_row_counter + 1
}
print(date())
save(FWER_output,file=paste("./First_Round_Computations/FWER_computation_with_gradient_output",task_id,".rdata",sep=""))
