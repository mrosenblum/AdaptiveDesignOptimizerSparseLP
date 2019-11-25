## Visualize solutions
source("discretized_optimization_problem_settings_adaptive_design_early_stopping_or_enrichment_two_subpopulations_normal_mixture_prior_can_set_dec_region.R")
#load("c.rdata")

load(paste("sln2M83.RData",sep=""))
z_solution <- sln2$z 
#Explore Dual
ncp_list[which(sln2$dual[1:length(ncp_list)]>0.01)]

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

save(list_of_rectangles_dec_with_decision_probs_merged,file=paste("list_of_rectangles_dec_with_decision_probs_merged.rdata",sep=""))

## To view results of rounding and merging

postscript(paste("new_rectangles.eps"),height=8,horizontal=FALSE,onefile=FALSE,width=8)
plot(0,type="n",xlim=c(-8,8),ylim=c(-8,8),main=paste("Decision Rule for Stage 2 enrollment"),xlab="Z_stage_1_subpop_1",ylab="Z_stage_1_subpop_2",cex=2)

for(counter in 1:length(list_of_rectangles_dec_with_decision_probs_merged))
{
  r <- list_of_rectangles_dec_with_decision_probs_merged[[counter]]
  if(r$preset_decision>0){color_value <- r$preset_decision}else{color_value <- 5}
  rect(max(r$lower_boundaries[1],-10),max(r$lower_boundaries[2],-10),min(r$upper_boundaries[1],10),min(r$upper_boundaries[2],10),col=color_value)	
}
dev.off()
