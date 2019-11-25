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
	(r$lower_boundaries[1]==r2$upper_boundaries[1] && max(c(r$lower_boundaries[2],r2$lower_boundaries[2]))< min(c(r$upper_boundaries[2],r2$upper_boundaries[2])))
 || (r2$lower_boundaries[1]==r$upper_boundaries[1] && max(c(r$lower_boundaries[2],r2$lower_boundaries[2]))< min(c(r$upper_boundaries[2],r2$upper_boundaries[2])))
 || (r$lower_boundaries[2]==r2$upper_boundaries[2] && max(c(r$lower_boundaries[1],r2$lower_boundaries[1]))< min(c(r$upper_boundaries[1],r2$upper_boundaries[1])))
 || (r2$lower_boundaries[2]==r$upper_boundaries[2] && max(c(r$lower_boundaries[1],r2$lower_boundaries[1]))< min(c(r$upper_boundaries[1],r2$upper_boundaries[1])))
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
		  if(!is.null(r$d_probs) &&  r$d_probs[d]>= 0.9){
				list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_preliminary <- 1
				list_of_rectangles_dec_with_decision_probs[[counter]]$preset_decision_value_preliminary <- as.numeric(d==c(1,2,3,4))
		  }

		  }
		  if (is.null(r$d_probs) || max(r$d_probs) < 0.9)
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

list_of_rectangles_dec_merged <- list_of_rectangles_dec_with_decision_probs_augmented

save(list_of_rectangles_dec_merged,file=paste("list_of_rectangles_dec_all_set.rdata",sep=""))


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

length(list_of_rectangles_dec_merged)

for( r in list_of_rectangles_dec_merged){if(r$lower_boundaries[1]==r$upper_boundaries[1]) print(r)}



