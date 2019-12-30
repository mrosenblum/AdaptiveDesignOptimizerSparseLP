#install.packages("cubature",repos='http://cran.us.r-project.org')
rm(list=ls())
ncp_fwer_check_list <- list()
#plot(0,xlim=c(-6,6),ylim=c(-6,6))
## FWER_output encoding is: (ncp1,ncp2,FWER,gradient1,grandient2,width_value/2,triangle_indicator)
for(previous_task_id in 1:1900){
	load(paste("First_Round_Computations/FWER_computation_with_gradient_output",previous_task_id,".rdata",sep=""))
        if(dim(FWER_output)[1]==0)print(previous_task_id) else{
	for(i in 1:dim(FWER_output)[1]){
	  	    w <- FWER_output[i,6] # extract square side-length, i.e., width_value
	    if((FWER_output[i,7]==0) && # case where region is square centered at (FWER_output[i,1], FWER_output[i,2])+/- width_value/2 in each dimension
	  	!(FWER_output[i,3]+(w/2)*(abs(FWER_output[i,4])+(w/2)*4)+ (w/2)*(abs(FWER_output[i,5])+(w/2)*4)<=0.05)){# case where FWER controlled within box defined by point +/- w/2 in each dimension, so not necessary to consider any further
	  	  
		# Consider upper right quadrant and if FWER bounded above by 0.05 then can omit it
		if(!(FWER_output[i,3] <=0.05 && (FWER_output[i,3]+(w/2)*(FWER_output[i,4]+FWER_output[i,5])+8*(w/2)^2<=0.05) &&
		          (FWER_output[i,3]+(w/2)*FWER_output[i,4]+4*(w/2)^2<=0.05) && (FWER_output[i,3]+(w/2)*FWER_output[i,5]+4*(w/2)^2<=0.05))){
		  ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(FWER_output[i,1]+(w/4),FWER_output[i,2]+(w/4),w/2,0))) # add upper right quadrant as new region to evaluate
		}
		  
		# Consider lower right quadrant
		if(!(FWER_output[i,3] <=0.05 && (FWER_output[i,3]+(w/2)*(FWER_output[i,4]-FWER_output[i,5])+8*(w/2)^2<=0.05) &&
		     (FWER_output[i,3]+(w/2)*FWER_output[i,4]+4*(w/2)^2<=0.05) && (FWER_output[i,3]+(w/2)*(-FWER_output[i,5])+4*(w/2)^2<=0.05))){ #check lower right quadrant and if FWER bounded above by 0.05 then can omit it
		  if((abs(FWER_output[i,1]+FWER_output[i,2])<1e-10) && (FWER_output[i,3] <= 0.05) && (FWER_output[i,3]+(w/2)*(FWER_output[i,4]-FWER_output[i,5])+8*(w/2)^2<=0.05) && (FWER_output[i,3]+(w/2)*(-FWER_output[i,5])+4*(w/2)^2<=0.05)){ #if on H0C null space boundary and can rule out lower triangle then add upper triangle to list
		  ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(FWER_output[i,1]+(w/2),FWER_output[i,2],w/2,1)))
		  } else {ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(FWER_output[i,1]+(w/4),FWER_output[i,2]-(w/4),w/2,0))) # add lower right quadrant as new region to evaluate}
		  }
		}
	      
		  # Consider lower left quadrant and if FWER bounded above by 0.05 then can omit it 
		  if(!(FWER_output[i,3] <=0.05 && (FWER_output[i,3]+(w/2)*(-FWER_output[i,4]-FWER_output[i,5])+8*(w/2)^2<=0.05) &&
		     (FWER_output[i,3]+(w/2)*(-FWER_output[i,4])+4*(w/2)^2<=0.05) && (FWER_output[i,3]+(w/2)*(-FWER_output[i,5])+4*(w/2)^2<=0.05))){
		    ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(FWER_output[i,1]-(w/4),FWER_output[i,2]-(w/4),w/2,0))) # add lower left quadrant as new region to evaluate
		  }
		  
		  # Consider upper left quadrant
		  if(!(FWER_output[i,3] <=0.05 && (FWER_output[i,3]+(w/2)*(-FWER_output[i,4]+FWER_output[i,5])+8*(w/2)^2<=0.05) &&
		     (FWER_output[i,3]+(w/2)*(-FWER_output[i,4])+4*(w/2)^2<=0.05) && (FWER_output[i,3]+(w/2)*FWER_output[i,5]+4*(w/2)^2<=0.05))){ #check upper left quadrant and if FWER bounded above by 0.05 then can omit it
		     if((abs(FWER_output[i,1]+FWER_output[i,2])<1e-10) && (FWER_output[i,3] <= 0.05) && (FWER_output[i,3]+(w/2)*(-FWER_output[i,4]+FWER_output[i,5])+8*(w/2)^2<=0.05) && (FWER_output[i,3]+(w/2)*(-FWER_output[i,4])+4*(w/2)^2<=0.05)){ #if on H0C null space boundary and can rule out lower triangle then add upper triangle to list
		      ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(FWER_output[i,1],FWER_output[i,2]+(w/2),w/2,1)))
		      } else {ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(FWER_output[i,1]-(w/4),FWER_output[i,2]+(w/4),w/2,0))) # add upper left quadrant as new region to evaluate}
		    }
		  }
	  } else if(FWER_output[i,7]==1){ # triangle
		  if(!(FWER_output[i,3] <=0.05 && (FWER_output[i,3]+w*(-FWER_output[i,4])+4*w^2<=0.05) && (FWER_output[i,3]+w*(-FWER_output[i,5])+4*w^2<=0.05))){# if can rule out triangular region then omit it, else cut it into 2 smaller triangles and square
		    ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(FWER_output[i,1]-(w/2),FWER_output[i,2],w/2,1)))
		    ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(FWER_output[i,1],FWER_output[i,2]+(w/2),w/2,1)))
		    ncp_fwer_check_list <- c(ncp_fwer_check_list,list(c(FWER_output[i,1]-(w/4),FWER_output[i,2]-(w/4),w/2,0)))
		  }
	  }
}}}
	 
print(length(ncp_fwer_check_list))
print(length(ncp_fwer_check_list)/190)
new_w <- w/2
save(new_w,ncp_fwer_check_list,file="first_round_ncp_list_output_to_check_next_round.rdata")
