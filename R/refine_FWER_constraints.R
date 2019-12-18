refine_FWER_constraints <- function(){
   # list of pairs of non-centrality parameters in G_{tau,w}
   ##need to input following information:
   #load("sln2M1.rdata") #TOFIX
   #load("ncp.list1.rdata") #TOFIX
   #load("list_of_rectangles_dec1.rdata") # Use same discretization of decision region as in iteration 1
   number_preset_decision_rectangles <- 0
   z_solution <- sln$z
   # Get active FWER constraints
   ncp_active_FWER_constraints <- ncp.list[which(sln$dual[1:length(ncp.list)]>0.00001)]

   #Create smaller set of FWER constraints for use at future iterations
   ncp.list <- list()
   for(z in seq(-9,9,length=18)) {ncp.list <- c(ncp.list,list(c(0,z)))};
   for(z in seq(-9,9,length=18)) {ncp.list <- c(ncp.list,list(c(z,0)))};
   for(z in seq(-9,9,length=18)) {ncp.list <- c(ncp.list,list(c(z,-(rho1/rho2)*z)))};

   for(ncp in ncp_active_FWER_constraints){
	if(abs(ncp[1]) < 1e-10 && abs(ncp[2])>1e-10){ # constraint on boundary of H01 null space
		       for(z in seq(ncp[2]-0.5,ncp[2]+0.5,length=5)) {ncp.list <- c(ncp.list,list(c(0,z)))}} else
	if(abs(ncp[1]) > 1e-10 && abs(ncp[2])<1e-10){ # constraint on boundary of H02 null space
		       for(z in seq(ncp[1]-0.5,ncp[1]+0.5,length=5)) {ncp.list <- c(ncp.list,list(c(z,0)))}} else
	if(abs(ncp[1]) > 1e-10 && abs(ncp[2])>1e-10){ # constraint on boundary of H0C null space
		       for(z in seq(ncp[1]-0.5,ncp[1]+0.5,length=5)) {ncp.list <- c(ncp.list,list(c(z,-(rho1/rho2)*z)))}}
   }
   ncp.list <- c(ncp.list,list(c(0,0)))
   ncp.list <- unique(ncp.list)
   save(ncp.list,file=paste("ncp.list_short2.rdata",sep=""))

   # Create new set of FWER constraints concentrated around active constraints from previous linear program solution
   ncp.list <- list()
   for(z in seq(-9,9,length=18)) {ncp.list <- c(ncp.list,list(c(0,z)))};
   for(z in seq(-9,9,length=18)) {ncp.list <- c(ncp.list,list(c(z,0)))};
   for(z in seq(-9,9,length=18)) {ncp.list <- c(ncp.list,list(c(z,-(rho1/rho2)*z)))};

   number_new_constraints_per_previous_active_constraint <- max(1,floor(500/length(ncp_active_FWER_constraints)))
   for(ncp in ncp_active_FWER_constraints){
	if(abs(ncp[1]) < 1e-10 && abs(ncp[2])>1e-10){ # constraint on boundary of H01 null space
		       for(z in seq(ncp[2]-0.5,ncp[2]+0.5,length=number_new_constraints_per_previous_active_constraint)) {ncp.list <- c(ncp.list,list(c(0,z)))}} else
	if(abs(ncp[1]) > 1e-10 && abs(ncp[2])<1e-10){ # constraint on boundary of H02 null space
		       for(z in seq(ncp[1]-0.5,ncp[1]+0.5,length=number_new_constraints_per_previous_active_constraint)) {ncp.list <- c(ncp.list,list(c(z,0)))}} else
	if(abs(ncp[1]) > 1e-10 && abs(ncp[2])>1e-10){ # constraint on boundary of H0C null space
		       for(z in seq(ncp[1]-0.5,ncp[1]+0.5,length=number_new_constraints_per_previous_active_constraint)) {ncp.list <- c(ncp.list,list(c(z,-(rho1/rho2)*z)))}}
   }
   ncp.list <- c(ncp.list,list(c(0,0)))
   ncp.list <- unique(ncp.list)
   constraints_per_A1_file <- 1

   save(list_of_rectangles_dec,file=paste("list_of_rectangles_dec",LP.iteration,".rdata",sep="")) # list of rectangles for decision region same as LP.iteration 1
   save(ncp.list,file=paste("ncp.list",LP.iteration,".rdata",sep=""))
} 

