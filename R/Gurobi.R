library(gurobi)
library(Matrix)


number_A1_files <- scan("number_A1_files.txt")
A1 = numeric(0)
for(i in 1:number_A1_files){
        tmp = load(paste("A1",i,".rdata",sep=""))
        tmp = constraint_list

        A1 = rbind(A1,tmp)
}


tmp = load("A2.rdata")
A2  = equality_constraints
A2  = sparseMatrix(i=A2[,1],j=A2[,2],x=A2[,3])

tmp = load("A3.rdata")
A3  = rbind(power_constraint_matrix_H01, power_constraint_matrix_H02,power_constraint_matrix_H0C)

tmp    = load("power_constraints.rdata")
a3     = power.constraints

tmp = load("Inequality_Constraints_to_Restrict_MTP_to_Sufficient_Statistics.rdata")
A4  = additional_inequality_constraints_part1


col = dim(A2)
r4  = max(A4[,1])
c4  = max(A4[,2])


if (c4<col[2]){
     A4 = rbind(A4,c(r4,col[2],0))
}

A4  = sparseMatrix(i=A4[,1],j=A4[,2],x=A4[,3])

tmp11 = read.table("number_equality_constraints_of_first_type.txt")
a21   = tmp11$V1
tmp   = dim(A2)[1]
a2    = c(rep(1,a21),rep(0,tmp-a21))

tmp = load("Inequality_Constraints_to_set_monotonicity_in_hypotheses_rejected.rdata")
A5 = additional_inequality_constraints_part2
col = dim(A2)
r5  = max(A5[,1])
c5  = max(A5[,2])


if (c5<col[2]){
     A5 = rbind(A5,c(r5,col[2],0))
}

A5  = sparseMatrix(i=A5[,1],j=A5[,2],x=A5[,3])

tmpc = load("c.rdata")
cc   = objective_function_vector

AA  = rbind(A1,A4,-A3,A5)

a1  = rep(0.05,dim(A1)[1])

a4  = rep(0,dim(A4)[1])

aa  = c(a1,a4,-a3,rep(0,r5))

tmp = length(cc)

number_inequality_constraint = dim(AA)[1]
number_equality_constraint   = dim(A2)[1]

A = rbind(AA,A2)

b = c(aa,a2)

model <- list()

model$obj   <- cc
model$A     <- A
model$rhs   <- b
model$sense <- c(rep('<',number_inequality_constraint),rep('=',number_equality_constraint))
model$lb    <- rep(0,tmp)
model$ub    <- rep(1,tmp)

result      <- gurobi(model)

sln        = list()
sln$z      = result$x
sln$status = result$status
sln$dual   = result$pi
sln$val    = result$objval

i = 60
save(sln=sln,file=paste("sln2M",i,".rdata",sep="")
