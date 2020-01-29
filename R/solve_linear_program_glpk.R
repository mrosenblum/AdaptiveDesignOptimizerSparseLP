solve_linear_program_glpk <- function(total.alpha){
number_A1_files <- scan("number_A1_files.txt")
A1 = numeric(0)
for(i in 1:number_A1_files){
  tmp = load(paste("A1",i,".rdata",sep=""))
  tmp = constraint_list
  A1 = rbind(A1,tmp)
  print(i)
}

tmp = load("A2.rdata")
A2 = equality_constraints

tmp = load("A3.rdata")
A3  = rbind(power_constraint_matrix_H01, power_constraint_matrix_H02,power_constraint_matrix_H0C)

tmp    = load("power_constraints.rdata")
a3     = power.constraints

tmp = load("Inequality_Constraints_to_Restrict_MTP_to_Sufficient_Statistics.rdata")
A4 = additional_inequality_constraints_part1

tmp11 = read.table("number_equality_constraints_of_first_type.txt")
a21   = tmp11$V1

tmp = load("Inequality_Constraints_to_set_monotonicity_in_hypotheses_rejected.rdata")
A5  = additional_inequality_constraints_part2

tmp = load("c.rdata")
obj = objective_function_vector

A2  = Matrix::sparseMatrix(i = A2[,1],j=A2[,2],x=A2[,3])

A4  = Matrix::sparseMatrix(i = A4[,1],j=A4[,2],x=A4[,3])

tmp  = dim(A2)
tmp4 = dim(A4)

if (tmp4[2]<tmp[2])
{
  A4[tmp4[1],tmp4[2]] = 0
}

tmp = dim(A1)
a1  = total.alpha * rep(1,tmp[1])

tmp = dim(A2)
a2 = c(rep(1,a21),rep(0,tmp[1]-a21))

tmp = dim(A4)
a4  = rep(0,tmp[1])

AA    = rbind(A1,A4,-A3)
aa    = c(a1,a4,-power.constraints)

tmp   = length(obj)
lb    = rep(0,tmp)
ub    = rep(1,tmp)


tmpr = max(A5[,1])
tmpc = max(A5[,2])

tmp2  = dim(A1)


if (tmpc<tmp2[2])
{
  A5 = rbind(A5,c(tmpr,tmp2[2],0))
}
A5   = Matrix::sparseMatrix(i=A5[,1],j=A5[,2],x=A5[,3])
aan  = c(aa,rep(0,tmpr))
AAn  = rbind(AA,A5)

d         = dim(AAn)[2]

tmp       = cbind(seq(from=1,to=d,by=1),seq(from=1,to=d,by=1),rep(1,d))
At        = Matrix::sparseMatrix(i=tmp[,1],j=tmp[,2],x=tmp[,3])

A         = rbind(AAn,At,-At,A2)
a         = c(aan,rep(1,d),rep(0,d),a2)
const.dir = c(rep("<=",dim(AAn)[1]+2*d),rep("==",dim(A2)[1]))

tmp      = Rglpk::Rglpk_solve_LP(obj,A,const.dir,a,max=FALSE,control = list("tm_limit" = 120000))

sln        = list()
sln$z      = tmp$solution
sln$status = tmp$status
sln$dual   = tmp$solution_dual
sln$val    = tmp$optimum

return(sln)
}

