library(Matrix)
library(Rglpk)

#setwd("/Users/ethanfangxy/Dropbox/Michael/RsymphonyTest2")

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
A3  = rbind(power_constraint_vector_subpopulation1,power_constraint_vector_subpopulation2,power_constraint_vector_combined_population)
#A3  = rbind(power_constraint_vector_treatment1,power_constraint_vector_treatment2,power_constraint_vector_both_treatments)


tmp = load("Inequality_Constraints_to_Restrict_MTP_to_Sufficient_Statistics.rdata")
A4 = additional_inequality_constraints_part1

tmp11 = read.table("number_equality_constraints_of_first_type.txt")
a21   = tmp11$V1

tmp = load("Inequality_Constraints_to_set_monotonicity_in_hypotheses_rejected.rdata")
A5  = additional_inequality_constraints_part2

tmp = load("c.rdata")
obj = objective_function_vector

A2  = sparseMatrix(i = A2[,1],j=A2[,2],x=A2[,3])

A4  = sparseMatrix(i = A4[,1],j=A4[,2],x=A4[,3])

tmp  = dim(A2)
tmp4 = dim(A4)

if (tmp4[2]<tmp[2])
{
  A4[tmp4[1],tmp4[2]] = 0
}

tmp = dim(A1)
a1  = 0.05 * rep(1,tmp[1])

tmp = dim(A2)
a2 = c(rep(1,a21),rep(0,tmp[1]-a21))

tmp = dim(A4)
a4  = rep(0,tmp[1])

AA    = rBind(A1,A4,-A3)
power = 0.80
aa    = c(a1,a4,-power*rep(1,3))

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
A5   = sparseMatrix(i=A5[,1],j=A5[,2],x=A5[,3])
aan  = c(aa,rep(0,tmpr))
AAn  = rBind(AA,A5)

d         = dim(AAn)[2]

tmp       = cbind(seq(from=1,to=d,by=1),seq(from=1,to=d,by=1),rep(1,d))
At        = sparseMatrix(i=tmp[,1],j=tmp[,2],x=tmp[,3])

A         = rBind(AAn,At,-At,A2)
a         = c(aan,rep(1,d),rep(0,d),a2)
const.dir = c(rep("<=",dim(AAn)[1]+2*d),rep("==",dim(A2)[1]))

tmp      = Rglpk_solve_LP(obj,A,const.dir,a,max=FALSE)
solution = tmp$solution

save(solution,file="solution.RData")


